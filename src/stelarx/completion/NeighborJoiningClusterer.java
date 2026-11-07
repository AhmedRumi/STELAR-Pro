package stelarx.completion;

import stelarx.Logging;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.Arrays;

/**
 * Neighbor-Joining and BIONJ guide trees built from the existing quartet-similarity
 * matrix, as additional candidate sources alongside the UPGMA guide tree.
 *
 * <h3>Why add these at all</h3>
 * UPGMA is an <em>ultrametric</em> method: it assumes a molecular clock. The guide
 * matrix is {@code dist = 1 - sim_Q}, which is neither ultrametric nor additive, so
 * UPGMA recovers only consensus-level clades that the gene trees already witness.
 * Measured on 200-taxon duplication data, its guide tree contributes 0 new clusters
 * and 0-6 new splits to X. NJ is instead consistent for <em>additive</em> distances,
 * and {@link Transform#NEG_LOG_SIM} is the transform under which the similarity
 * becomes additive ({@code S = exp(-D)}), so these two changes together are what
 * make a different topology reachable.
 *
 * <h3>Rooting</h3>
 * Textbook NJ stops with three active nodes and emits an unrooted trifurcation. The
 * rooted DP needs a binary rooted tree, so agglomeration continues to a single
 * remaining pair: {@code n-1} merges, the last join becoming the root. The rooting
 * is therefore NJ's own join order, not an arbitrary outgroup choice.
 *
 * <h3>Complexity</h3>
 * O(n^3) time, O(n^2) memory - the pair search is O(r^2) per merge. At n = 201 that
 * is ~2.7M operations against 10^4-10^5 scored splits downstream, so the cost is not
 * measurable. Only the dense (non-packed) matrix is supported; a packed matrix means
 * n is far past the point where an O(n^3) guide tree is appropriate.
 */
public final class NeighborJoiningClusterer {

    /** Reduction formula: plain NJ, or BIONJ's variance-weighted average. */
    public enum Variant { NJ, BIONJ }

    /**
     * Similarity-to-distance transform.
     *
     * <p>{@link #ONE_MINUS_SIM} is what the UPGMA path effectively uses; it is bounded
     * but not additive. {@link #NEG_LOG_SIM} inverts {@code S = exp(-D)} and is the
     * additive choice that NJ's consistency proof actually requires.
     */
    public enum Transform {
        /** 1 - s. Bounded, not additive; the UPGMA path's implicit choice. */
        ONE_MINUS_SIM,
        /** -log s. Additive under S = exp(-D), which is what NJ's consistency needs. */
        NEG_LOG_SIM,
        /** sqrt(1 - s). Bounded like ONE_MINUS_SIM but with different curvature. */
        SQRT_ONE_MINUS_SIM,
        /**
         * -log((3s - 1) / 2). Three quartet topologies exist, so random agreement is
         * s = 1/3, not 0. Rescaling to (s - 1/3)/(1 - 1/3) before the log sends
         * random-level similarity to infinite distance instead of the finite 1.1 that
         * NEG_LOG_SIM assigns it, which is the behaviour an additive model implies.
         */
        NEG_LOG_CORRECTED,
        /** (1 - s)/s. Unbounded with a heavier tail than -log. */
        INVERSE
    }

    /** Floor applied before any division or log, so degenerate pairs stay finite. */
    private static final double SIM_FLOOR = 1e-6;

    /** Dense n*n distance matrix under the chosen transform, shared with BME search. */
    static double[] distanceMatrix(SimilarityMatrix sim, Transform transform) {
        int n = sim.n;
        double[] d = new double[n * n];
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double v = distance(sim.getSim(i, j), transform);
                d[i * n + j] = v;
                d[j * n + i] = v;
            }
        }
        return d;
    }

    /** Similarity to distance under the chosen transform. */
    private static double distance(double s, Transform t) {
        return switch (t) {
            case ONE_MINUS_SIM -> 1.0 - s;
            case SQRT_ONE_MINUS_SIM -> Math.sqrt(Math.max(0.0, 1.0 - s));
            case NEG_LOG_SIM -> -Math.log(Math.max(s, SIM_FLOOR));
            case NEG_LOG_CORRECTED ->
                -Math.log(Math.max((3.0 * s - 1.0) / 2.0, SIM_FLOOR));
            case INVERSE -> (1.0 - s) / Math.max(s, SIM_FLOOR);
        };
    }

    /**
     * Build a rooted binary guide tree over all n taxa.
     *
     * @return the tree, or {@code null} when the matrix is packed (n too large),
     *         in which case a warning is logged and the caller should skip it.
     */
    public static Tree build(SimilarityMatrix sim, int treeIndex,
                             Variant variant, Transform transform) {
        if (sim.isPacked()) {
            Logging.warn("%s/%s guide tree skipped: packed similarity matrix (n=%d) is "
                + "past the scale where an O(n^3) distance tree is appropriate",
                variant, transform, sim.n);
            return null;
        }
        int n = sim.n;
        if (n == 1) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId = 0;
            leaf.rangeStart = 0;
            leaf.rangeEnd = 1;
            return new Tree(treeIndex, leaf, new int[]{0}, new int[]{0}, 1, 1);
        }
        return new NeighborJoiningClusterer(sim, n, variant, transform).run(treeIndex);
    }

    // ── State ────────────────────────────────────────────────────────────────

    private final int n;
    private final Variant variant;
    private final Transform transform;
    private final double[] d;          // dense n*n, symmetric, active rows only meaningful
    private final double[] v;          // BIONJ variance matrix; null for plain NJ
    private final TreeNode[] clusterRoot;
    private final int[] activeArr;
    private int activeCount;
    private int zeroSimPairs;          // off-diagonal pairs with sim == 0 (never co-occur)
    private int belowRandomPairs;      // pairs at or below the 1/3 random-quartet level

    private NeighborJoiningClusterer(SimilarityMatrix sim, int n,
                                     Variant variant, Transform transform) {
        this.n = n;
        this.variant = variant;
        this.transform = transform;
        this.d = new double[n * n];
        this.v = (variant == Variant.BIONJ) ? new double[n * n] : null;
        this.clusterRoot = new TreeNode[n];
        this.activeArr = new int[n];
        this.activeCount = n;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double s = sim.getSim(i, j);
                if (s <= 0.0) zeroSimPairs++;
                if (s <= 1.0 / 3.0) belowRandomPairs++;
                double dist = distance(s, transform);
                d[i * n + j] = dist;
                d[j * n + i] = dist;
            }
            d[i * n + i] = 0.0;
        }
        if (v != null) System.arraycopy(d, 0, v, 0, d.length);

        for (int i = 0; i < n; i++) {
            activeArr[i] = i;
            TreeNode leaf = new TreeNode();
            leaf.taxonId = i;
            clusterRoot[i] = leaf;
        }
    }

    // ── Core loop ────────────────────────────────────────────────────────────

    private Tree run(int treeIndex) {
        double[] rowSum = new double[n];

        while (activeCount > 2) {
            int r = activeCount;

            // R_i = sum of d(i,k) over active k
            for (int ai = 0; ai < r; ai++) {
                int i = activeArr[ai];
                double s = 0.0;
                for (int bk = 0; bk < r; bk++) s += d[i * n + activeArr[bk]];
                rowSum[i] = s;
            }

            // Pick the pair minimising Q(i,j) = (r-2)*d(i,j) - R_i - R_j
            int bi = -1, bj = -1;
            double bestQ = Double.POSITIVE_INFINITY;
            for (int ai = 0; ai < r; ai++) {
                int i = activeArr[ai];
                for (int bk = ai + 1; bk < r; bk++) {
                    int j = activeArr[bk];
                    double q = (r - 2) * d[i * n + j] - rowSum[i] - rowSum[j];
                    if (q < bestQ) { bestQ = q; bi = i; bj = j; }
                }
            }

            double dij = d[bi * n + bj];
            // Branch lengths from bi/bj to their new parent; l_i + l_j == d_ij.
            double li = 0.5 * dij + (rowSum[bi] - rowSum[bj]) / (2.0 * (r - 2));
            double lj = dij - li;

            double lambda = 0.5;
            if (variant == Variant.BIONJ) {
                double vij = v[bi * n + bj];
                if (vij > 0.0) {
                    double acc = 0.0;
                    for (int bk = 0; bk < r; bk++) {
                        int k = activeArr[bk];
                        if (k == bi || k == bj) continue;
                        acc += v[bj * n + k] - v[bi * n + k];
                    }
                    lambda = 0.5 + acc / (2.0 * (r - 2) * vij);
                    if (lambda < 0.0) lambda = 0.0;
                    else if (lambda > 1.0) lambda = 1.0;
                }
            }

            join(bi, bj);

            // Reduce: bi becomes the new internal node u.
            for (int bk = 0; bk < activeCount; bk++) {
                int k = activeArr[bk];
                if (k == bi) continue;
                double dik = d[bi * n + k], djk = d[bj * n + k];
                double duk = (variant == Variant.BIONJ)
                    ? lambda * (dik - li) + (1.0 - lambda) * (djk - lj)
                    : 0.5 * (dik + djk - dij);
                if (duk < 0.0) duk = 0.0;   // negative branches are an artifact, not a distance
                d[bi * n + k] = duk;
                d[k * n + bi] = duk;
                if (v != null) {
                    double vuk = lambda * v[bi * n + k] + (1.0 - lambda) * v[bj * n + k]
                               - lambda * (1.0 - lambda) * v[bi * n + bj];
                    if (vuk < 0.0) vuk = 0.0;
                    v[bi * n + k] = vuk;
                    v[k * n + bi] = vuk;
                }
            }
        }

        // Two active clusters remain: join them as the root.
        int a = activeArr[0], b = activeArr[1];
        join(a, b);

        Logging.info("%s guide tree built on %s: %d of %d pairs never co-occur, "
            + "%d at or below the 1/3 random-quartet level (distance floored)",
            variant, transform, zeroSimPairs, n * (n - 1) / 2, belowRandomPairs);
        return buildTree(treeIndex);
    }

    /** Attach clusters i and j under a fresh parent; j leaves the active set. */
    private void join(int i, int j) {
        TreeNode parent = new TreeNode();
        parent.left = clusterRoot[i];
        parent.right = clusterRoot[j];
        clusterRoot[i].parent = parent;
        clusterRoot[j].parent = parent;
        clusterRoot[i] = parent;
        // Swap-and-shrink removal of j.
        for (int p = 0; p < activeCount; p++) {
            if (activeArr[p] == j) {
                activeArr[p] = activeArr[activeCount - 1];
                break;
            }
        }
        activeCount--;
    }

    // ── Tree assembly (mirrors UPGMAClusterer) ───────────────────────────────

    private Tree buildTree(int treeIndex) {
        TreeNode root = clusterRoot[activeArr[0]];
        int[] postArr = new int[n];
        int[] cursor = {0};
        assignRanges(root, postArr, cursor);
        int[] posMap = new int[n];
        Arrays.fill(posMap, -1);
        for (int pos = 0; pos < n; pos++) posMap[postArr[pos]] = pos;
        return new Tree(treeIndex, root, postArr, posMap, n, n);
    }

    private void assignRanges(TreeNode node, int[] postArr, int[] cursor) {
        if (node.isLeaf()) {
            node.rangeStart = cursor[0];
            postArr[cursor[0]++] = node.taxonId;
            node.rangeEnd = cursor[0];
        } else {
            node.rangeStart = cursor[0];
            assignRanges(node.left, postArr, cursor);
            assignRanges(node.right, postArr, cursor);
            node.rangeEnd = cursor[0];
        }
    }
}
