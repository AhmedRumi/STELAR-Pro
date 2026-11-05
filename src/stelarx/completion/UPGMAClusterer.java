package stelarx.completion;

import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.Threading;

import java.util.Arrays;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;

/**
 * Efficient UPGMA clusterer using flat primitive arrays + pointer-based tree construction.
 *
 * <h3>Algorithm summary</h3>
 * Standard UPGMA agglomerative clustering (n−1 merges):
 * <ol>
 *   <li>Find the active pair (I, J) with highest similarity.</li>
 *   <li>Merge J into I: build a new TreeNode parent, update weights.</li>
 *   <li>Update sim(I, k) for every other active k via weighted average.</li>
 *   <li>Repair per-row best-neighbor caches; recompute stale rows.</li>
 * </ol>
 * After n−1 merges the single remaining cluster root is the UPGMA dendrogram.
 * A post-order walk assigns rangeStart/rangeEnd/postorderArray and returns a Tree.
 *
 * <h3>Data structures</h3>
 * <ul>
 *   <li>{@code SegmentedDoubleArray mat} — upper triangle (i ≤ j), n(n+1)/2 doubles.
 *       Either a private copy of the similarity triangle or, when the caller no
 *       longer needs the matrix, the matrix's own triangle consumed in place
 *       (no second n²/2 allocation).</li>
 *   <li>{@code double[] bestSim}, {@code int[] bestJ} — per-row best-neighbor cache,
 *       O(n) total.</li>
 *   <li>{@code int[] activeArr} + {@code int[] posInActive} — swap-and-shrink active set:
 *       removing a deactivated cluster J from activeArr takes O(1).</li>
 *   <li>{@code TreeNode[] clusterRoot} — pointer to current cluster's subtree root.
 *       Bipartitions are extracted after the fact via rangeStart/rangeEnd.</li>
 * </ul>
 *
 * <h3>Parallelization (intra-iteration)</h3>
 * <ul>
 *   <li><b>Find-best</b>: parallel reduction over {@code bestSim[]}.</li>
 *   <li><b>Row update</b>: parallel-for over active k ≠ I.  Each k is owned by exactly
 *       one thread; writes to {@code mat[idx(I,k)]}, {@code bestSim[k]}, {@code bestJ[k]}
 *       are disjoint across threads.</li>
 *   <li><b>Stale-row recompute</b>: stale k indices collected in a
 *       {@code ConcurrentLinkedQueue} during the parallel update, then dispatched as a
 *       second parallel batch.</li>
 *   <li><b>Recompute row I</b>: always serial (I's entire row changed).</li>
 * </ul>
 *
 * <h3>Complexity</h3>
 * <ul>
 *   <li>Time: O(n²) expected; O(n³) worst case (every row stale every iteration).</li>
 *   <li>Memory: O(n²) for the upper triangle; O(n) for all other structures.</li>
 * </ul>
 */
public class UPGMAClusterer {

    // ── Instance state ────────────────────────────────────────────────────────

    private final int n;
    private final SegmentedDoubleArray mat;   // upper triangle, n*(n+1)/2
    private final double[] bestSim;           // per-row best similarity
    private final int[] bestJ;                // per-row best neighbor index
    private final int[] weight;
    private final TreeNode[] clusterRoot;     // current subtree root per cluster

    // Active-set maintained via swap-and-shrink (O(1) removal)
    private final int[] activeArr;            // activeArr[0..activeCount-1] = active indices
    private final int[] posInActive;          // posInActive[i] = position of i in activeArr
    private int activeCount;

    // ── Construction ─────────────────────────────────────────────────────────

    private UPGMAClusterer(SegmentedDoubleArray triangle, int n) {
        this.n = n;
        this.mat = triangle;
        this.bestSim = new double[n];
        this.bestJ = new int[n];
        this.weight = new int[n];
        this.clusterRoot = new TreeNode[n];
        this.activeArr = new int[n];
        this.posInActive = new int[n];
        this.activeCount = n;

        // Initialise active set (identity ordering) and leaf TreeNodes
        for (int i = 0; i < n; i++) {
            activeArr[i] = i;
            posInActive[i] = i;
            weight[i] = 1;
            TreeNode leaf = new TreeNode();
            leaf.taxonId = i;
            clusterRoot[i] = leaf;
        }

        // Build per-row best-neighbor cache (i-then-j scan, strict comparison)
        for (int i = 0; i < n; i++) {
            double bs = Double.NEGATIVE_INFINITY;
            int bj = -1;
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                double s = get(i, j);
                if (s > bs) { bs = s; bj = j; }
            }
            bestSim[i] = bs;
            bestJ[i] = bj;
        }
    }

    // ── Public entry points ───────────────────────────────────────────────────

    /**
     * Run UPGMA on a private copy of the similarity triangle; {@code sim} stays
     * readable afterwards.
     */
    public static Tree build(SimilarityMatrix sim, int treeIndex) {
        return build(sim, treeIndex, false);
    }

    /**
     * Run UPGMA on {@code sim} and return the dendrogram as a {@link Tree} with
     * all n taxa.
     *
     * @param consume when true the matrix's own triangle is mutated in place and
     *                {@code sim} becomes unreadable (saves one n²/2 double
     *                allocation); when false a copy is clustered.
     */
    public static Tree build(SimilarityMatrix sim, int treeIndex, boolean consume) {
        int n = sim.n;
        if (n == 1) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId = 0;
            leaf.rangeStart = 0;
            leaf.rangeEnd = 1;
            return new Tree(treeIndex, leaf, new int[]{0}, new int[]{0}, 1, 1);
        }
        SegmentedDoubleArray triangle = consume ? sim.consumeSimilarity() : sim.copySimilarity();
        return new UPGMAClusterer(triangle, n).run(treeIndex);
    }

    // ── Core UPGMA loop ───────────────────────────────────────────────────────

    private Tree run(int treeIndex) {
        int threads = Threading.getNumThreads();
        for (int iter = 0; iter < n - 1; iter++) {

            // ── Step 1: Find global best I (parallel reduction) ───────────────
            int I = findBestI(threads);
            int J = bestJ[I];

            // ── Step 2: Build new internal TreeNode (merge J into I) ──────────
            TreeNode newNode = new TreeNode();
            newNode.left = clusterRoot[I];
            newNode.right = clusterRoot[J];
            clusterRoot[I].parent = newNode;
            clusterRoot[J].parent = newNode;
            clusterRoot[I] = newNode;

            final double wI = weight[I];
            final double wJ = weight[J];
            weight[I] += weight[J];

            // ── Step 3: Deactivate J in O(1) via swap-and-shrink ─────────────
            removeActive(J);

            // ── Step 4: Parallel update of sim(I, k) and row-k caches ─────────
            final int snapCount = activeCount;
            final int finalI = I;
            final int finalJ = J;
            ConcurrentLinkedQueue<Integer> staleQueue = new ConcurrentLinkedQueue<>();
            Threading.processRangeParallel(snapCount, ai -> {
                int k = activeArr[ai];
                if (k == finalI) return;

                double oldIK = get(finalI, k);
                double oldJK = get(finalJ, k);
                double newIK = (oldIK * wI + oldJK * wJ) / (wI + wJ);
                set(finalI, k, newIK);

                if (bestJ[k] == finalJ) {
                    // J gone: forced full recompute
                    staleQueue.add(k);
                } else if (newIK > bestSim[k]) {
                    // I is now the best neighbor of k
                    bestJ[k] = finalI;
                    bestSim[k] = newIK;
                } else if (bestJ[k] == finalI && newIK < bestSim[k]) {
                    // cached best got worse; something else may now be better
                    staleQueue.add(k);
                }
                // else: bestJ[k] is some active m ≠ I,J; sim(k,m) unchanged; valid
            });

            // ── Step 5: Parallel recompute of stale rows ──────────────────────
            Integer[] staleArr = staleQueue.toArray(new Integer[0]);
            if (staleArr.length > 0) {
                Threading.processRangeParallel(staleArr.length, si -> recomputeBest(staleArr[si]));
            }

            // ── Step 6: Recompute row I (serial — its entire row changed) ─────
            recomputeBest(I);
        }
        return buildTree(treeIndex);
    }

    // ── Parallel find-best ────────────────────────────────────────────────────

    /**
     * Parallel reduction: find active i with maximum bestSim[i].
     * T threads each scan a disjoint chunk; one serial reduce over T local winners.
     */
    private int findBestI(int threads) {
        final int snap = activeCount;
        int actual = Math.min(threads, snap);
        int chunk = Math.max(1, (snap + actual - 1) / actual);
        int[] localI = new int[actual];
        double[] localS = new double[actual];
        Arrays.fill(localI, -1);
        Arrays.fill(localS, Double.NEGATIVE_INFINITY);

        CountDownLatch latch = new CountDownLatch(actual);
        for (int t = 0; t < actual; t++) {
            int tid = t;
            int lo = t * chunk;
            int hi = Math.min(lo + chunk, snap);
            Threading.submit(() -> {
                int bI = -1;
                double bS = Double.NEGATIVE_INFINITY;
                for (int ai = lo; ai < hi; ai++) {
                    int i = activeArr[ai];
                    if (bestSim[i] > bS) { bS = bestSim[i]; bI = i; }
                }
                localI[tid] = bI;
                localS[tid] = bS;
                latch.countDown();
            });
        }
        try { latch.await(); }
        catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException(e);
        }

        int I = -1;
        double best = Double.NEGATIVE_INFINITY;
        for (int t = 0; t < actual; t++) {
            if (localS[t] > best) { best = localS[t]; I = localI[t]; }
        }
        return I;
    }

    // ── Active-set management (O(1) removal via swap-and-shrink) ─────────────

    private void removeActive(int j) {
        int pos = posInActive[j];
        int last = activeArr[activeCount - 1];
        activeArr[pos] = last;
        posInActive[last] = pos;
        activeCount--;
    }

    // ── Row best-neighbor recompute ───────────────────────────────────────────

    /** Full linear scan of active row k to recompute bestSim[k] / bestJ[k]. */
    private void recomputeBest(int k) {
        double bs = Double.NEGATIVE_INFINITY;
        int bj = -1;
        for (int ai = 0; ai < activeCount; ai++) {
            int m = activeArr[ai];
            if (m == k) continue;
            double s = get(k, m);
            if (s > bs) { bs = s; bj = m; }
        }
        bestSim[k] = bs;
        bestJ[k] = bj;
    }

    // ── Upper-triangle access ─────────────────────────────────────────────────

    private double get(int i, int j) { return mat.get(SimilarityMatrix.packedIndex(n, i, j)); }
    private void set(int i, int j, double value) { mat.set(SimilarityMatrix.packedIndex(n, i, j), value); }

    // ── Tree construction ─────────────────────────────────────────────────────

    /**
     * After all merges, the single remaining active cluster's root is the UPGMA
     * dendrogram root.  Assigns rangeStart/rangeEnd + postorderArray via DFS,
     * then wraps everything in a {@link Tree}.
     */
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
