import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SimilarityMatrixBuilder;
import stelarx.completion.SortedRows;
import stelarx.completion.SortedRowsBuilder;
import stelarx.completion.TreeCompleter;
import stelarx.completion.UPGMAClusterer;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.Arrays;
import java.util.List;

/**
 * Bit-exact parity of the segmented upper-triangle similarity matrix against
 * an independent dense n×n reference (the original ASTRAL-MP-style scatter
 * with mirrored writes), plus nearest-neighbour order parity against a boxed
 * comparator sort, UPGMA copy-vs-consume parity, and segment-boundary stress.
 */
public final class PackedSimilarityParityTest {
    private static void equalBits(double expected, double actual, String label) {
        if (Double.doubleToRawLongBits(expected) != Double.doubleToRawLongBits(actual)) {
            throw new AssertionError(label + ": expected=" + expected + " actual=" + actual);
        }
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass one incomplete tree file");
        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        try {
            System.clearProperty("stelarpro.similarity.segmentShift");
            run(args[0], "default segments");
            // Four doubles per segment: every row of a 6+-taxon triangle straddles
            // segment boundaries, exercising the segmented index arithmetic.
            System.setProperty("stelarpro.similarity.segmentShift", "2");
            run(args[0], "4-double segments");
        } finally {
            System.clearProperty("stelarpro.similarity.segmentShift");
            Threading.shutdown();
        }
    }

    private static void run(String path, String label) throws Exception {
        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(path, registry);
        int n = registry.size();

        // ── 1. Similarity values: segmented triangle vs dense reference ──────
        DenseReference dense = new DenseReference(trees, n);
        SimilarityMatrix packed = SimilarityMatrixBuilder.buildCPU(trees, n);
        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                equalBits(dense.sim[a * n + b], packed.getSim(a, b), label + " similarity " + a + "," + b);
                equalBits(dense.dist[a * n + b], packed.getDist(a, b), label + " distance " + a + "," + b);
            }
        }

        // ── 2. Nearest-neighbour order vs boxed comparator sort (old builder) ─
        int[] referenceRows = referenceSortedRows(dense.dist, n);
        SortedRows rows = SortedRowsBuilder.build(packed);
        SortedRows denseRows = SortedRowsBuilder.buildFromDense(dense.dist, n);
        for (int a = 0; a < n; a++) {
            for (int rank = 0; rank < n; rank++) {
                check(referenceRows[a * n + rank] == rows.get(a, rank),
                    label + " sorted-row mismatch at row=" + a + " rank=" + rank);
                check(referenceRows[a * n + rank] == denseRows.get(a, rank),
                    label + " dense sorted-row mismatch at row=" + a + " rank=" + rank);
            }
        }

        // ── 3. Completion through rows equals completion through the matrix ──
        List<Tree> viaRows = TreeCompleter.completeAll(trees, rows, n);
        List<Tree> viaMatrix = TreeCompleter.completeAll(trees, packed, n);
        for (int i = 0; i < trees.size(); i++) {
            String expected = viaRows.get(i).toNewick(registry);
            String actual = viaMatrix.get(i).toNewick(registry);
            check(expected.equals(actual), label + " completion mismatch for tree " + i
                + ":\n" + expected + "\n" + actual);
            check(viaRows.get(i).isComplete && viaRows.get(i).leafCount == n,
                label + " tree " + i + " not completed");
        }

        // ── 4. UPGMA on a copy equals UPGMA consuming the triangle in place ──
        Tree copied = UPGMAClusterer.build(packed, trees.size());
        check(!packed.isConsumed(), label + " copy-mode UPGMA consumed the matrix");
        equalBits(dense.sim[0 * n + (n - 1)], packed.getSim(0, n - 1), label + " matrix changed by copy-mode UPGMA");
        Tree consumed = UPGMAClusterer.build(packed, trees.size(), true);
        check(packed.isConsumed(), label + " consume-mode UPGMA left the matrix readable");
        String copiedNewick = copied.toNewick(registry);
        String consumedNewick = consumed.toNewick(registry);
        check(copiedNewick.equals(consumedNewick),
            label + " UPGMA mismatch:\n" + copiedNewick + "\n" + consumedNewick);
        try {
            packed.getSim(0, 1);
            throw new AssertionError(label + " consumed matrix still answers getSim");
        } catch (IllegalStateException expected) {
            // a stale reader must fail loudly, never observe UPGMA's mutations
        }
        try {
            UPGMAClusterer.build(packed, trees.size());
            throw new AssertionError(label + " consumed matrix was clustered again");
        } catch (IllegalStateException expected) {
            // same guard for a second clustering
        }
        System.out.println("Packed similarity/rows/UPGMA/completion parity (" + label + "): PASS");
    }

    /** Old SortedRowsBuilder.buildCPU: boxed Integer sort with the ASTRAL-MP comparator. */
    private static int[] referenceSortedRows(double[] dist, int n) {
        int[] sortedRows = new int[n * n];
        for (int x = 0; x < n; x++) {
            Integer[] indices = new Integer[n];
            for (int j = 0; j < n; j++) indices[j] = j;
            final int base = x * n;
            Arrays.sort(indices, (a, b) -> {
                int c = Double.compare(dist[base + a], dist[base + b]);
                return c != 0 ? c : Integer.compare(b, a);
            });
            for (int rank = 0; rank < n; rank++) sortedRows[base + rank] = indices[rank];
        }
        return sortedRows;
    }

    /**
     * The original dense CPU accumulation (flat n×n numerator/denominator with
     * mirrored writes), kept here verbatim as the independent oracle.
     */
    private static final class DenseReference {
        final double[] sim;
        final double[] dist;

        DenseReference(List<Tree> trees, int n) {
            double[] numSum = new double[n * n];
            double[] denSum = new double[n * n];
            for (Tree tree : trees) accumulate(tree, n, numSum, denSum);
            sim = new double[n * n];
            dist = new double[n * n];
            for (int i = 0; i < sim.length; i++) {
                sim[i] = (denSum[i] > 0.0) ? numSum[i] / denSum[i] : 0.0;
                dist[i] = 1.0 - sim[i];
            }
            for (int a = 0; a < n; a++) {
                sim[a * n + a] = 1.0;
                dist[a * n + a] = 0.0;
            }
        }

        private static long c2(int x) { return (x < 2) ? 0L : (long)x * (x - 1) / 2; }

        private static void accumulate(Tree tree, int n, double[] numSum, double[] denSum) {
            int kt = tree.leafCount;
            if (kt < 4) return;
            long denPerPair = c2(kt - 2);
            int[] postArr = tree.postorderArray;
            scatterAtNode(tree.root, tree, kt, n, numSum);
            for (int i = 0; i < kt; i++) {
                int a = postArr[i];
                long rowOff = (long)a * n;
                for (int j = 0; j < kt; j++) {
                    if (i == j) continue;
                    denSum[(int)(rowOff + postArr[j])] += denPerPair;
                }
            }
        }

        private static void scatterAtNode(TreeNode node, Tree tree, int kt, int n, double[] numSum) {
            if (node.isLeaf()) return;
            scatterAtNode(node.left, tree, kt, n, numSum);
            scatterAtNode(node.right, tree, kt, n, numSum);
            int subL = node.left.rangeEnd - node.left.rangeStart;
            int subR = node.right.rangeEnd - node.right.rangeStart;
            int subU = node.rangeEnd - node.rangeStart;
            int subO = kt - subU;
            long cL = c2(subL), cR = c2(subR), cO = c2(subO);
            long totalPairs = cL + cR + cO;
            int[] postArr = tree.postorderArray;
            long simLR = totalPairs - cL - cR;
            if (simLR != 0) {
                scatterRangeRange(postArr, node.left.rangeStart, node.left.rangeEnd,
                    node.right.rangeStart, node.right.rangeEnd, simLR, numSum, n);
            }
            if (subO > 0) {
                long simLO = totalPairs - cL - cO;
                if (simLO != 0) {
                    scatterRangeOthers(postArr, kt, node.left.rangeStart, node.left.rangeEnd,
                        node.rangeStart, node.rangeEnd, simLO, numSum, n);
                }
                long simRO = totalPairs - cR - cO;
                if (simRO != 0) {
                    scatterRangeOthers(postArr, kt, node.right.rangeStart, node.right.rangeEnd,
                        node.rangeStart, node.rangeEnd, simRO, numSum, n);
                }
            }
        }

        private static void scatterRangeRange(int[] postArr, int aLo, int aHi, int bLo, int bHi,
                                              long sim, double[] numSum, int n) {
            double sd = (double)sim;
            for (int pi = aLo; pi < aHi; pi++) {
                int a = postArr[pi];
                for (int pj = bLo; pj < bHi; pj++) {
                    int b = postArr[pj];
                    numSum[a * n + b] += sd;
                    numSum[b * n + a] += sd;
                }
            }
        }

        private static void scatterRangeOthers(int[] postArr, int kt, int aLo, int aHi,
                                               int subLo, int subHi, long sim,
                                               double[] numSum, int n) {
            double sd = (double)sim;
            for (int pi = aLo; pi < aHi; pi++) {
                int a = postArr[pi];
                for (int pj = 0; pj < subLo; pj++) {
                    int b = postArr[pj];
                    numSum[a * n + b] += sd;
                    numSum[b * n + a] += sd;
                }
                for (int pj = subHi; pj < kt; pj++) {
                    int b = postArr[pj];
                    numSum[a * n + b] += sd;
                    numSum[b * n + a] += sd;
                }
            }
        }
    }
}
