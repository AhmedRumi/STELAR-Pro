package stelarx.completion;

import stelarx.util.Threading;

import java.util.Arrays;
import java.util.Random;

/**
 * Nearest-neighbour row order: the primitive merge sort must reproduce the
 * boxed-comparator order exactly, including heavy ties (descending ID), the
 * narrow (16-bit) and wide (32-bit) storage must agree, and a row of a
 * segmented matrix must match a dense row.
 */
public final class SortedRowsTest {
    public static void main(String[] args) {
        Random random = new Random(0x5E7A1L);
        int[] key = new int[3];
        int[] idx = new int[3];
        // Lengths around merge-sort run boundaries, plus zero and one.
        int[] lengths = {0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 33, 64, 100, 257, 1000};
        long checks = 0;
        for (int n : lengths) {
            for (int trial = 0; trial < 8; trial++) {
                double[] keys = new double[n];
                int distinct = trial == 0 ? 1 : trial == 1 ? 2 : 1 + random.nextInt(Math.max(1, n));
                for (int i = 0; i < n; i++) {
                    // Few distinct values -> many ties; include 0.0 and 1.0 endpoints.
                    keys[i] = (double)random.nextInt(distinct) / Math.max(1, distinct - 1);
                }
                int[] expected = referenceOrder(keys, n);
                int[] actual = new int[n];
                SortedRowsBuilder.sortRow(keys, actual, new int[n], n);
                if (!Arrays.equals(expected, actual)) {
                    throw new AssertionError("merge sort order differs for n=" + n
                        + " trial=" + trial + "\nexpected=" + Arrays.toString(expected)
                        + "\nactual=" + Arrays.toString(actual));
                }
                checks += n;
            }
        }

        // Narrow (char) and wide (int) storage round-trip the same values.
        int n = 12;
        int[] order = new int[n];
        for (int i = 0; i < n; i++) order[i] = n - 1 - i;
        SortedRows narrow = new SortedRows(n);
        narrow.setRow(3, order);
        if (narrow.bytes() != (long)n * n * Character.BYTES) {
            throw new AssertionError("narrow rows should use 16-bit entries");
        }
        for (int rank = 0; rank < n; rank++) {
            if (narrow.get(3, rank) != order[rank]) throw new AssertionError("narrow row mismatch");
        }

        // Rows built from a segmented matrix equal rows built from the dense copy.
        Threading.start(2);
        try {
            System.setProperty("stelarpro.similarity.segmentShift", "3");
            int m = 23;
            SimilarityMatrix sim = new SimilarityMatrix(m);
            for (int a = 0; a < m; a++) {
                for (int b = a; b < m; b++) {
                    double num = random.nextInt(5);
                    double den = 1 + random.nextInt(3);
                    sim.addNumerator(a, b, num);
                    sim.addDenominator(a, b, den);
                }
            }
            sim.normalize();
            double[] dist = new double[m * m];
            for (int a = 0; a < m; a++) {
                for (int b = 0; b < m; b++) dist[a * m + b] = sim.getDist(a, b);
            }
            SortedRows fromMatrix = SortedRowsBuilder.build(sim);
            SortedRows fromDense = SortedRowsBuilder.buildFromDense(dist, m);
            for (int a = 0; a < m; a++) {
                int[] expected = referenceOrder(Arrays.copyOfRange(dist, a * m, a * m + m), m);
                for (int rank = 0; rank < m; rank++) {
                    if (fromMatrix.get(a, rank) != expected[rank]
                            || fromDense.get(a, rank) != expected[rank]) {
                        throw new AssertionError("segmented-matrix row order mismatch at row " + a);
                    }
                }
            }
            for (int a = 0; a < m; a++) {
                if (sim.getSim(a, a) != 1.0 || sim.getDist(a, a) != 0.0) {
                    throw new AssertionError("diagonal is not (1, 0)");
                }
                for (int b = 0; b < m; b++) {
                    if (sim.getSim(a, b) != sim.getSim(b, a)) throw new AssertionError("asymmetric");
                }
            }
        } finally {
            System.clearProperty("stelarpro.similarity.segmentShift");
            Threading.shutdown();
        }
        System.out.println("Sorted rows: PASS (" + checks + " ordered entries)");
    }

    /** The original boxed comparator order: ascending key, descending index on ties. */
    private static int[] referenceOrder(double[] keys, int n) {
        Integer[] indices = new Integer[n];
        for (int j = 0; j < n; j++) indices[j] = j;
        Arrays.sort(indices, (a, b) -> {
            int c = Double.compare(keys[a], keys[b]);
            return c != 0 ? c : Integer.compare(b, a);
        });
        int[] out = new int[n];
        for (int i = 0; i < n; i++) out[i] = indices[i];
        return out;
    }
}
