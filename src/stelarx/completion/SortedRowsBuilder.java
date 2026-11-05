package stelarx.completion;

import stelarx.util.Threading;

/**
 * Builds {@link SortedRows} (nearest-neighbour order per taxon) from a distance
 * source.
 *
 * Ordering of row x: ascending distance D(x, j); ties broken by descending
 * taxon ID, matching ASTRAL-MP's AbstractMatrix.sortColumn
 * ({@code comp == 0 ? -o1.compareTo(o2) : comp}).  The key (distance, −ID) is
 * unique per element, so the order is fully determined and independent of the
 * sorting algorithm.  A primitive merge sort is used instead of sorting boxed
 * Integers through a comparator: same order, no per-element allocation.
 *
 * Reference: mirrors AbstractMatrix.sortByDistance / assureOrderedTaxa in
 * astral-mp-legacy-codebase/AbstractMatrix.java.
 */
public class SortedRowsBuilder {

    /** Per-thread scratch so a 100k-row build does not churn 3·n doubles per row. */
    private static final class Scratch {
        double[] key;
        int[] idx;
        int[] tmp;

        void ensure(int n) {
            if (key == null || key.length != n) {
                key = new double[n];
                idx = new int[n];
                tmp = new int[n];
            }
        }
    }

    private static final ThreadLocal<Scratch> SCRATCH = ThreadLocal.withInitial(Scratch::new);

    /** Rows ordered by {@code sim.getDist(x, j)}; the matrix is only read. */
    public static SortedRows build(SimilarityMatrix sim) {
        int n = sim.n;
        SortedRows rows = new SortedRows(n);
        Threading.processRangeParallel(n, x -> {
            Scratch s = SCRATCH.get();
            s.ensure(n);
            for (int j = 0; j < n; j++) s.key[j] = sim.getDist(x, j);
            sortRow(s.key, s.idx, s.tmp, n);
            rows.setRow(x, s.idx);
        });
        return rows;
    }

    /**
     * Rows ordered by a flat row-major n×n distance array
     * ({@code dist[x*n + j] = D(x, j)}), e.g. {@link DistanceMatrix#dist}.
     */
    public static SortedRows buildFromDense(double[] dist, int n) {
        if ((long)dist.length != (long)n * n) {
            throw new IllegalArgumentException("distance array length " + dist.length
                + " does not match n=" + n);
        }
        SortedRows rows = new SortedRows(n);
        Threading.processRangeParallel(n, x -> {
            Scratch s = SCRATCH.get();
            s.ensure(n);
            System.arraycopy(dist, x * n, s.key, 0, n);
            sortRow(s.key, s.idx, s.tmp, n);
            rows.setRow(x, s.idx);
        });
        return rows;
    }

    /**
     * Sorts the indices 0..n−1 into {@code idx} by (key ascending, index descending).
     * Bottom-up merge sort; {@code tmp} is scratch of the same length.
     */
    static void sortRow(double[] key, int[] idx, int[] tmp, int n) {
        for (int i = 0; i < n; i++) idx[i] = i;
        int[] src = idx, dst = tmp;
        for (int width = 1; width < n; width <<= 1) {
            for (int lo = 0; lo < n; lo += width << 1) {
                int mid = Math.min(lo + width, n);
                int hi = Math.min(lo + (width << 1), n);
                int i = lo, j = mid, k = lo;
                while (i < mid && j < hi) {
                    if (before(key, src[j], src[i])) dst[k++] = src[j++];
                    else dst[k++] = src[i++];
                }
                while (i < mid) dst[k++] = src[i++];
                while (j < hi) dst[k++] = src[j++];
            }
            int[] t = src; src = dst; dst = t;
        }
        if (src != idx) System.arraycopy(src, 0, idx, 0, n);
    }

    /** True when element a sorts strictly before element b. */
    private static boolean before(double[] key, int a, int b) {
        int c = Double.compare(key[a], key[b]);
        return c != 0 ? c < 0 : a > b;   // descending ID on tie
    }
}
