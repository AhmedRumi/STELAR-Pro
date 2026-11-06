package stelarx.completion;

/**
 * Per-taxon nearest-neighbour order derived from a distance matrix.
 *
 * {@code get(x, rank)} is the taxon ID of x's rank-th nearest neighbour
 * (rank 0 is x itself at distance 0 unless another taxon ties at distance 0,
 * in which case the usual descending-ID tie rule applies).  This is the only
 * matrix-derived information tree completion needs: the first ranked taxon
 * that is present in a tree is the insertion anchor.
 *
 * Rows are stored as one array per taxon so no single Java array has to hold
 * n² entries, and as unsigned 16-bit values whenever every taxon ID fits
 * (n ≤ 65,536), halving the footprint to 2·n² bytes.
 */
public final class SortedRows {
    public final int n;
    private final char[][] narrow;
    private final int[][] wide;

    SortedRows(int n) {
        this.n = n;
        if (n <= Character.MAX_VALUE + 1) {
            this.narrow = new char[n][];
            this.wide = null;
        } else {
            this.narrow = null;
            this.wide = new int[n][];
        }
    }

    /** Taxon ID at the given rank in x's nearest-neighbour order. */
    public int get(int x, int rank) {
        return narrow != null ? narrow[x][rank] : wide[x][rank];
    }

    void setRow(int x, int[] order) {
        if (order.length != n) throw new IllegalArgumentException("row length " + order.length);
        if (narrow != null) {
            char[] row = new char[n];
            for (int rank = 0; rank < n; rank++) row[rank] = (char)order[rank];
            narrow[x] = row;
        } else {
            wide[x] = order.clone();
        }
    }

    /** Heap bytes used by the row payload. */
    public long bytes() {
        return (long)n * n * (narrow != null ? Character.BYTES : Integer.BYTES);
    }
}
