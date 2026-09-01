package stelarx.tree;

/**
 * Compact per-taxon lists of leaf positions for one gene tree.
 *
 * <p>The CSR row for taxon {@code v} is
 * {@code positions[offsets[v]..offsets[v + 1])}. Positions are sorted, so
 * membership in a subtree range is answered by binary search. Unlike the
 * legacy one-position map, this representation retains every gene copy.</p>
 */
public final class TaxonPositionIndex {
    private final int[] offsets;
    private final int[] positions;
    private final int distinctTaxonCount;

    private TaxonPositionIndex(int[] offsets, int[] positions,
                               int distinctTaxonCount) {
        this.offsets = offsets;
        this.positions = positions;
        this.distinctTaxonCount = distinctTaxonCount;
    }

    /** Build sorted position vectors in O(n + L) time for n taxa and L leaves. */
    public static TaxonPositionIndex build(int[] postorderArray, int totalTaxa) {
        int[] offsets = new int[totalTaxa + 1];
        for (int taxon : postorderArray) {
            if (taxon < 0 || taxon >= totalTaxa) {
                throw new IllegalArgumentException("Taxon ID outside registry: " + taxon);
            }
            offsets[taxon + 1]++;
        }

        int distinct = 0;
        for (int taxon = 0; taxon < totalTaxa; taxon++) {
            if (offsets[taxon + 1] != 0) distinct++;
            offsets[taxon + 1] += offsets[taxon];
        }

        int[] positions = new int[postorderArray.length];
        int[] cursors = offsets.clone();
        // Scanning positions in ascending order makes every taxon row sorted.
        for (int pos = 0; pos < postorderArray.length; pos++) {
            int taxon = postorderArray[pos];
            positions[cursors[taxon]++] = pos;
        }
        return new TaxonPositionIndex(offsets, positions, distinct);
    }

    public int distinctTaxonCount() { return distinctTaxonCount; }
    public int positionCount() { return positions.length; }
    public int taxonCount() { return offsets.length - 1; }

    /** CSR start offset for one taxon; used when flattening indexes for CUDA. */
    public int startOffset(int taxon) { return offsets[taxon]; }

    /** CSR end offset for one taxon; used when flattening indexes for CUDA. */
    public int endOffset(int taxon) { return offsets[taxon + 1]; }

    /** Copy all sorted leaf positions into a packed CUDA host array. */
    public void copyPositionsTo(int[] destination, int destinationOffset) {
        System.arraycopy(positions, 0, destination, destinationOffset, positions.length);
    }

    /** Return whether at least one copy of {@code taxon} occurs in [lo, hi). */
    public boolean containsInRange(int taxon, int lo, int hi) {
        int rowStart = offsets[taxon];
        int rowEnd = offsets[taxon + 1];
        int found = lowerBound(positions, rowStart, rowEnd, lo);
        return found < rowEnd && positions[found] < hi;
    }

    /**
     * Return whether {@code position} is the first copy of its taxon in [lo, hi).
     * This suppresses duplicate copies while scanning a subtree's leaf range.
     */
    public boolean isFirstInRange(int taxon, int position, int lo) {
        int rowStart = offsets[taxon];
        int rowEnd = offsets[taxon + 1];
        int first = lowerBound(positions, rowStart, rowEnd, lo);
        return first < rowEnd && positions[first] == position;
    }

    /** Count copies of one taxon in [lo, hi); primarily useful for verification. */
    public int countInRange(int taxon, int lo, int hi) {
        int rowStart = offsets[taxon];
        int rowEnd = offsets[taxon + 1];
        return lowerBound(positions, rowStart, rowEnd, hi)
             - lowerBound(positions, rowStart, rowEnd, lo);
    }

    private static int lowerBound(int[] values, int from, int to, int key) {
        int lo = from;
        int hi = to;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (values[mid] < key) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
}
