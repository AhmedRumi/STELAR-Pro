package stelarx.completion;

/**
 * Taxon similarity matrix for quartet (S2/S3) or rooted-triplet (S4) guides.
 *
 * For taxa a and b:
 *   M[a][b] = Σ_{t: a,b present} num_t(a,b)
 *             ──────────────────────────────
 *             Σ_{t: a,b present} den_t(a,b)
 *
 * The builder supplies the quartet or triplet numerator and denominator;
 * this container normalizes either measure in the same way.
 *
 * Storage.  The matrix is symmetric, so every pair is stored exactly once in
 * the upper triangle (a ≤ b), held in {@link SegmentedDoubleArray} segments so
 * that no single Java array ever has to address n² cells.  During accumulation
 * two triangles exist (numerator and denominator, 8·n(n+1)/2 bytes in total);
 * {@link #normalize()} divides in place and releases the denominator, leaving
 * one triangle (4·n(n+1)/2 ≈ 4·n² bytes) that holds the similarity.
 *
 * After normalize():
 *   sim[a][b]  ∈ [0,1], diagonal = 1, symmetric
 *   dist[a][b] = 1 − sim[a][b]   (used to rank completion anchors)
 *
 * Pairs that never co-occur in any tree have sim = 0, dist = 1.
 *
 * Lifetime.  UPGMA needs a mutable copy of the triangle; when nothing reads the
 * matrix afterwards, {@link UPGMAClusterer#build(SimilarityMatrix, int, boolean)}
 * can consume the triangle in place instead of copying it.  A consumed matrix
 * rejects every further read with {@link IllegalStateException}.
 */
public class SimilarityMatrix {

    public final int n;

    /** Σ num_t(a,b) per pair while accumulating; the similarity after normalize(). */
    private SegmentedDoubleArray num;

    /** Σ den_t(a,b) per pair; released by normalize(). */
    private SegmentedDoubleArray den;

    private boolean normalized;
    /** True only when similarity preprocessing crossed the wide Java-array boundary. */
    private boolean streamedHostBatches;

    public SimilarityMatrix(int n) {
        if (n < 0) throw new IllegalArgumentException("negative taxon count: " + n);
        this.n = n;
        long cells = triangleCellCount(n);
        this.num = new SegmentedDoubleArray(cells);
        this.den = new SegmentedDoubleArray(cells);
    }

    /**
     * Finalize: sim[a][b] = numSum / denSum in place, then release the denominator.
     * Pairs that never co-occur get sim = 0 (dist = 1).  Diagonal is set to sim = 1.
     */
    public void normalize() {
        requireLive();
        if (normalized) throw new IllegalStateException("similarity matrix already normalized");
        double[][] nums = num.segments();
        double[][] dens = den.segments();
        for (int s = 0; s < nums.length; s++) {
            double[] numSeg = nums[s];
            double[] denSeg = dens[s];
            for (int i = 0; i < numSeg.length; i++) {
                numSeg[i] = denSeg[i] > 0.0 ? numSeg[i] / denSeg[i] : 0.0;
            }
        }
        for (int a = 0; a < n; a++) num.set(index(a, a), 1.0);
        // Make the (potentially multi-GiB) denominator immediately collectible.
        den = null;
        normalized = true;
    }

    /** Returns M[a][b]. Call after normalize(). */
    public double getSim(int a, int b) {
        return similarity().get(index(a, b));
    }

    /** Returns (1 − M[a][b]). Call after normalize(). */
    public double getDist(int a, int b) {
        return 1.0 - getSim(a, b);
    }

    public boolean isNormalized() { return normalized; }

    /** True once the similarity triangle has been handed to an in-place consumer. */
    public boolean isConsumed() { return num == null; }

    void markStreamedHostBatches() { streamedHostBatches = true; }
    public boolean usedStreamedHostBatches() { return streamedHostBatches; }

    // ── Storage arithmetic ───────────────────────────────────────────────────

    static long triangleCellCount(int n) {
        return (long)n * (n + 1L) / 2L;
    }

    /** Heap bytes held by both accumulators during a build of an n-taxon matrix. */
    public static long accumulatorBytes(int n) {
        return 2L * triangleCellCount(n) * Double.BYTES;
    }

    /** Heap bytes held by the normalized similarity triangle. */
    public static long similarityBytes(int n) {
        return triangleCellCount(n) * Double.BYTES;
    }

    static long packedIndex(int n, int a, int b) {
        if (a > b) { int t = a; a = b; b = t; }
        return (long)a * n - (long)a * (a + 1L) / 2L + b;
    }

    long index(int a, int b) {
        return packedIndex(n, a, b);
    }

    // ── Accumulation (builder-only) ──────────────────────────────────────────

    void addNumerator(int a, int b, double value) {
        num.add(index(a, b), value);
    }

    void addDenominator(int a, int b, double value) {
        den.add(index(a, b), value);
    }

    double[][] numeratorSegments() { requireAccumulating(); return num.segments(); }
    double[][] denominatorSegments() { requireAccumulating(); return den.segments(); }
    int segmentShift() { return num.segmentShift(); }

    // ── Triangle hand-over ───────────────────────────────────────────────────

    /** Normalized triangle for read-only use. */
    private SegmentedDoubleArray similarity() {
        requireLive();
        if (!normalized) throw new IllegalStateException("similarity matrix is not normalized");
        return num;
    }

    /** Fresh mutable copy of the normalized triangle; this matrix stays readable. */
    SegmentedDoubleArray copySimilarity() {
        SegmentedDoubleArray src = similarity();
        SegmentedDoubleArray copy = new SegmentedDoubleArray(src.length(), src.segmentShift());
        double[][] from = src.segments();
        double[][] to = copy.segments();
        for (int s = 0; s < from.length; s++) {
            System.arraycopy(from[s], 0, to[s], 0, from[s].length);
        }
        return copy;
    }

    /**
     * Hands the normalized triangle to a caller that will mutate it.  No copy is
     * made; this matrix becomes unreadable so a stale reference cannot observe
     * the mutation silently.
     */
    SegmentedDoubleArray consumeSimilarity() {
        SegmentedDoubleArray triangle = similarity();
        num = null;
        return triangle;
    }

    private void requireLive() {
        if (num == null) {
            throw new IllegalStateException("similarity matrix was consumed by its last reader");
        }
    }

    private void requireAccumulating() {
        requireLive();
        if (normalized) throw new IllegalStateException("similarity matrix already normalized");
    }
}
