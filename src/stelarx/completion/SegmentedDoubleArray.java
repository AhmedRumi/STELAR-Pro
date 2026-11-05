package stelarx.completion;

/**
 * A long-indexed double array split into fixed-size Java arrays.
 *
 * Java primitive arrays are int-indexed and therefore cannot represent the
 * large symmetric matrices needed above roughly 65k taxa.  Segmentation keeps
 * every physical allocation comfortably below that limit while retaining
 * exact double precision and O(1) access.
 */
final class SegmentedDoubleArray {
    static final int SEGMENT_SHIFT = 26;             // 2^26 doubles = 512 MiB
    private final long length;
    private final int segmentShift;
    private final int segmentSize;
    private final long segmentMask;
    private final double[][] segments;

    /**
     * Segment size used by production matrices.  The system property
     * {@code stelarpro.similarity.segmentShift} (1..30) exists so tests can force
     * many small segments and exercise segment-boundary arithmetic on tiny inputs.
     */
    static int defaultSegmentShift() {
        int shift = Integer.getInteger("stelarpro.similarity.segmentShift", SEGMENT_SHIFT);
        if (shift < 1 || shift > 30) {
            throw new IllegalArgumentException("stelarpro.similarity.segmentShift must be 1..30: " + shift);
        }
        return shift;
    }

    SegmentedDoubleArray(long length) {
        this(length, defaultSegmentShift());
    }

    SegmentedDoubleArray(long length, int segmentShift) {
        if (length < 0) throw new IllegalArgumentException("negative segmented-array length");
        if (segmentShift < 1 || segmentShift > 30) {
            throw new IllegalArgumentException("invalid segment shift: " + segmentShift);
        }
        this.segmentShift = segmentShift;
        this.segmentSize = 1 << segmentShift;
        this.segmentMask = segmentSize - 1L;
        long countLong = (length + segmentSize - 1L) >>> segmentShift;
        if (countLong > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("too many matrix segments: " + countLong);
        }
        this.length = length;
        this.segments = new double[(int)countLong][];
        for (int s = 0; s < segments.length; s++) {
            long remaining = length - ((long)s << segmentShift);
            segments[s] = new double[(int)Math.min(segmentSize, remaining)];
        }
    }

    long length() { return length; }
    int segmentShift() { return segmentShift; }
    double[][] segments() { return segments; }
    long bytes() { return length * Double.BYTES; }

    double get(long index) {
        return segments[(int)(index >>> segmentShift)][(int)(index & segmentMask)];
    }

    void set(long index, double value) {
        segments[(int)(index >>> segmentShift)][(int)(index & segmentMask)] = value;
    }

    void add(long index, double value) {
        segments[(int)(index >>> segmentShift)][(int)(index & segmentMask)] += value;
    }
}
