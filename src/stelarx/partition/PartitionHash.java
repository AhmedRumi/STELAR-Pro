package stelarx.partition;

import stelarx.cluster.ClusterHash;

/**
 * Order-invariant hash key for a rooted gene-tree bipartition/partition.
 *
 * <p>STELAR-Pro keys a resolved rooted bipartition only by its unordered child
 * pair {@code {M1,M2}}. The legacy constructors retain the outside-subtree slot
 * for STELAR-Pro paths.</p>
 *
 * BINARY (d=3): a tripartition (M1|M2|M3) is identified by the unordered pair {M1,M2}
 * plus M3.  For complete trees M3 = S\M1\M2 is determined by (M1,M2); for incomplete
 * trees two nodes can share M1,M2 but differ in M3, so h3 is included.  This binary
 * path is kept BYTE-IDENTICAL to the pre-polytomy implementation.
 *
 * POLYTOMOUS (d≥4): a rooted d-partition M0|…|M_{d-1} is identified by the
 * unordered multiset of its d-1 child hashes plus its distinguished final
 * complement hash.  Only the children contribute to the rooted-triplet weight,
 * so exchanging a child with the complement is not a symmetry.  We sort the
 * child fingerprints lexicographically and append the complement fingerprint.
 *
 * The two representations never collide: binary partitions (d=3) and polytomous ones
 * (d≥4) come from disjoint parser paths and {@code equals} short-circuits on {@code d}.
 */
public final class PartitionHash {

    private final int cachedHashCode;
    private final int d;       // number of parts (3 for binary, k+1 for polytomous)

    /** Flattened canonical fingerprints for both binary and general partitions. */
    private final long[] data;

    /** STELAR-Pro rooted-bipartition identity: unordered child pair {A,B}. */
    public PartitionHash(ClusterHash a, ClusterHash b) {
        boolean aFirst = compare(a, b) <= 0;
        ClusterHash first  = aFirst ? a : b;
        ClusterHash second = aFirst ? b : a;

        int m = a.sums.length;
        this.d = 2;
        data = new long[4 * m];
        for (int s = 0; s < m; s++) {
            data[s]         = first.sums[s];
            data[s + m]     = first.xors[s];
            data[s + 2 * m] = second.sums[s];
            data[s + 3 * m] = second.xors[s];
        }
        int h = 1;
        for (long v : data) h = 31 * h + Long.hashCode(v);
        this.cachedHashCode = h;
    }

    public PartitionHash(ClusterHash a, ClusterHash b, ClusterHash c) {
        // Decide ordering of the two explicit parts (a,b are interchangeable)
        boolean aFirst = compare(a, b) <= 0;
        ClusterHash first  = aFirst ? a : b;
        ClusterHash second = aFirst ? b : a;

        int m = a.sums.length;
        this.d = 3;
        data = new long[6 * m];
        for (int s = 0; s < m; s++) {
            data[s]         = first.sums[s];
            data[s + m]     = first.xors[s];
            data[s + 2 * m] = second.sums[s];
            data[s + 3 * m] = second.xors[s];
            data[s + 4 * m] = c.sums[s];
            data[s + 5 * m] = c.xors[s];
        }

        int h = 1;
        for (long v : data) h = 31 * h + Long.hashCode(v);
        this.cachedHashCode = h;
    }

    /**
     * General rooted d-partition key (used for polytomous nodes, d≥4):
     * order-invariant over the d-1 children, with {@code parts[d-1]} retained as
     * the distinguished complement.
     */
    public PartitionHash(ClusterHash[] parts) {
        int dd = parts.length;
        int m = parts[0].sums.length;
        long[][] fps = new long[dd][2 * m];
        for (int i = 0; i < dd; i++) {
            for (int s = 0; s < m; s++) {
                fps[i][s]     = parts[i].sums[s];
                fps[i][s + m] = parts[i].xors[s];
            }
        }
        java.util.Arrays.sort(fps, 0, dd - 1, (x, y) -> {
            for (int s = 0; s < x.length; s++) {
                int c = Long.compareUnsigned(x[s], y[s]);
                if (c != 0) return c;
            }
            return 0;
        });
        long[] flat = new long[dd * 2 * m];
        int p = 0;
        for (long[] fp : fps) for (long v : fp) flat[p++] = v;

        this.d = dd;
        this.data = flat;
        int h = 1;
        for (long v : flat) h = 31 * h + Long.hashCode(v);
        this.cachedHashCode = h;
    }

    /** STELAR-Pro identity for an unordered collection of rooted child sets. */
    public PartitionHash(ClusterHash[] parts, int childCount) {
        if (childCount < 2 || childCount > parts.length) {
            throw new IllegalArgumentException("Invalid rooted child count");
        }
        int m = parts[0].sums.length;
        long[][] fingerprints = new long[childCount][2 * m];
        for (int i = 0; i < childCount; i++) {
            for (int s = 0; s < m; s++) {
                fingerprints[i][s] = parts[i].sums[s];
                fingerprints[i][s + m] = parts[i].xors[s];
            }
        }
        java.util.Arrays.sort(fingerprints, (x, y) -> {
            for (int s = 0; s < x.length; s++) {
                int comparison = Long.compareUnsigned(x[s], y[s]);
                if (comparison != 0) return comparison;
            }
            return 0;
        });

        this.d = childCount;
        this.data = new long[childCount * 2 * m];
        int cursor = 0;
        for (long[] fingerprint : fingerprints) {
            for (long value : fingerprint) data[cursor++] = value;
        }
        int h = 1;
        for (long value : data) h = 31 * h + Long.hashCode(value);
        this.cachedHashCode = h;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PartitionHash p)) return false;
        if (d != p.d) return false;
        return java.util.Arrays.equals(data, p.data);
    }

    @Override
    public int hashCode() { return cachedHashCode; }

    /** Lexicographic comparison of two ClusterHash objects (sums first, then xors). */
    private static int compare(ClusterHash a, ClusterHash b) {
        int m = a.sums.length;
        for (int s = 0; s < m; s++) {
            int c = Long.compareUnsigned(a.sums[s], b.sums[s]);
            if (c != 0) return c;
        }
        for (int s = 0; s < m; s++) {
            int c = Long.compareUnsigned(a.xors[s], b.xors[s]);
            if (c != 0) return c;
        }
        return 0;
    }
}
