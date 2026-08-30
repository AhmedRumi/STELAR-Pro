package stelarx.partition;

import stelarx.cluster.ClusterHash;

/**
 * A compact rooted gene-tree child partition. The first d-1 groups are the
 * actual children of a node; the last group is a legacy complement slot that
 * may be empty at the supplied root and is ignored by the STELAR-X objective.
 *
 * For an internal node u of a rooted gene tree g with k children:
 *   d = k + 1
 *   M_i = subtree of child i  (range [partStarts[i], partEnds[i]))  for i = 0..d-2
 *   M_{d-1} = Lg \ sub(u)     -- the "everything else in this tree" (complement)
 *
 * Binary nodes are the special case k=2, d=3.  They use the named scalar fields
 * directly and leave the general-polytomy arrays null, avoiding four redundant
 * per-partition arrays on the overwhelmingly common binary path.
 *
 * For a general polytomy the complement part is not stored as a range (it is
 * non-contiguous in general); its hash and size occupy the last array slot.
 * Binary partitions retain those values in {@code hash3}/{@code size3} instead.
 */
public final class Partition {

    /** Number of parts: 3 for binary, k+1 for a polytomous node with k children. */
    public final int d;

    /** Part hashes [0..d-1]; {@code hashes[d-1]} is the complement. */
    public final ClusterHash[] hashes;

    /** Part sizes [0..d-1]; {@code sizes[d-1]} is the complement size. */
    public final int[] sizes;

    /** Ranges of the d-1 child parts (contiguous postorder spans), length d-1. */
    public final int[] partStarts;
    public final int[] partEnds;

    /** Index of the gene tree this partition came from (of the first exemplar). */
    public final int treeIndex;

    // ── Binary (d=3) convenience aliases — valid/used only for binary partitions ──
    public final ClusterHash hash1, hash2, hash3;
    public final int size1, size2, size3;
    public final int leftStart, leftEnd, rightStart, rightEnd;

    /** Binary constructor (d=3). */
    public Partition(ClusterHash h1, ClusterHash h2, ClusterHash h3,
                     int sz1, int sz2, int sz3,
                     int treeIndex,
                     int leftStart, int leftEnd,
                     int rightStart, int rightEnd) {
        this.d = 3;
        this.hashes = null;
        this.sizes  = null;
        this.partStarts = null;
        this.partEnds   = null;
        this.treeIndex = treeIndex;
        // aliases
        this.hash1 = h1; this.hash2 = h2; this.hash3 = h3;
        this.size1 = sz1; this.size2 = sz2; this.size3 = sz3;
        this.leftStart = leftStart; this.leftEnd = leftEnd;
        this.rightStart = rightStart; this.rightEnd = rightEnd;
    }

    /**
     * General d-partition constructor (d ≥ 3).  {@code hashes}/{@code sizes} have
     * length d ({@code [d-1]} = complement); {@code partStarts}/{@code partEnds}
     * have length d-1 (the child ranges).
     */
    public Partition(ClusterHash[] hashes, int[] sizes,
                     int[] partStarts, int[] partEnds, int treeIndex) {
        this.d = hashes.length;
        this.hashes = hashes;
        this.sizes  = sizes;
        this.partStarts = partStarts;
        this.partEnds   = partEnds;
        this.treeIndex = treeIndex;
        // aliases (meaningful only when d==3; for d>3 they mirror first/second/last
        // parts and must NOT be used — correct consumers branch on d)
        this.hash1 = hashes[0];
        this.hash2 = hashes[1];
        this.hash3 = hashes[d - 1];
        this.size1 = sizes[0];
        this.size2 = sizes[1];
        this.size3 = sizes[d - 1];
        this.leftStart  = partStarts[0];
        this.leftEnd    = partEnds[0];
        this.rightStart = partStarts[d - 2];
        this.rightEnd   = partEnds[d - 2];
    }

    public boolean isPolytomous() { return d > 3; }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Partition{sz=");
        if (d == 3) {
            sb.append(size1).append('|').append(size2).append('|').append(size3);
        } else {
            for (int i = 0; i < d; i++) {
                if (i > 0) sb.append('|');
                sb.append(sizes[i]);
            }
        }
        sb.append(", d=").append(d).append(", t=").append(treeIndex).append('}');
        return sb.toString();
    }
}
