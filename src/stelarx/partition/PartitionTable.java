package stelarx.partition;

import stelarx.Logging;
import stelarx.cluster.ClusterHash;
import stelarx.hash.PrefixHashArrays;
import stelarx.pro.UniqueTaxonSubtreeHashes;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.ProgressBar;

import java.util.*;

/**
 * Table of unique rooted gene-tree child bipartitions with their frequencies.
 *
 * For each biological speciation node u (root included) we extract:
 *   part1 = sub(left(u))    range [L.start, L.end)   -- left subtree
 *   part2 = sub(right(u))   range [R.start, R.end)   -- right subtree
 *   part3 = Lg \ sub(u)     complement of [u.start, u.end) w.r.t. Lg
 *
 * The legacy third part stores the taxa outside u only to preserve the compact
 * data ABI used by the optimized intersection engines. STELAR-X ignores it.
 *
 * In STELAR-Pro, every part hash contains each species once. The subtree and
 * outside-subtree hashes come from the shared small-to-large tree index.
 * Deduplication: PartitionHash is order-invariant over (part1, part2).
 */
public class PartitionTable {

    public static final class Entry {
        public final PartitionHash hash;
        public final Partition     exemplar;
        public int                 frequency;

        Entry(PartitionHash h, Partition p) { this.hash = h; this.exemplar = p; this.frequency = 1; }
    }

    private final Map<PartitionHash, Entry> table = new HashMap<>();
    private final int m;
    private final UniqueTaxonSubtreeHashes uniqueTaxonHashes;
    private boolean hasPoly = false;   // true once any d>3 (polytomous) partition is stored

    /** True iff any extracted partition is polytomous (d > 3). */
    public boolean hasPolytomousPartitions() { return hasPoly; }

    // -------------------------------------------------------------------------

    public PartitionTable(List<Tree> trees, PrefixHashArrays pref) {
        this(trees, pref, null);
    }

    /** STELAR-Pro entry point with duplicate-invariant partition hashes. */
    public PartitionTable(List<Tree> trees, PrefixHashArrays pref,
                          UniqueTaxonSubtreeHashes uniqueTaxonHashes) {
        long t0 = System.nanoTime();
        this.uniqueTaxonHashes = uniqueTaxonHashes;
        this.m = uniqueTaxonHashes == null
            ? pref.numSeeds() : uniqueTaxonHashes.numSeeds();

        int totalCandidates = 0;
        int treesDone = 0;
        ProgressBar bar = new ProgressBar("Tripartition extraction", trees.size());
        for (Tree tree : trees) {
            totalCandidates += extractFromTree(tree, pref);
            bar.update(++treesDone);
        }
        bar.done();

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Rooted partition extraction: %d candidates -> %d unique child bipartitions in %d ms",
            totalCandidates, table.size(), ms);
    }

    private int extractFromTree(Tree tree, PrefixHashArrays pref) {
        int ti = tree.treeIndex;
        int L  = tree.leafCount;
        int[] count = {0};
        extractNode(tree.root, ti, L, pref, count);
        return count[0];
    }

    /**
     * Recurse post-order and register partitions rooted at speciation nodes.
     */
    private void extractNode(TreeNode node, int ti, int L,
                              PrefixHashArrays pref, int[] count) {
        if (node.isLeaf()) return;
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) extractNode(child, ti, L, pref, count);
        } else {
            extractNode(node.left,  ti, L, pref, count);
            extractNode(node.right, ti, L, pref, count);
        }

        // Descendant speciations remain valid even below a skipped duplication.
        if (!node.isSpeciation()) return;

        // ── Polytomous node: d = k+1 partition (k child subtrees + complement) ──
        if (node.isPolytomous()) {
            int k = node.children.length;
            int d = k + 1;
            ClusterHash[] hashes = new ClusterHash[d];
            int[] sizes      = new int[d];
            int[] partStarts = new int[k];
            int[] partEnds   = new int[k];
            for (int i = 0; i < k; i++) {
                TreeNode c = node.children[i];
                int cs = c.rangeStart, ce = c.rangeEnd;
                partStarts[i] = cs;
                partEnds[i] = ce;
                hashes[i] = childHash(ti, node, i, c, pref);
                sizes[i] = hashes[i].size;
            }
            ClusterHash complement = complementHash(ti, node, L, pref);
            int szC = complement.size;
            if (szC == 0) return;                               // root only — already skipped
            sizes[d - 1]  = szC;
            hashes[d - 1] = complement;

            PartitionHash ph = new PartitionHash(hashes);
            Entry existing = table.get(ph);
            if (existing != null) {
                existing.frequency++;
            } else {
                Partition p = new Partition(hashes, sizes, partStarts, partEnds, ti);
                table.put(ph, new Entry(ph, p));
                hasPoly = true;
            }
            count[0]++;
            return;
        }

        // ── Binary node ──
        int lStart = node.left.rangeStart,  lEnd = node.left.rangeEnd;
        int rStart = node.right.rangeStart, rEnd = node.right.rangeEnd;

        ClusterHash h1 = childHash(ti, node, 0, node.left, pref);
        ClusterHash h2 = childHash(ti, node, 1, node.right, pref);
        ClusterHash h3 = complementHash(ti, node, L, pref);
        int sz1 = h1.size;
        int sz2 = h2.size;
        int sz3 = h3.size;

        PartitionHash ph = new PartitionHash(h1, h2, h3);

        Entry existing = table.get(ph);
        if (existing != null) {
            existing.frequency++;
        } else {
            Partition p = new Partition(h1, h2, h3, sz1, sz2, sz3,
                                        ti, lStart, lEnd, rStart, rEnd);
            table.put(ph, new Entry(ph, p));
        }
        count[0]++;
    }

    /** Compute a ClusterHash for a subtree range or its complement. */
    private ClusterHash buildHash(int ti, int lo, int hi, boolean complement,
                                   int size, PrefixHashArrays pref) {
        long[] rawSums = new long[m], rawXors = new long[m];
        for (int s = 0; s < m; s++) {
            rawSums[s] = complement ? pref.compSum(ti, s, lo, hi)
                                    : pref.rangeSum(ti, s, lo, hi);
            rawXors[s] = complement ? pref.compXor(ti, s, lo, hi)
                                    : pref.rangeXor(ti, s, lo, hi);
        }
        return new ClusterHash(rawSums, rawXors, size, m);
    }

    private ClusterHash childHash(int treeIndex, TreeNode speciationRoot,
                                  int childIndex, TreeNode child,
                                  PrefixHashArrays pref) {
        return uniqueTaxonHashes == null
            ? buildHash(treeIndex, child.rangeStart, child.rangeEnd, false,
                child.rangeSize(), pref)
            : uniqueTaxonHashes.getChild(treeIndex, speciationRoot, childIndex);
    }

    private ClusterHash complementHash(int treeIndex, TreeNode node, int leafCount,
                                       PrefixHashArrays pref) {
        return uniqueTaxonHashes == null
            ? buildHash(treeIndex, node.rangeStart, node.rangeEnd, true,
                leafCount - node.rangeSize(), pref)
            : uniqueTaxonHashes.getComplement(treeIndex, node);
    }

    // -------------------------------------------------------------------------

    public Entry get(PartitionHash ph) { return table.get(ph); }
    public int size()                  { return table.size(); }
    public Collection<Entry> entries() { return table.values(); }
    public boolean isDuplicateInvariant() { return uniqueTaxonHashes != null; }
}
