package stelarx.pro;

import stelarx.Logging;
import stelarx.cluster.ClusterHash;
import stelarx.hash.TaxonHasher;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;

/**
 * Indexes duplicate-invariant rooted bipartitions at biological speciation nodes.
 *
 * <p>Each tree is processed in postorder. Every node temporarily carries a taxon
 * bag so its descendants can be merged upward, but only a speciation node creates
 * a persistent hash record. That record contains the candidate subtree hash and
 * the child/complement parts needed to represent that speciation-driven rooted
 * bipartition. Node origin is irrelevant: a refinement introduced during
 * polytomy resolution is indexed whenever the root/tag stage marks it as
 * speciation. Duplication nodes, leaves, and untagged fallback nodes receive no
 * index entries. Hash sums and XORs change only on first insertion, so each
 * species contributes once even when the subtree contains several gene copies.</p>
 *
 * <p>For a tree with {@code L} leaf occurrences, {@code V} nodes, and {@code m}
 * hash seeds, expected construction time is {@code O(m * (V + L log L))} with
 * open-addressing maps. Working sets contain {@code O(L)} entries. For resolved
 * binary trees, immutable speciation records require {@code O(m * S)} space for
 * {@code S} speciation nodes until their consumers have been built.</p>
 */
public final class UniqueTaxonSubtreeHashes {
    private List<IdentityHashMap<TreeNode, SpeciationBipartitionHashes>> treeHashes;
    private final ClusterHash allTaxaHash;
    private final int m;

    public UniqueTaxonSubtreeHashes(List<Tree> trees, TaxonHasher hasher) {
        long started = System.nanoTime();
        this.m = hasher.numSeeds();
        this.treeHashes = new ArrayList<>(trees.size());
        for (int i = 0; i < trees.size(); i++) {
            treeHashes.add(null);
        }

        long bipartitionCount = 0;
        int[] occurrenceScratch = new int[hasher.numTaxa()];
        for (Tree tree : trees) {
            if (tree.treeIndex < 0 || tree.treeIndex >= trees.size()
                    || treeHashes.get(tree.treeIndex) != null) {
                throw new IllegalArgumentException(
                    "Gene-tree indices must be unique and contiguous from zero");
            }
            TreeIndex index = buildTreeIndex(tree, hasher, occurrenceScratch);
            treeHashes.set(tree.treeIndex, index.bipartitions());
            bipartitionCount += index.bipartitions().size();
        }
        this.allTaxaHash = hashAllTaxa(hasher);

        long elapsedMs = (System.nanoTime() - started) / 1_000_000;
        Logging.info("Built %d duplicate-invariant speciation-rooted bipartition "
            + "hash records in %d gene tree(s) by small-to-large merging in %d ms",
            bipartitionCount, trees.size(), elapsedMs);
    }

    /** Candidate subtree hash; valid only for a biological speciation node. */
    public ClusterHash get(int treeIndex, TreeNode node) {
        return requireSpeciationHash(treeIndex, node).subtree();
    }

    /** Hash of one child part of a speciation-rooted bipartition. */
    public ClusterHash getChild(int treeIndex, TreeNode node, int childIndex) {
        ClusterHash[] children = requireSpeciationHash(treeIndex, node).children();
        if (childIndex < 0 || childIndex >= children.length) {
            throw new IllegalArgumentException("Invalid child index " + childIndex);
        }
        return children[childIndex];
    }

    /** Distinct-taxon hash of the leaves outside a speciation-rooted subtree. */
    public ClusterHash getComplement(int treeIndex, TreeNode node) {
        return requireSpeciationHash(treeIndex, node).complement();
    }

    /** Used by verification to assert that non-speciation nodes were not indexed. */
    public boolean contains(int treeIndex, TreeNode node) {
        return treeHashes.get(treeIndex).containsKey(node);
    }

    public ClusterHash allTaxaHash() { return allTaxaHash; }
    public int numSeeds() { return m; }

    /**
     * Release the persistent node-to-hash indexes after all S1 consumers finish
     * construction. Canonical hashes already copied into X, partitions, and DP
     * transitions remain valid.
     */
    public void release() {
        treeHashes = List.of();
    }

    private SpeciationBipartitionHashes requireSpeciationHash(
            int treeIndex, TreeNode node) {
        SpeciationBipartitionHashes hash = treeHashes.get(treeIndex).get(node);
        if (hash == null) {
            throw new IllegalArgumentException(
                "Node is not an indexed speciation root in tree " + treeIndex);
        }
        return hash;
    }

    private TreeIndex buildTreeIndex(
            Tree tree, TaxonHasher hasher, int[] treeOccurrences) {
        IdentityHashMap<TreeNode, SpeciationBipartitionHashes> hashes =
            new IdentityHashMap<>();
        IdentityHashMap<TreeNode, TaxonBag> liveBags = new IdentityHashMap<>();

        // Counts are transient. They distinguish a species wholly contained in
        // sub(u) from one that also has a copy in Lg \ sub(u), which is required
        // for the third part of a rooted gene-tree partition.
        int[] touchedTaxa = new int[Math.min(tree.leafCount, treeOccurrences.length)];
        int touchedCount = 0;
        for (int taxon : tree.postorderArray) {
            if (treeOccurrences[taxon]++ == 0) touchedTaxa[touchedCount++] = taxon;
        }
        ClusterHash treeTaxaHash = hashPresentTaxa(touchedTaxa, touchedCount, hasher);

        // Build an iterative postorder so hashing does not add recursion-depth
        // limits on highly unbalanced gene trees.
        ArrayDeque<TreeNode> pending = new ArrayDeque<>();
        ArrayDeque<TreeNode> postorder = new ArrayDeque<>();
        pending.push(tree.root);
        while (!pending.isEmpty()) {
            TreeNode node = pending.pop();
            postorder.push(node);
            if (node.isLeaf()) continue;
            if (node.isPolytomous()) {
                for (TreeNode child : node.children) pending.push(child);
            } else {
                pending.push(node.left);
                pending.push(node.right);
            }
        }

        while (!postorder.isEmpty()) {
            TreeNode node = postorder.pop();
            TaxonBag bag;
            ClusterHash[] childHashes = null;
            if (node.isLeaf()) {
                bag = new TaxonBag(node.taxonId, treeOccurrences, hasher);
            } else if (node.isPolytomous()) {
                bag = null;
                if (node.isSpeciation()) {
                    childHashes = new ClusterHash[node.children.length];
                }
                for (int childIndex = 0; childIndex < node.children.length; childIndex++) {
                    TreeNode child = node.children[childIndex];
                    TaxonBag childBag = removeChildBag(liveBags, child, tree.treeIndex);
                    if (childHashes != null) {
                        childHashes[childIndex] = childBag.snapshot(m);
                    }
                    bag = bag == null ? childBag
                        : mergeSmallIntoLarge(bag, childBag, treeOccurrences, hasher);
                }
            } else {
                TaxonBag left = removeChildBag(liveBags, node.left, tree.treeIndex);
                TaxonBag right = removeChildBag(liveBags, node.right, tree.treeIndex);
                if (node.isSpeciation()) {
                    childHashes = new ClusterHash[]{left.snapshot(m), right.snapshot(m)};
                }
                bag = mergeSmallIntoLarge(left, right, treeOccurrences, hasher);
            }

            // The event tag is the only admission rule. Thus a resolver-created
            // node tagged as speciation enters the index exactly like any other
            // speciation. Every non-speciation node retains only its mutable bag;
            // if it is a child part, its hash belongs solely to the speciation
            // parent's record and is not a standalone candidate.
            if (node.isSpeciation()) {
                hashes.put(node, new SpeciationBipartitionHashes(
                    bag.snapshot(m), childHashes,
                    bag.complementSnapshot(treeTaxaHash, m)));
            }
            liveBags.put(node, bag);
        }

        TaxonBag rootBag = liveBags.remove(tree.root);
        if (rootBag == null || !liveBags.isEmpty()) {
            throw new IllegalStateException(
                "Incomplete small-to-large merge for gene tree " + tree.treeIndex);
        }
        // No working taxon set survives beyond this tree. Only compact immutable
        // ClusterHash values remain in the persistent index.
        rootBag.release();
        // Reuse one O(n) occurrence array across all trees without paying O(n)
        // clearing time for each incomplete tree.
        for (int i = 0; i < touchedCount; i++) treeOccurrences[touchedTaxa[i]] = 0;
        return new TreeIndex(hashes);
    }

    private static TaxonBag removeChildBag(
            IdentityHashMap<TreeNode, TaxonBag> liveBags,
            TreeNode child, int treeIndex) {
        TaxonBag bag = liveBags.remove(child);
        if (bag == null) {
            throw new IllegalStateException(
                "Missing child taxon set in gene tree " + treeIndex);
        }
        return bag;
    }

    /** Merge the smaller distinct-taxon set into the larger and release it. */
    private static TaxonBag mergeSmallIntoLarge(
            TaxonBag first, TaxonBag second, int[] treeOccurrences,
            TaxonHasher hasher) {
        TaxonBag large = first;
        TaxonBag small = second;
        if (large.size() < small.size()) {
            large = second;
            small = first;
        }
        large.absorbAndRelease(small, treeOccurrences, hasher);
        return large;
    }

    private ClusterHash hashPresentTaxa(int[] touchedTaxa, int touchedCount,
                                        TaxonHasher hasher) {
        long[] sums = new long[m];
        long[] xors = new long[m];
        for (int i = 0; i < touchedCount; i++) {
            int taxon = touchedTaxa[i];
            for (int seed = 0; seed < m; seed++) {
                long value = hasher.get(seed, taxon);
                sums[seed] += value;
                xors[seed] ^= value;
            }
        }
        return new ClusterHash(sums, xors, touchedCount, m);
    }

    private ClusterHash hashAllTaxa(TaxonHasher hasher) {
        long[] sums = new long[m];
        long[] xors = new long[m];
        for (int taxon = 0; taxon < hasher.numTaxa(); taxon++) {
            for (int seed = 0; seed < m; seed++) {
                long value = hasher.get(seed, taxon);
                sums[seed] += value;
                xors[seed] ^= value;
            }
        }
        return new ClusterHash(sums, xors, hasher.numTaxa(), m);
    }

    /** Mutable distinct-taxon map plus its incrementally maintained hashes. */
    private static final class TaxonBag {
        private IntCountMap taxa;
        private long[] sums;
        private long[] xors;
        private long[] completeSums;
        private long[] completeXors;
        private int completeCount;

        TaxonBag(int taxon, int[] treeOccurrences, TaxonHasher hasher) {
            taxa = new IntCountMap();
            taxa.addTo(taxon, 1);
            sums = new long[hasher.numSeeds()];
            xors = new long[hasher.numSeeds()];
            completeSums = new long[hasher.numSeeds()];
            completeXors = new long[hasher.numSeeds()];
            addHash(sums, xors, taxon, hasher);
            if (treeOccurrences[taxon] == 1) markComplete(taxon, hasher);
        }

        int size() { return taxa.size(); }

        void absorbAndRelease(TaxonBag donor, int[] treeOccurrences,
                              TaxonHasher hasher) {
            int[] donorKeys = donor.taxa.keys();
            int[] donorCounts = donor.taxa.counts();
            for (int slot = 0; slot < donorKeys.length; slot++) {
                int encodedTaxon = donorKeys[slot];
                if (encodedTaxon == IntCountMap.EMPTY) continue;
                int taxon = encodedTaxon - 1;
                int donorCount = donorCounts[slot];
                int oldCount = taxa.addTo(taxon, donorCount);
                int mergedCount = oldCount + donorCount;
                if (mergedCount > treeOccurrences[taxon]) {
                    throw new IllegalStateException(
                        "Subtree copy count exceeds its gene-tree total");
                }
                if (oldCount == 0) addHash(sums, xors, taxon, hasher);
                if (oldCount < treeOccurrences[taxon]
                        && mergedCount == treeOccurrences[taxon]) {
                    markComplete(taxon, hasher);
                }
            }
            donor.release();
        }

        ClusterHash snapshot(int seeds) {
            return new ClusterHash(sums, xors, taxa.size(), seeds);
        }

        ClusterHash complementSnapshot(ClusterHash treeTaxaHash, int seeds) {
            long[] outsideSums = new long[seeds];
            long[] outsideXors = new long[seeds];
            for (int seed = 0; seed < seeds; seed++) {
                outsideSums[seed] = treeTaxaHash.sums[seed] - completeSums[seed];
                outsideXors[seed] = treeTaxaHash.xors[seed] ^ completeXors[seed];
            }
            return new ClusterHash(outsideSums, outsideXors,
                treeTaxaHash.size - completeCount, seeds);
        }

        private static void addHash(long[] targetSums, long[] targetXors,
                                    int taxon, TaxonHasher hasher) {
            for (int seed = 0; seed < targetSums.length; seed++) {
                long value = hasher.get(seed, taxon);
                targetSums[seed] += value;
                targetXors[seed] ^= value;
            }
        }

        private void markComplete(int taxon, TaxonHasher hasher) {
            addHash(completeSums, completeXors, taxon, hasher);
            completeCount++;
        }

        /** Drop all references owned by a consumed working set. */
        void release() {
            if (taxa != null) taxa.release();
            taxa = null;
            sums = null;
            xors = null;
            completeSums = null;
            completeXors = null;
            completeCount = 0;
        }
    }

    /** Primitive open-addressing taxon-count map; keys store {@code id + 1}. */
    private static final class IntCountMap {
        static final int EMPTY = 0;
        private static final int MIN_CAPACITY = 4;

        private int[] keys = new int[MIN_CAPACITY];
        private int[] counts = new int[MIN_CAPACITY];
        private int size;
        private int resizeThreshold = threshold(MIN_CAPACITY);

        int size() { return size; }
        int[] keys() { return keys; }
        int[] counts() { return counts; }

        /** Add a positive count and return the previous count. */
        int addTo(int taxon, int count) {
            int encoded = taxon + 1;
            int slot = findSlot(encoded, keys);
            if (keys[slot] == encoded) {
                int previous = counts[slot];
                counts[slot] += count;
                return previous;
            }
            if (size + 1 > resizeThreshold) {
                resize();
                slot = findSlot(encoded, keys);
            }
            keys[slot] = encoded;
            counts[slot] = count;
            size++;
            return 0;
        }

        void release() {
            keys = null;
            counts = null;
            size = 0;
            resizeThreshold = 0;
        }

        private void resize() {
            int[] oldKeys = keys;
            int[] oldCounts = counts;
            keys = new int[oldKeys.length << 1];
            counts = new int[keys.length];
            resizeThreshold = threshold(keys.length);
            for (int oldSlot = 0; oldSlot < oldKeys.length; oldSlot++) {
                int encoded = oldKeys[oldSlot];
                if (encoded == EMPTY) continue;
                int newSlot = findSlot(encoded, keys);
                keys[newSlot] = encoded;
                counts[newSlot] = oldCounts[oldSlot];
            }
        }

        private static int findSlot(int encoded, int[] values) {
            int mask = values.length - 1;
            int slot = mix(encoded) & mask;
            while (values[slot] != EMPTY && values[slot] != encoded) {
                slot = (slot + 1) & mask;
            }
            return slot;
        }

        private static int threshold(int capacity) {
            return capacity * 2 / 3;
        }

        private static int mix(int value) {
            value ^= value >>> 16;
            value *= 0x7feb352d;
            value ^= value >>> 15;
            value *= 0x846ca68b;
            return value ^ (value >>> 16);
        }
    }

    private record SpeciationBipartitionHashes(
        ClusterHash subtree, ClusterHash[] children, ClusterHash complement) {}

    private record TreeIndex(
        IdentityHashMap<TreeNode, SpeciationBipartitionHashes> bipartitions) {}
}
