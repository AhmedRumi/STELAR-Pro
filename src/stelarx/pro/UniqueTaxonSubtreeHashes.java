package stelarx.pro;

import stelarx.cluster.ClusterHash;
import stelarx.hash.TaxonHasher;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.List;

/**
 * Indexes each gene-tree node by the hash and size of its distinct species set.
 * An offline last-occurrence scan plus Fenwick range queries avoids materializing
 * a species bit set at every node while ensuring repeated copies contribute once.
 */
public final class UniqueTaxonSubtreeHashes {
    private final List<IdentityHashMap<TreeNode, ClusterHash>> treeHashes;
    private final ClusterHash allTaxaHash;
    private final int m;

    public UniqueTaxonSubtreeHashes(List<Tree> trees, TaxonHasher hasher) {
        this.m = hasher.numSeeds();
        this.treeHashes = new ArrayList<>(trees.size());
        for (int i = 0; i < trees.size(); i++) treeHashes.add(null);
        for (Tree tree : trees) {
            if (tree.treeIndex < 0 || tree.treeIndex >= trees.size()
                    || treeHashes.get(tree.treeIndex) != null) {
                throw new IllegalArgumentException(
                    "Gene-tree indices must be unique and contiguous from zero");
            }
            treeHashes.set(tree.treeIndex, buildTreeIndex(tree, hasher));
        }
        this.allTaxaHash = hashAllTaxa(hasher);
    }

    public ClusterHash get(int treeIndex, TreeNode node) {
        ClusterHash hash = treeHashes.get(treeIndex).get(node);
        if (hash == null) {
            throw new IllegalArgumentException(
                "Node is absent from duplicate-invariant tree index " + treeIndex);
        }
        return hash;
    }

    public ClusterHash allTaxaHash() { return allTaxaHash; }
    public int numSeeds() { return m; }

    private IdentityHashMap<TreeNode, ClusterHash> buildTreeIndex(
            Tree tree, TaxonHasher hasher) {
        int leafCount = tree.leafCount;
        @SuppressWarnings("unchecked")
        List<TreeNode>[] endingAt = new List[leafCount + 1];
        IdentityHashMap<TreeNode, ClusterHash> hashes = new IdentityHashMap<>();

        // Collect subtree-range queries top-down; all ranges are contiguous.
        ArrayDeque<TreeNode> pending = new ArrayDeque<>();
        pending.push(tree.root);
        while (!pending.isEmpty()) {
            TreeNode node = pending.pop();
            if (endingAt[node.rangeEnd] == null) endingAt[node.rangeEnd] = new ArrayList<>();
            endingAt[node.rangeEnd].add(node);
            if (!node.isLeaf()) {
                if (node.isPolytomous()) {
                    for (TreeNode child : node.children) pending.push(child);
                } else {
                    pending.push(node.right);
                    pending.push(node.left);
                }
            }
        }

        FenwickLong[] sums = new FenwickLong[m];
        FenwickLong[] xors = new FenwickLong[m];
        for (int seed = 0; seed < m; seed++) {
            sums[seed] = new FenwickLong(leafCount, false);
            xors[seed] = new FenwickLong(leafCount, true);
        }
        FenwickInt distinctCounts = new FenwickInt(leafCount);
        int[] lastPosition = new int[hasher.numTaxa()];
        Arrays.fill(lastPosition, -1);

        for (int position = 0; position < leafCount; position++) {
            int taxon = tree.postorderArray[position];
            int previous = lastPosition[taxon];
            if (previous >= 0) distinctCounts.add(previous, -1);
            distinctCounts.add(position, 1);
            for (int seed = 0; seed < m; seed++) {
                long value = hasher.get(seed, taxon);
                if (previous >= 0) {
                    sums[seed].add(previous, -value);
                    xors[seed].add(previous, value);
                }
                sums[seed].add(position, value);
                xors[seed].add(position, value);
            }
            lastPosition[taxon] = position;

            List<TreeNode> queries = endingAt[position + 1];
            if (queries == null) continue;
            for (TreeNode node : queries) {
                long[] nodeSums = new long[m];
                long[] nodeXors = new long[m];
                for (int seed = 0; seed < m; seed++) {
                    nodeSums[seed] = sums[seed].range(node.rangeStart, node.rangeEnd);
                    nodeXors[seed] = xors[seed].range(node.rangeStart, node.rangeEnd);
                }
                int size = distinctCounts.range(node.rangeStart, node.rangeEnd);
                hashes.put(node, new ClusterHash(nodeSums, nodeXors, size, m));
            }
        }
        return hashes;
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

    private static final class FenwickLong {
        private final long[] values;
        private final boolean xor;

        FenwickLong(int size, boolean xor) {
            this.values = new long[size + 1];
            this.xor = xor;
        }

        void add(int position, long value) {
            for (int index = position + 1; index < values.length; index += index & -index) {
                if (xor) values[index] ^= value;
                else values[index] += value;
            }
        }

        long range(int start, int end) {
            return xor ? prefix(end) ^ prefix(start) : prefix(end) - prefix(start);
        }

        private long prefix(int end) {
            long result = 0L;
            for (int index = end; index > 0; index -= index & -index) {
                if (xor) result ^= values[index];
                else result += values[index];
            }
            return result;
        }
    }

    private static final class FenwickInt {
        private final int[] values;

        FenwickInt(int size) { this.values = new int[size + 1]; }

        void add(int position, int value) {
            for (int index = position + 1; index < values.length; index += index & -index) {
                values[index] += value;
            }
        }

        int range(int start, int end) { return prefix(end) - prefix(start); }

        private int prefix(int end) {
            int result = 0;
            for (int index = end; index > 0; index -= index & -index) result += values[index];
            return result;
        }
    }
}
