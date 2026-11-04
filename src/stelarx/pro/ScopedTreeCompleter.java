package stelarx.pro;

import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SortedRowsBuilder;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.Threading;

import java.util.Arrays;
import java.util.List;

/** Completes DISCO components only within each original pre-cut lineage. */
public final class ScopedTreeCompleter {
    private final int n;
    private final int[] rows;
    private final int[][] packedRows;

    private ScopedTreeCompleter(SimilarityMatrix similarity, int n,
                                boolean anyMissing) {
        this.n = n;
        this.rows = anyMissing && !similarity.isPacked()
            ? SortedRowsBuilder.buildCPU(similarity.dist, n) : null;
        this.packedRows = anyMissing && similarity.isPacked()
            ? SortedRowsBuilder.buildPackedCPU(similarity) : null;
    }

    public static ScopedTreeCompleter prepare(SimilarityMatrix similarity,
                                              List<Tree> trees, int n) {
        boolean anyMissing = false;
        for (Tree tree : trees) {
            if (!(tree instanceof DiscoDecomposer.ScopedTree scoped)) {
                throw new IllegalArgumentException("Expected scoped DISCO component");
            }
            if (scoped.missingCount() > 0) {
                anyMissing = true;
                break;
            }
        }
        return new ScopedTreeCompleter(similarity, n, anyMissing);
    }

    /** Complete independent components in parallel within the live batch. */
    public List<Tree> completeBatch(List<Tree> batch) {
        Tree[] result = batch.toArray(new Tree[0]);
        Threading.processRangeParallel(result.length, i -> {
            DiscoDecomposer.ScopedTree scoped =
                (DiscoDecomposer.ScopedTree) result[i];
            if (scoped.missingCount() > 0) result[i] = complete(scoped);
        });
        return Arrays.asList(result);
    }

    private Tree complete(DiscoDecomposer.ScopedTree source) {
        TreeNode root = source.root;
        TreeNode[] taxonNode = new TreeNode[n];
        indexLeaves(root, taxonNode);

        // The packed bits represent taxa lost to cuts below this root, not
        // taxa in ancestral speciation siblings or absent from the gene tree.
        long[] missing = source.missingWords();
        for (int wi = 0; wi < missing.length; wi++) {
            long word = missing[wi];
            while (word != 0L) {
                int x = (wi << 6) + Long.numberOfTrailingZeros(word);
                word &= word - 1;
                if (x >= n || taxonNode[x] != null) {
                    throw new IllegalStateException("Invalid within-lineage missing taxon " + x);
                }
                int anchor = findAnchor(x, taxonNode);
                TreeNode oldLeaf = taxonNode[anchor];
                TreeNode parent = oldLeaf.parent;
                TreeNode newLeaf = new TreeNode();
                newLeaf.taxonId = x;
                TreeNode join = new TreeNode();
                join.left = oldLeaf;
                join.right = newLeaf;
                join.parent = parent;
                oldLeaf.parent = join;
                newLeaf.parent = join;
                if (parent == null) root = join;
                else if (parent.left == oldLeaf) parent.left = join;
                else if (parent.right == oldLeaf) parent.right = join;
                else throw new IllegalStateException("Broken component parent link");
                taxonNode[x] = newLeaf;
            }
        }

        int expected = source.leafCount + source.missingCount();
        int[] order = new int[expected];
        int[] cursor = {0};
        assignRanges(root, order, cursor);
        if (cursor[0] != expected) {
            throw new IllegalStateException("Scoped completion leaf count mismatch");
        }
        int[] positions = new int[n];
        Arrays.fill(positions, -1);
        for (int i = 0; i < expected; i++) {
            if (positions[order[i]] != -1) {
                throw new IllegalStateException("Scoped completion duplicated taxon");
            }
            positions[order[i]] = i;
        }
        return new Tree(source.treeIndex, root, order, positions, expected, n);
    }

    private int findAnchor(int taxon, TreeNode[] present) {
        if (packedRows != null) {
            for (int candidate : packedRows[taxon]) {
                if (present[candidate] != null) return candidate;
            }
        } else {
            int offset = taxon * n;
            for (int rank = 0; rank < n; rank++) {
                int candidate = rows[offset + rank];
                if (present[candidate] != null) return candidate;
            }
        }
        throw new IllegalStateException("No anchor in DISCO component");
    }

    private static void indexLeaves(TreeNode node, TreeNode[] byTaxon) {
        if (node.isLeaf()) {
            byTaxon[node.taxonId] = node;
        } else {
            indexLeaves(node.left, byTaxon);
            indexLeaves(node.right, byTaxon);
        }
    }

    private static void assignRanges(TreeNode node, int[] order, int[] cursor) {
        node.rangeStart = cursor[0];
        if (node.isLeaf()) {
            order[cursor[0]++] = node.taxonId;
        } else {
            assignRanges(node.left, order, cursor);
            assignRanges(node.right, order, cursor);
        }
        node.rangeEnd = cursor[0];
    }
}
