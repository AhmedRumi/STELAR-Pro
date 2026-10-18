package stelarx.tree;

import stelarx.taxon.TaxonRegistry;

/**
 * A rooted gene tree with postorder leaves and taxon indexes. Streamed S2
 * candidate exemplars retain these indexes but intentionally have no node graph.
 */
public class Tree {
    /** Index in the gene tree list (0..k-1). */
    public final int treeIndex;

    /** Root of the tree; null only for a compact candidate-membership exemplar. */
    public final TreeNode root;

    /**
     * postorderArray[pos] = taxon ID at left-to-right position pos.
     * Length = leafCount.
     */
    public final int[] postorderArray;

    /**
     * Legacy representative position per taxon (-1 if absent). Multicopy code
     * must use {@link #taxonPositions}, which retains every position.
     */
    public final int[] positionMap;

    /** All leaf-copy positions, grouped by taxon in sorted CSR rows. */
    public final TaxonPositionIndex taxonPositions;

    /** Number of leaves in this tree. */
    public final int leafCount;

    /** Number of distinct species represented by the leaves. */
    public final int distinctTaxonCount;

    /** True when this tree contains all n taxa. */
    public final boolean isComplete;

    /** True when at least one internal node has three or more rooted children. */
    public final boolean hasPolytomy;

    public Tree(int treeIndex, TreeNode root,
                int[] postorderArray, int[] positionMap,
                int leafCount, int totalTaxa) {
        this(treeIndex, root, postorderArray, positionMap, leafCount, totalTaxa, false);
    }

    public Tree(int treeIndex, TreeNode root,
                int[] postorderArray, int[] positionMap,
                int leafCount, int totalTaxa, boolean hasPolytomy) {
        this.treeIndex = treeIndex;
        this.root = root;
        this.postorderArray = postorderArray;
        this.positionMap = positionMap;
        this.leafCount = leafCount;
        this.taxonPositions = TaxonPositionIndex.build(postorderArray, totalTaxa);
        this.distinctTaxonCount = taxonPositions.distinctTaxonCount();
        this.isComplete = (distinctTaxonCount == totalTaxa);
        this.hasPolytomy = hasPolytomy;
    }

    /** Keep only the leaf indexes needed to interpret a streamed candidate exemplar. */
    public static Tree candidateExemplar(int index, Tree completed, int totalTaxa) {
        if (completed.positionMap.length != totalTaxa) {
            throw new IllegalArgumentException("Candidate exemplar taxon count mismatch");
        }
        return new Tree(index, completed);
    }

    private Tree(int index, Tree completed) {
        this.treeIndex = index;
        this.root = null;
        this.postorderArray = completed.postorderArray;
        this.positionMap = completed.positionMap;
        this.taxonPositions = completed.taxonPositions;
        this.leafCount = completed.leafCount;
        this.distinctTaxonCount = completed.distinctTaxonCount;
        this.isComplete = completed.isComplete;
        this.hasPolytomy = completed.hasPolytomy;
    }

    /** Reconstruct Newick string (no branch lengths). */
    public String toNewick(TaxonRegistry reg) {
        if (root == null) {
            throw new IllegalStateException("Compact candidate exemplar has no tree topology");
        }
        return (hasPolytomy ? nodeToNewickPolytomy(root, reg) : nodeToNewick(root, reg)) + ";";
    }

    private String nodeToNewick(TreeNode n, TaxonRegistry reg) {
        if (n.isLeaf()) return reg.getName(n.taxonId);
        return "(" + nodeToNewick(n.left, reg) + "," + nodeToNewick(n.right, reg) + ")";
    }

    private String nodeToNewickPolytomy(TreeNode n, TaxonRegistry reg) {
        if (n.isLeaf()) return reg.getName(n.taxonId);
        if (!n.isPolytomous()) {
            return "(" + nodeToNewickPolytomy(n.left, reg) + ","
                + nodeToNewickPolytomy(n.right, reg) + ")";
        }
        StringBuilder out = new StringBuilder("(");
        for (int i = 0; i < n.children.length; i++) {
            if (i > 0) out.append(',');
            out.append(nodeToNewickPolytomy(n.children[i], reg));
        }
        return out.append(')').toString();
    }
}
