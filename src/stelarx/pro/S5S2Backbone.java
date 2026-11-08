package stelarx.pro;

import stelarx.Logging;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Induce strict-lineage DISCO/UPGMA proposals on an S5 backbone subset. */
final class S5S2Backbone {
    private S5S2Backbone() {}

    static List<Tree> restrict(List<Tree> completed, Tree globalUpgma,
                               String[] globalNames, TaxonRegistry subset) {
        Map<String, Integer> globalIds = new HashMap<>(globalNames.length * 2);
        for (int i = 0; i < globalNames.length; i++) globalIds.put(globalNames[i], i);
        int[] toSubset = new int[globalNames.length];
        Arrays.fill(toSubset, -1);
        for (int i = 0; i < subset.size(); i++) {
            Integer global = globalIds.get(subset.getName(i));
            if (global == null) throw new IllegalArgumentException(
                "Backbone taxon missing from the global S2 matrix: " + subset.getName(i));
            toSubset[global] = i;
        }

        List<Tree> proposals = new ArrayList<>(completed.size() + 1);
        for (Tree source : completed) {
            Tree induced = restrictOne(source, toSubset, proposals.size(), subset.size());
            // One- and two-taxon trees cannot provide a rooted triplet split.
            if (induced != null && induced.leafCount >= 3) proposals.add(induced);
        }
        Tree upgma = restrictOne(globalUpgma, toSubset, proposals.size(), subset.size());
        if (upgma == null || !upgma.isComplete) {
            throw new IllegalStateException("Restricted UPGMA does not cover the backbone");
        }
        proposals.add(upgma);
        Logging.info("S5/S2 backbone subset: %d taxa, %d DISCO proposals + 1 UPGMA",
            subset.size(), proposals.size() - 1);
        return proposals;
    }

    private static Tree restrictOne(Tree source, int[] toSubset, int index, int n) {
        TreeNode root = prune(source.root, toSubset);
        if (root == null) return null;
        int[] order = new int[n];
        int[] count = {0};
        assignRanges(root, order, count);
        int[] leaves = Arrays.copyOf(order, count[0]);
        int[] position = new int[n];
        Arrays.fill(position, -1);
        for (int i = 0; i < leaves.length; i++) {
            if (position[leaves[i]] != -1)
                throw new IllegalStateException("Restricted DISCO tree duplicated a taxon");
            position[leaves[i]] = i;
        }
        return new Tree(index, root, leaves, position, leaves.length, n);
    }

    private static TreeNode prune(TreeNode source, int[] toSubset) {
        if (source.isLeaf()) {
            int mapped = toSubset[source.taxonId];
            if (mapped < 0) return null;
            TreeNode leaf = new TreeNode();
            leaf.taxonId = mapped;
            return leaf;
        }
        if (source.isPolytomous())
            throw new IllegalStateException("S5/S2 backbone requires binary DISCO trees");
        TreeNode left = prune(source.left, toSubset);
        TreeNode right = prune(source.right, toSubset);
        if (left == null) return right;
        if (right == null) return left;
        TreeNode parent = new TreeNode();
        parent.left = left;
        parent.right = right;
        parent.isSpeciationNode = true; // completed S2 nodes are proposal nodes
        left.parent = parent;
        right.parent = parent;
        return parent;
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
