package stelarx.pro;

import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;

/**
 * Decomposes rooted, ASTRAL-Pro-tagged multicopy trees into single-copy trees.
 *
 * <p>This is the DISCO postorder rule used by {@code scripts/disco.py}: at each
 * duplication node, detach the child containing fewer distinct species (the
 * right child on a tie), retain every detached component, and suppress the
 * resulting unary duplication nodes. The input trees are copied first because
 * STELAR-Pro must keep them unchanged for its rooted-triplet scoring model.</p>
 */
public final class DiscoDecomposer {
    private DiscoDecomposer() {}

    /** Summary of one decomposition pass. */
    public record Result(List<Tree> trees, int duplicationCuts,
                         int discardedSmallTrees) {}

    /**
     * Decompose all input trees and assign output indices consecutively.
     * Components smaller than {@code minimumLeaves} are omitted, matching
     * DISCO's default filtering of trees that cannot display a quartet.
     */
    public static Result decomposeAll(List<Tree> inputTrees, int totalTaxa,
                                      int firstTreeIndex, int minimumLeaves) {
        if (minimumLeaves < 1) {
            throw new IllegalArgumentException("DISCO minimum leaf count must be positive");
        }

        List<Tree> output = new ArrayList<>();
        int nextTreeIndex = firstTreeIndex;
        int cuts = 0;
        int discarded = 0;

        for (Tree input : inputTrees) {
            if (input.hasPolytomy) {
                throw new IllegalArgumentException(
                    "DISCO requires binary rooted/tagged trees; unresolved tree "
                        + input.treeIndex + " was supplied");
            }

            TreeNode root = copyTree(input.root, null);
            List<TreeNode> postorder = binaryPostorder(root);
            IdentityHashMap<TreeNode, Boolean> detachLeft =
                chooseDuplicationCuts(postorder);
            List<TreeNode> componentRoots = new ArrayList<>(detachLeft.size() + 1);

            // Apply cuts in postorder, as in DISCO. The saved decisions use the
            // original child species sets, before any descendant edge is cut.
            for (TreeNode node : postorder) {
                Boolean chooseLeft = detachLeft.get(node);
                if (chooseLeft == null) continue;
                TreeNode detached = chooseLeft ? node.left : node.right;
                if (chooseLeft) node.left = null;
                else            node.right = null;
                detached.parent = null;
                componentRoots.add(detached);
                cuts++;
            }
            componentRoots.add(root);

            for (TreeNode componentRoot : componentRoots) {
                TreeNode simplified = suppressUnary(componentRoot, null);
                int leafCount = countLeaves(simplified);
                if (leafCount < minimumLeaves) {
                    discarded++;
                    continue;
                }
                output.add(buildSingleCopyTree(
                    simplified, nextTreeIndex++, totalTaxa, leafCount));
            }
        }

        return new Result(List.copyOf(output), cuts, discarded);
    }

    /** Record each DISCO cut while retaining only the live frontier taxon sets. */
    private static IdentityHashMap<TreeNode, Boolean> chooseDuplicationCuts(
            List<TreeNode> postorder) {
        IdentityHashMap<TreeNode, Set<Integer>> liveTaxa = new IdentityHashMap<>();
        IdentityHashMap<TreeNode, Boolean> detachLeft = new IdentityHashMap<>();

        for (TreeNode node : postorder) {
            if (node.isLeaf()) {
                Set<Integer> singleton = new HashSet<>(1);
                singleton.add(node.taxonId);
                liveTaxa.put(node, singleton);
                continue;
            }

            Set<Integer> left = liveTaxa.remove(node.left);
            Set<Integer> right = liveTaxa.remove(node.right);
            if (left == null || right == null) {
                throw new IllegalStateException("Incomplete DISCO postorder taxon sets");
            }
            if (node.isDuplicationNode) {
                // scripts/disco.py uses '<'; an equal-size tie therefore cuts right.
                detachLeft.put(node, left.size() < right.size());
            }
            if (left.size() < right.size()) {
                Set<Integer> swap = left;
                left = right;
                right = swap;
            }
            left.addAll(right);
            liveTaxa.put(node, left);
        }
        return detachLeft;
    }

    private static List<TreeNode> binaryPostorder(TreeNode root) {
        ArrayDeque<TreeNode> pending = new ArrayDeque<>();
        ArrayDeque<TreeNode> reverse = new ArrayDeque<>();
        pending.push(root);
        while (!pending.isEmpty()) {
            TreeNode node = pending.pop();
            reverse.push(node);
            if (node.isLeaf()) continue;
            if (node.isPolytomous() || node.left == null || node.right == null) {
                throw new IllegalArgumentException("DISCO received a non-binary tree");
            }
            pending.push(node.left);
            pending.push(node.right);
        }
        return new ArrayList<>(reverse);
    }

    /** Remove unary duplication nodes left behind by the selected edge cuts. */
    private static TreeNode suppressUnary(TreeNode node, TreeNode parent) {
        if (node == null) return null;
        TreeNode left = node.left == null ? null : suppressUnary(node.left, node);
        TreeNode right = node.right == null ? null : suppressUnary(node.right, node);
        if (left == null && right == null) {
            node.parent = parent;
            return node;
        }
        if (left == null || right == null) {
            TreeNode child = left != null ? left : right;
            child.parent = parent;
            return child;
        }
        node.left = left;
        node.right = right;
        node.children = null;
        node.parent = parent;
        left.parent = node;
        right.parent = node;
        return node;
    }

    private static Tree buildSingleCopyTree(TreeNode root, int treeIndex,
                                             int totalTaxa, int leafCount) {
        int[] postorderArray = new int[leafCount];
        int[] cursor = {0};
        assignRanges(root, postorderArray, cursor);

        int[] positionMap = new int[totalTaxa];
        Arrays.fill(positionMap, -1);
        for (int position = 0; position < postorderArray.length; position++) {
            int taxon = postorderArray[position];
            if (positionMap[taxon] != -1) {
                throw new IllegalStateException(
                    "DISCO component still contains duplicate taxon " + taxon);
            }
            positionMap[taxon] = position;
        }
        return new Tree(treeIndex, root, postorderArray, positionMap,
            leafCount, totalTaxa);
    }

    private static void assignRanges(TreeNode node, int[] postorderArray,
                                     int[] cursor) {
        node.rangeStart = cursor[0];
        if (node.isLeaf()) {
            postorderArray[cursor[0]++] = node.taxonId;
        } else {
            assignRanges(node.left, postorderArray, cursor);
            assignRanges(node.right, postorderArray, cursor);
        }
        node.rangeEnd = cursor[0];
    }

    private static int countLeaves(TreeNode root) {
        int count = 0;
        ArrayDeque<TreeNode> pending = new ArrayDeque<>();
        pending.push(root);
        while (!pending.isEmpty()) {
            TreeNode node = pending.pop();
            if (node.isLeaf()) count++;
            else {
                pending.push(node.left);
                pending.push(node.right);
            }
        }
        return count;
    }

    private static TreeNode copyTree(TreeNode source, TreeNode parent) {
        TreeNode copy = new TreeNode();
        copy.taxonId = source.taxonId;
        copy.isDuplicationNode = source.isDuplicationNode;
        copy.isSpeciationNode = source.isSpeciationNode;
        copy.parent = parent;
        if (!source.isLeaf()) {
            if (source.isPolytomous()) {
                throw new IllegalArgumentException("DISCO received a non-binary tree");
            }
            copy.left = copyTree(source.left, copy);
            copy.right = copyTree(source.right, copy);
        }
        return copy;
    }
}
