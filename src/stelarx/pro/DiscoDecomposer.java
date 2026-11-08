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

    private static final long[] NO_MISSING_TAXA = new long[0];

    /** A DISCO component plus taxa lost to cuts below its original root. */
    public static final class ScopedTree extends Tree {
        private final long[] missingWords;
        private final int missingCount;

        private ScopedTree(int index, TreeNode root, int[] leaves, int[] positions,
                           int totalTaxa, long[] missingWords) {
            super(index, root, leaves, positions, leaves.length, totalTaxa);
            this.missingWords = missingWords;
            int count = 0;
            for (long word : missingWords) count += Long.bitCount(word);
            this.missingCount = count;
        }

        public long[] missingWords() { return missingWords; }
        public int missingCount() { return missingCount; }
    }

    /** Build a subset-scoped component after restricting its observed leaves and missing scope. */
    public static ScopedTree restrictedScopedTree(int index, TreeNode root,
                                                   int totalTaxa, long[] missingWords) {
        int leafCount = countLeaves(root);
        int[] leaves = new int[leafCount];
        int[] cursor = {0};
        assignRanges(root, leaves, cursor);
        int[] positions = new int[totalTaxa];
        Arrays.fill(positions, -1);
        for (int i = 0; i < leafCount; i++) {
            if (positions[leaves[i]] >= 0)
                throw new IllegalStateException("Restricted DISCO component duplicated a taxon");
            positions[leaves[i]] = i;
            if (missingWords.length > (leaves[i] >>> 6))
                missingWords[leaves[i] >>> 6] &= ~(1L << leaves[i]);
        }
        return new ScopedTree(index, root, leaves, positions, totalTaxa, missingWords);
    }

    private record CutPlan(IdentityHashMap<TreeNode, Boolean> detachLeft,
                           IdentityHashMap<TreeNode, long[]> originalScopes) {}

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
        return decomposeAll(inputTrees, totalTaxa, firstTreeIndex, minimumLeaves, false);
    }

    /** Scoped mode stores only pre-cut component scopes, not per-node taxon sets. */
    public static Result decomposeAll(List<Tree> inputTrees, int totalTaxa,
                                      int firstTreeIndex, int minimumLeaves,
                                      boolean lineageScoped) {
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
            CutPlan plan = chooseDuplicationCuts(postorder, totalTaxa, lineageScoped);
            IdentityHashMap<TreeNode, Boolean> detachLeft = plan.detachLeft();
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
                // Capture the pre-cut root before unary suppression replaces it.
                long[] scope = lineageScoped
                    ? plan.originalScopes().remove(componentRoot) : null;
                if (lineageScoped && scope == null) {
                    throw new IllegalStateException("Missing pre-cut DISCO component scope");
                }
                TreeNode simplified = suppressUnary(componentRoot, null);
                int leafCount = countLeaves(simplified);
                if (leafCount < minimumLeaves) {
                    discarded++;
                    continue;
                }
                output.add(buildSingleCopyTree(
                    simplified, nextTreeIndex++, totalTaxa, leafCount, scope));
            }
        }

        return new Result(List.copyOf(output), cuts, discarded);
    }

    /** Record each DISCO cut while retaining only the live frontier taxon sets. */
    private static CutPlan chooseDuplicationCuts(
            List<TreeNode> postorder, int totalTaxa, boolean lineageScoped) {
        IdentityHashMap<TreeNode, Set<Integer>> liveTaxa = new IdentityHashMap<>();
        IdentityHashMap<TreeNode, Boolean> detachLeft = new IdentityHashMap<>();
        IdentityHashMap<TreeNode, long[]> scopes = lineageScoped
            ? new IdentityHashMap<>() : null;

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
                boolean cutLeft = left.size() < right.size();
                detachLeft.put(node, cutLeft);
                if (lineageScoped) {
                    // One packed scope per detached component survives merging.
                    scopes.put(cutLeft ? node.left : node.right,
                        packTaxa(cutLeft ? left : right, totalTaxa));
                }
            }
            if (left.size() < right.size()) {
                Set<Integer> swap = left;
                left = right;
                right = swap;
            }
            left.addAll(right);
            liveTaxa.put(node, left);
        }
        if (lineageScoped) {
            TreeNode root = postorder.get(postorder.size() - 1);
            scopes.put(root, packTaxa(liveTaxa.get(root), totalTaxa));
        }
        return new CutPlan(detachLeft, scopes);
    }

    private static long[] packTaxa(Set<Integer> taxa, int totalTaxa) {
        long[] words = new long[(totalTaxa + 63) >>> 6];
        for (int taxon : taxa) words[taxon >>> 6] |= 1L << taxon;
        return words;
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
                                             int totalTaxa, int leafCount,
                                             long[] originalScope) {
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
            if (originalScope != null) originalScope[taxon >>> 6] &= ~(1L << taxon);
        }
        if (originalScope == null) {
            return new Tree(treeIndex, root, postorderArray, positionMap,
                leafCount, totalTaxa);
        }
        // Only taxa absent after cuts within this lineage need completion.
        int used = originalScope.length;
        while (used > 0 && originalScope[used - 1] == 0L) used--;
        long[] missing = used == 0 ? NO_MISSING_TAXA
            : Arrays.copyOf(originalScope, used);
        return new ScopedTree(treeIndex, root, postorderArray, positionMap,
            totalTaxa, missing);
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
