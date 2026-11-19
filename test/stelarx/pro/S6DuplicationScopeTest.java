package stelarx.pro;

import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Scope selection must not alter DISCO cuts, scoring input, or S4's scope. */
public final class S6DuplicationScopeTest {
    private static TreeNode leaf(int taxon) {
        TreeNode n = new TreeNode(); n.taxonId = taxon; return n;
    }

    private static TreeNode node(boolean duplication, TreeNode left, TreeNode right) {
        TreeNode n = new TreeNode();
        n.left = left; n.right = right; left.parent = right.parent = n;
        n.isDuplicationNode = duplication; n.isSpeciationNode = !duplication;
        return n;
    }

    private static void order(TreeNode n, List<Integer> leaves) {
        n.rangeStart = leaves.size();
        if (n.isLeaf()) leaves.add(n.taxonId);
        else { order(n.left, leaves); order(n.right, leaves); }
        n.rangeEnd = leaves.size();
    }

    private static Tree tree(TreeNode root, int taxa) {
        List<Integer> leaves = new ArrayList<>(); order(root, leaves);
        int[] postorder = leaves.stream().mapToInt(Integer::intValue).toArray();
        int[] positions = new int[taxa]; Arrays.fill(positions, -1);
        for (int i = 0; i < postorder.length; i++) positions[postorder[i]] = i;
        return new Tree(0, root, postorder, positions, postorder.length, taxa);
    }

    private static void check(boolean ok, String message) {
        if (!ok) throw new AssertionError(message);
    }

    private static boolean missing(Tree t, int taxon) {
        long[] words = ((DiscoDecomposer.ScopedTree) t).missingWords();
        return (taxon >>> 6) < words.length && (words[taxon >>> 6] & (1L << taxon)) != 0;
    }

    public static void main(String[] args) {
        // The right (0,c) locus is detached on a tie. b is observed only in
        // the sister locus; outside branched off before the duplication.
        for (int b : new int[]{1, 65}) {
            int c = b + 1, outside = b + 2, absent = b + 3;
            Tree gene = tree(node(false,
                node(true, node(false, leaf(0), leaf(b)),
                           node(false, leaf(0), leaf(c))), leaf(outside)), absent + 1);
            List<Tree> strict = DiscoDecomposer.decomposeAll(
                List.of(gene), absent + 1, 0, 1, true).trees();
            List<Tree> dup = DiscoDecomposer.decomposeAll(
                List.of(gene), absent + 1, 0, 1,
                DiscoDecomposer.CompletionScope.DUPLICATION).trees();
            check(strict.size() == 2 && dup.size() == 2, "cut count");
            check(!missing(strict.get(0), b), "S4 strict scope broadened");
            check(missing(dup.get(0), b), "sister-locus species missing from S6 scope");
            check(!missing(dup.get(0), outside), "ancestral sibling entered scope");
            check(!missing(dup.get(0), absent), "unobserved species entered scope");
            check(missing(dup.get(1), c), "residual top scope changed");
            for (int i = 0; i < 2; i++)
                check(Arrays.equals(strict.get(i).postorderArray, dup.get(i).postorderArray),
                    "completion policy changed observed DISCO trees");
            check(gene.leafCount == 5 && gene.root.left.isDuplicationNode,
                "original scoring evidence was modified");
        }
        // A younger duplication must recompute its own observed union rather
        // than inherit species 2/3 from the older duplication's sister locus.
        Tree nested = tree(node(true,
            node(true, node(false, leaf(0), leaf(1)),
                       node(false, leaf(0), leaf(1))),
            node(false, node(false, leaf(0), leaf(1)),
                        node(false, leaf(2), leaf(3)))), 4);
        List<Tree> components = DiscoDecomposer.decomposeAll(
            List.of(nested), 4, 0, 1, DiscoDecomposer.CompletionScope.DUPLICATION).trees();
        check(!missing(components.get(0), 2) && !missing(components.get(0), 3),
            "ancestral scope propagated into younger duplication");
        System.out.println("S6 duplication scope and unchanged S4 scope: PASS");
    }
}
