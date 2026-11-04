package stelarx.pro;

import stelarx.completion.SimilarityMatrix;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.Threading;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Checks that scoped DISCO completion never imports ancestral siblings. */
public final class LineageScopeTest {
    private static TreeNode leaf(int taxon) {
        TreeNode node = new TreeNode();
        node.taxonId = taxon;
        return node;
    }

    private static TreeNode node(boolean duplication, TreeNode left, TreeNode right) {
        TreeNode node = new TreeNode();
        node.isDuplicationNode = duplication;
        node.isSpeciationNode = !duplication;
        node.left = left;
        node.right = right;
        left.parent = right.parent = node;
        return node;
    }

    private static void order(TreeNode node, List<Integer> leaves) {
        node.rangeStart = leaves.size();
        if (node.isLeaf()) leaves.add(node.taxonId);
        else {
            order(node.left, leaves);
            order(node.right, leaves);
        }
        node.rangeEnd = leaves.size();
    }

    private static Tree tree(TreeNode root, int totalTaxa) {
        List<Integer> leaves = new ArrayList<>();
        order(root, leaves);
        int[] order = leaves.stream().mapToInt(Integer::intValue).toArray();
        int[] positions = new int[totalTaxa];
        Arrays.fill(positions, -1);
        for (int i = 0; i < order.length; i++) positions[order[i]] = i;
        return new Tree(0, root, order, positions, order.length, totalTaxa);
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }

    public static void main(String[] args) {
        Threading.start(2);
        try {
            Tree noLoss = tree(node(false,
                node(true, node(false, leaf(0), leaf(1)),
                           node(false, leaf(0), leaf(1))),
                node(false, leaf(2), leaf(3))), 4);
            DiscoDecomposer.Result first = DiscoDecomposer.decomposeAll(
                List.of(noLoss), 4, 0, 1, true);
            check(first.trees().size() == 2, "no-loss component count");
            for (Tree candidate : first.trees()) {
                DiscoDecomposer.ScopedTree scoped =
                    (DiscoDecomposer.ScopedTree) candidate;
                check(scoped.missingCount() == 0,
                    "no-loss component must span its original lineage");
                if (scoped.leafCount == 2) {
                    check(scoped.positionMap[2] == -1 && scoped.positionMap[3] == -1,
                        "ancestral siblings entered detached lineage");
                }
            }

            Tree nested = tree(node(false,
                node(true, node(false, leaf(0), leaf(1)),
                           node(false, leaf(2), leaf(3))),
                node(false, leaf(4), leaf(5))), 6);
            DiscoDecomposer.Result second = DiscoDecomposer.decomposeAll(
                List.of(nested), 6, 0, 1, true);
            check(second.trees().size() == 2, "nested-cut component count");
            Tree retained = second.trees().stream().filter(t -> t.leafCount == 4)
                .findFirst().orElseThrow();
            Tree detached = second.trees().stream().filter(t -> t.leafCount == 2)
                .findFirst().orElseThrow();
            check(((DiscoDecomposer.ScopedTree) retained).missingCount() == 2,
                "retained lineage must recover taxa lost below its root");
            check(((DiscoDecomposer.ScopedTree) detached).missingCount() == 0,
                "detached lineage must exclude taxa outside its root");

            SimilarityMatrix guide = new SimilarityMatrix(6);
            guide.normalize();
            List<Tree> completed = ScopedTreeCompleter.prepare(
                guide, second.trees(), 6).completeBatch(second.trees());
            check(completed.get(second.trees().indexOf(retained)).leafCount == 6,
                "retained lineage completion failed");
            check(completed.get(second.trees().indexOf(detached)).leafCount == 2,
                "outside taxa were added to detached lineage");

            DiscoDecomposer.Result legacy = DiscoDecomposer.decomposeAll(
                List.of(nested), 6, 0, 1, false);
            check(legacy.trees().stream().noneMatch(t ->
                t instanceof DiscoDecomposer.ScopedTree),
                "global mode should not retain scope bitsets");
            System.out.println("Lineage-scoped DISCO tests passed");
        } finally {
            Threading.shutdown();
        }
    }
}
