package stelarx.completion;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.List;

/** Regression coverage for completion of incomplete, natively polytomous trees. */
public final class TreeCompleterPolytomyTest {
    /**
     * Root-preserving nearest-neighbour completion of the fixture (the input
     * root and every existing relationship are kept; each missing taxon joins
     * its most similar present taxon).  Regenerated 2026-11-05 when the fixture
     * was given explicitly binary roots, which the parser now requires.
     */
    private static final String[] EXPECTED = {
        "((A,B,C),((D,E),(F,G)));",
        "(A,((B,C,D),((F,E),G)));",
        "((A,C,E),(B,(D,(G,F))));",
        "(((A,C),D),((B,E,G),F));"
    };

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass the incomplete-polytomy fixture");
        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        try {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> trees = TreeParser.parseGeneTrees(args[0], registry, true);
            for (Tree tree : trees) {
                check(tree.hasPolytomy, "input tree " + tree.treeIndex + " lost its polytomy marker");
                check(!tree.isComplete, "fixture tree " + tree.treeIndex + " must be incomplete");
            }

            int n = registry.size();
            System.clearProperty("stelarpro.similarity.segmentShift");
            SimilarityMatrix matrix = SimilarityMatrixBuilder.buildCPU(trees, n);
            List<Tree> completed = TreeCompleter.completeAll(trees, matrix, n);

            // Tiny segments must not change anchors.
            System.setProperty("stelarpro.similarity.segmentShift", "2");
            SimilarityMatrix segmentedMatrix;
            try {
                segmentedMatrix = SimilarityMatrixBuilder.buildCPU(trees, n);
            } finally {
                System.clearProperty("stelarpro.similarity.segmentShift");
            }
            List<Tree> segmented = TreeCompleter.completeAll(trees, segmentedMatrix, n);

            for (int i = 0; i < trees.size(); i++) {
                validateCompleted(completed.get(i), n);
                validateCompleted(segmented.get(i), n);
                String newick = completed.get(i).toNewick(registry);
                String segmentedNewick = segmented.get(i).toNewick(registry);
                check(newick.equals(segmentedNewick),
                    "default/segmented completion mismatch for tree " + i);
                check(EXPECTED[i].equals(newick),
                    "completion topology changed for tree " + i + ": " + newick);
            }
            System.out.println("Incomplete polytomy completion: PASS");
        } finally {
            Threading.shutdown();
            System.clearProperty("stelarpro.similarity.segmentShift");
        }
    }

    private static void validateCompleted(Tree tree, int n) {
        check(tree.isComplete, "completed tree " + tree.treeIndex + " is still incomplete");
        check(tree.leafCount == n, "completed tree " + tree.treeIndex + " has wrong leaf count");
        check(tree.hasPolytomy, "completed tree " + tree.treeIndex + " lost all polytomies");
        boolean[] seen = new boolean[n];
        int[] position = {0};
        int[] polytomies = {0};
        validateNode(tree.root, null, seen, position, polytomies);
        check(position[0] == n, "tree traversal did not reach every taxon");
        check(polytomies[0] > 0, "tree contains no canonical polytomy node");
        for (int taxon = 0; taxon < n; taxon++) {
            check(seen[taxon], "missing taxon " + taxon);
            check(tree.positionMap[taxon] >= 0, "position map missing taxon " + taxon);
        }
    }

    private static void validateNode(TreeNode node, TreeNode parent, boolean[] seen,
                                     int[] position, int[] polytomies) {
        check(node.parent == parent, "broken parent pointer");
        if (node.isLeaf()) {
            check(node.taxonId >= 0 && node.taxonId < seen.length, "invalid leaf taxon ID");
            check(!seen[node.taxonId], "duplicate leaf taxon " + node.taxonId);
            seen[node.taxonId] = true;
            check(node.rangeStart == position[0] && node.rangeEnd == position[0] + 1,
                "invalid leaf range");
            position[0]++;
            return;
        }

        int start = position[0];
        if (node.isPolytomous()) {
            polytomies[0]++;
            check(node.children.length >= 3, "polytomy has fewer than three children");
            check(node.left == node.children[0], "left alias is not first polytomy child");
            check(node.right == node.children[node.children.length - 1],
                "right alias is not last polytomy child");
            for (TreeNode child : node.children) validateNode(child, node, seen, position, polytomies);
        } else {
            check(node.left != null && node.right != null, "binary node is missing a child");
            validateNode(node.left, node, seen, position, polytomies);
            validateNode(node.right, node, seen, position, polytomies);
        }
        check(node.rangeStart == start && node.rangeEnd == position[0], "invalid internal range");
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
