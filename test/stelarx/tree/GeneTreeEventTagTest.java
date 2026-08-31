package stelarx.tree;

import stelarx.taxon.TaxonRegistry;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Verifies ASTRAL-Pro D labels and unlabeled speciation nodes survive parsing. */
public final class GeneTreeEventTagTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected work directory");
        Path work = Path.of(args[0]).toAbsolutePath();
        Files.createDirectories(work);
        Path input = work.resolve("tagged-multicopy.tre");
        Files.writeString(input,
            "((A,A)D,(B,C));\n((A,B,C),(D,E));\n", StandardCharsets.UTF_8);

        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(input.toString(), registry, false);
        check(trees.size() == 2, "tree count");

        Tree first = trees.get(0);
        check(first.leafCount == 4, "multi-copy occurrence array");
        check(first.root.isSpeciation(), "unlabeled root must be speciation");
        check(first.root.left.isDuplication(), "D node must be duplication");
        check(first.root.right.isSpeciation(), "unlabeled internal must be speciation");

        Tree second = trees.get(1);
        check(second.root.left.isSpeciation(), "original polytomy event lost");
        TreeNode artificial = second.root.left.right;
        check(!artificial.isSpeciation() && !artificial.isDuplication(),
            "artificial binary refinement was tagged as biological");

        System.out.println("STELAR-Pro gene-tree event tags: PASS");
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
