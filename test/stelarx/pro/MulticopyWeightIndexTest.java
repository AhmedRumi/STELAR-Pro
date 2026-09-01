package stelarx.pro;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.weight.IntersectionCounter;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Verifies duplicate-invariant range intersections through position vectors. */
public final class MulticopyWeightIndexTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected work directory");
        Path work = Path.of(args[0]).toAbsolutePath();
        Files.createDirectories(work);
        Path input = work.resolve("multicopy-index.tre");
        Files.writeString(input,
            "(((A,A)D,B),(C,(A,D))D);\n"
          + "((D,A),(B,(A,C))D);\n"
          + "((A,A)D,(B,C));\n",
            StandardCharsets.UTF_8);

        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(input.toString(), registry, false);
        Tree first = trees.get(0);
        Tree second = trees.get(1);
        Tree incomplete = trees.get(2);
        int a = registry.getId("A");

        check(first.leafCount == 6, "leaf-copy count");
        check(first.distinctTaxonCount == 4, "distinct species count");
        check(first.isComplete, "multicopy complete tree was marked incomplete");
        check(!incomplete.isComplete,
            "duplicate copies masked a missing species in completeness test");

        check(first.taxonPositions.countInRange(a, 0, 6) == 3,
            "all A-copy positions were not retained");
        check(first.taxonPositions.countInRange(a, 0, 3) == 2,
            "A-copy subrange count");
        check(first.taxonPositions.containsInRange(a, 4, 5),
            "binary-search membership missed a later copy");

        check(IntersectionCounter.coreIntersect(first, 0, 3, second, 1, 5) == 2,
            "duplicate-invariant core intersection");
        check(IntersectionCounter.coreIntersect(first, 0, 2, second, 0, 4) == 1,
            "duplicate source copies were counted repeatedly");
        check(IntersectionCounter.intersect(first, 0, 3, second, 0, 2,
                true, 2) == 1,
            "complement intersection did not use distinct range size");
        check(IntersectionCounter.intersectWithFullTree(
                first, second, 1, 4, true) == 2,
            "full-tree complement did not use distinct species count");
        check(IntersectionCounter.coreIntersectMulti(first, 0, 6, second,
                new int[]{0, 3}, new int[]{2, 5}) == 3,
            "taxon shared by multiple target ranges was counted more than once");

        System.out.println("STELAR-Pro multicopy weight index: PASS");
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
