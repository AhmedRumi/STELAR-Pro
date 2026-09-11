package stelarx.pro;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/** Verifies DISCO cut selection, single-copy output, and input immutability. */
public final class DiscoDecomposerTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected work directory");
        Path work = Path.of(args[0]).toAbsolutePath();
        Files.createDirectories(work);
        Path input = work.resolve("tagged-multicopy.tre");

        // Both children of D contain four species, so scripts/disco.py's tie
        // rule detaches the right child. The outer root remains a speciation.
        Files.writeString(input,
            "((((A,B),(C,D)),((A,E),(F,G)))D,(H,I));\n",
            StandardCharsets.UTF_8);
        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> originals = TreeParser.parseGeneTrees(
            input.toString(), registry, false);
        Tree original = originals.get(0);

        DiscoDecomposer.Result result = DiscoDecomposer.decomposeAll(
            originals, registry.size(), originals.size(), 4);
        check(result.duplicationCuts() == 1, "duplication cut count");
        check(result.discardedSmallTrees() == 0, "discard count");
        check(result.trees().size() == 2, "component count");

        // Detached components precede the retained main component.
        check(taxa(result.trees().get(0), registry).equals(Set.of("A", "E", "F", "G")),
            "equal-size tie must detach the right child");
        check(taxa(result.trees().get(1), registry).equals(
                Set.of("A", "B", "C", "D", "H", "I")),
            "retained main component");

        for (int i = 0; i < result.trees().size(); i++) {
            Tree tree = result.trees().get(i);
            check(tree.treeIndex == originals.size() + i, "contiguous tree index");
            check(tree.leafCount == tree.distinctTaxonCount, "single-copy component");
            check(!containsDuplication(tree.root), "duplication node survived suppression");
        }

        // Decomposition operates on a copy; the GDL scoring tree must not change.
        check(original.leafCount == 10, "original leaf occurrences changed");
        check(original.root.left.isDuplication(), "original duplication tag/topology changed");

        System.out.println("STELAR-Pro DISCO decomposition: PASS");
    }

    private static Set<String> taxa(Tree tree, TaxonRegistry registry) {
        Set<String> names = new HashSet<>();
        for (int taxon : tree.postorderArray) names.add(registry.getName(taxon));
        return names;
    }

    private static boolean containsDuplication(TreeNode node) {
        if (node.isLeaf()) return false;
        return node.isDuplicationNode
            || containsDuplication(node.left)
            || containsDuplication(node.right);
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError("failed: " + message);
    }
}
