package stelarx.pro;

import stelarx.cluster.ClusterTable;
import stelarx.dp.DPTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.partition.PartitionTable;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Verifies arbitrary resolution and candidate treatment of inserted nodes. */
public final class GeneTreePolytomyResolverTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected work directory");
        Path work = Path.of(args[0]).toAbsolutePath();
        Files.createDirectories(work);

        Path input = work.resolve("polytomies.tre");
        Path output = work.resolve("resolved.tre");
        Files.writeString(input,
            "(A,B,C,D);\n((A:1,B:2),(D,E,F)90:3,G);\n", StandardCharsets.UTF_8);

        GeneTreePolytomyResolver.Result result =
            GeneTreePolytomyResolver.run(input.toString(), output.toString());
        check(result.treeCount() == 2, "tree count");
        List<String> lines = Files.readAllLines(output, StandardCharsets.UTF_8);
        check(lines.get(0).equals("((A,B):0.0,(C,D):0.0);"),
            "root polytomy resolution");
        check(lines.get(1).equals("(G,((A:1,B:2),(F,(D,E):0.0)90:3):0.0);"),
            "nested polytomy resolution or metadata preservation");

        Path duplicateInput = work.resolve("duplicate-polytomy.tre");
        Path duplicateOutput = work.resolve("duplicate-resolved.tre");
        Files.writeString(duplicateInput, "(A,A,B,C);\n", StandardCharsets.UTF_8);
        GeneTreePolytomyResolver.run(
            duplicateInput.toString(), duplicateOutput.toString());
        check(Files.readString(duplicateOutput, StandardCharsets.UTF_8)
                .equals("((A,A):0.0,(B,C):0.0);\n"),
            "duplicate species labels were not restored after resolution");

        Path binaryInput = work.resolve("binary-input.tre");
        Path binaryOutput = work.resolve("binary-output.tre");
        String binaryTree = "(('A a':1[copy],B:2)95:3,(C:4,D:5):6);\n";
        Files.writeString(binaryInput, binaryTree, StandardCharsets.UTF_8);
        GeneTreePolytomyResolver.run(binaryInput.toString(), binaryOutput.toString());
        check(Files.readString(binaryOutput, StandardCharsets.UTF_8).equals(binaryTree),
            "binary tree was changed");

        Path oneTree = work.resolve("one-polytomy.tre");
        Path oneResolved = work.resolve("one-resolved.tre");
        Files.writeString(oneTree, "(A,B,C,D);\n", StandardCharsets.UTF_8);
        GeneTreePolytomyResolver.run(oneTree.toString(), oneResolved.toString());

        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(oneResolved.toString(), registry, false);
        checkAllInternalSpeciation(trees.get(0).root);

        TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 23L);
        PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);
        UniqueTaxonSubtreeHashes unique = new UniqueTaxonSubtreeHashes(trees, hasher);
        ClusterTable clusters = new ClusterTable(
            trees, pref, registry.size(), unique);
        PartitionTable partitions = new PartitionTable(trees, pref, unique);
        DPTable dp = new DPTable(trees, pref, clusters, unique);
        check(clusters.size() == 6, "inserted speciation cluster count");
        check(partitions.size() == 3, "inserted speciation partition count");
        check(dp.numEmitted() == 3, "inserted speciation transition count");

        System.out.println("STELAR-Pro arbitrary polytomy resolution: PASS");
    }

    private static void checkAllInternalSpeciation(TreeNode node) {
        if (node.isLeaf()) return;
        check(node.isSpeciation(), "resolved unlabeled node was not treated as speciation");
        checkAllInternalSpeciation(node.left);
        checkAllInternalSpeciation(node.right);
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
