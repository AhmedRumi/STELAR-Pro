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

/** Verifies that only biological speciation nodes emit internal candidates. */
public final class SpeciationCandidateFilterTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected work directory");
        Path work = Path.of(args[0]).toAbsolutePath();
        Files.createDirectories(work);

        checkTree(work.resolve("duplication.tre"),
            "(((A,B)D,(C,D)),(E,F));\n", 10, 4, 4);
        // This directly exercises TreeParser's legacy in-memory fallback. Normal
        // Pro inference resolves before tagging, so those serialized refinement
        // nodes receive event tags and are covered by GeneTreePolytomyResolverTest.
        checkTree(work.resolve("untagged-parser-fallback.tre"),
            "((A,B,C),(D,E));\n", 8, 3, 3);

        System.out.println("STELAR-Pro speciation candidate filter: PASS");
    }

    private static void checkTree(Path input, String newick,
                                  int expectedClusters, int expectedPartitions,
                                  int expectedTransitions)
            throws Exception {
        Files.writeString(input, newick, StandardCharsets.UTF_8);
        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(input.toString(), registry, false);
        TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 17L);
        PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);

        UniqueTaxonSubtreeHashes unique = new UniqueTaxonSubtreeHashes(trees, hasher);
        ClusterTable clusters = new ClusterTable(trees, pref, registry.size(), unique);
        check(clusters.size() == expectedClusters, "candidate cluster count");
        verifyNonSpeciationChildSides(trees.get(0).root, clusters);

        PartitionTable partitions = new PartitionTable(trees, pref, unique);
        int partitionOccurrences = partitions.entries().stream()
            .mapToInt(entry -> entry.frequency).sum();
        check(partitions.size() == expectedPartitions, "partition count");
        check(partitionOccurrences == expectedPartitions, "partition occurrence count");

        DPTable dp = new DPTable(trees, pref, clusters, unique);
        check(dp.numEmitted() == expectedTransitions, "transition emission count");
        check(dp.numUniqueSplits() == expectedTransitions, "unique transition count");
    }

    private static void verifyNonSpeciationChildSides(
            TreeNode node, ClusterTable clusters) {
        if (node.isLeaf()) return;
        TreeNode[] children = node.isPolytomous()
            ? node.children : new TreeNode[]{node.left, node.right};
        if (node.isSpeciation()) {
            for (TreeNode child : children) {
                if (child.isLeaf() || child.isSpeciation()) continue;
                boolean present = clusters.entries().stream().anyMatch(entry ->
                    entry.exemplar.treeIndex == 0
                        && entry.exemplar.left == child.rangeStart
                        && entry.exemplar.right == child.rangeEnd);
                check(present,
                    "non-speciation child set of a speciation split was unavailable");
            }
        }
        for (TreeNode child : children) verifyNonSpeciationChildSides(child, clusters);
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
