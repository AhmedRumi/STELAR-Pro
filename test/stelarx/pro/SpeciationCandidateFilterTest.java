package stelarx.pro;

import stelarx.cluster.ClusterTable;
import stelarx.dp.DPTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.partition.PartitionTable;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
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
            "(((A,B)D,(C,D)),(E,F));\n", 9, 4, 4, 0, 2);
        checkTree(work.resolve("artificial.tre"),
            "((A,B,C),(D,E));\n", 7, 3, 3, 0, 2);

        System.out.println("STELAR-Pro speciation candidate filter: PASS");
    }

    private static void checkTree(Path input, String newick,
                                  int expectedClusters, int expectedPartitions,
                                  int expectedTransitions, int skippedStart, int skippedEnd)
            throws Exception {
        Files.writeString(input, newick, StandardCharsets.UTF_8);
        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(input.toString(), registry, false);
        TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 17L);
        PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);

        ClusterTable clusters = new ClusterTable(trees, pref, registry.size());
        check(clusters.size() == expectedClusters, "candidate cluster count");
        boolean skippedClusterPresent = clusters.entries().stream().anyMatch(entry ->
            !entry.exemplar.complement
                && entry.exemplar.left == skippedStart
                && entry.exemplar.right == skippedEnd);
        check(!skippedClusterPresent, "non-speciation cluster was retained");

        PartitionTable partitions = new PartitionTable(trees, pref);
        int partitionOccurrences = partitions.entries().stream()
            .mapToInt(entry -> entry.frequency).sum();
        check(partitions.size() == expectedPartitions, "partition count");
        check(partitionOccurrences == expectedPartitions, "partition occurrence count");

        DPTable dp = new DPTable(trees, pref, clusters);
        check(dp.numEmitted() == expectedTransitions, "transition emission count");
        check(dp.numUniqueSplits() == expectedTransitions, "unique transition count");
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
