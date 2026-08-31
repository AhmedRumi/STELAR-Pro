package stelarx.pro;

import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.dp.BipartitionSplit;
import stelarx.dp.DPTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Checks duplicate-invariant S1 cluster and split deduplication. */
public final class DuplicateAwareCandidateTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected work directory");
        Path work = Path.of(args[0]).toAbsolutePath();
        Files.createDirectories(work);
        Path input = work.resolve("tagged-multicopy.tre");
        Files.writeString(input,
            "((((A,A)D,B),(C,D)),(E,F));\n(((A,B),(C,D)),(E,F));\n",
            StandardCharsets.UTF_8);

        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(input.toString(), registry, false);
        TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 31L);
        PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);

        ClusterTable oldClusters = new ClusterTable(trees, pref, registry.size());
        DPTable oldDP = new DPTable(trees, pref, oldClusters);

        UniqueTaxonSubtreeHashes unique = new UniqueTaxonSubtreeHashes(trees, hasher);
        ClusterTable newClusters = new ClusterTable(trees, pref, registry.size(), unique);
        DPTable newDP = new DPTable(trees, pref, newClusters, unique);

        check(oldClusters.size() == 12, "occurrence-sensitive cluster baseline");
        check(newClusters.size() == 10, "duplicate-invariant cluster count");
        check(oldDP.numUniqueSplits() == 8, "occurrence-sensitive split baseline");
        check(newDP.numUniqueSplits() == 5, "duplicate-invariant split count");
        check(newDP.numOverlappingSpeciationNodesSkipped() == 0,
            "correctly tagged speciation was rejected");
        check(newDP.getRootHash().size == registry.size(), "root distinct-taxon size");

        Tree first = trees.get(0);
        TreeNode duplication = first.root.left.left.left;
        ClusterHash singletonA = unique.get(0, duplication);
        check(singletonA.size == 1, "duplicate subtree counted A more than once");
        check(newClusters.get(singletonA).frequency == 3,
            "duplication-rooted subtree was added as a candidate");

        TreeNode aabSpeciation = first.root.left.left;
        ClusterHash ab = unique.get(0, aabSpeciation);
        check(ab.size == 2, "speciation subtree did not deduplicate A");
        check(newClusters.get(ab).exemplar.size == 2,
            "cluster exemplar retained occurrence-sensitive size");

        for (var entry : newDP.entries()) {
            for (BipartitionSplit split : entry.getValue()) {
                check(split.lo.size + split.hi.size == entry.getKey().size,
                    "DP split contains overlapping species");
            }
        }

        System.out.println("STELAR-Pro duplicate-aware S1 candidates: PASS");
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
