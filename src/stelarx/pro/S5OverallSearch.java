package stelarx.pro;

import stelarx.Config;
import stelarx.Logging;
import stelarx.cluster.ClusterTable;
import stelarx.dp.DPTable;
import stelarx.dp.Inference;
import stelarx.gpu.GPUDPBuilder;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.partition.PartitionTable;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;
import stelarx.weight.WeightTable;

import java.util.List;

/** Re-runs rooted-triplet DP over accumulated greedy and UPGMA proposals. */
public final class S5OverallSearch {
    private S5OverallSearch() {}

    public record Result(String newick, String tripletScore) {}

    public static Result infer(List<Tree> originalTrees, List<Tree> candidates,
                               TaxonRegistry registry) {
        Config cfg = Config.getInstance();
        TaxonHasher hasher = new TaxonHasher(registry.size(),
            cfg.getNumHashSeeds(), cfg.getBaseSeed());
        PrefixHashArrays candidatePrefix = new PrefixHashArrays(candidates, hasher);
        PrefixHashArrays genePrefix = new PrefixHashArrays(originalTrees, hasher);
        UniqueTaxonSubtreeHashes hashes = new UniqueTaxonSubtreeHashes(candidates, hasher);
        UniqueTaxonSubtreeHashes geneHashes =
            new UniqueTaxonSubtreeHashes(originalTrees, hasher);
        int greedyCount = candidates.size() - 1; // UPGMA remains last
        ClusterTable clusters = new ClusterTable(candidates, candidatePrefix,
            registry.size(), hashes, 0, greedyCount);
        PartitionTable partitions = new PartitionTable(originalTrees, genePrefix, geneHashes);
        DPTable dp = new DPTable(candidates, candidatePrefix, clusters, hashes,
            0, greedyCount);
        hashes.release();
        geneHashes.release();
        if (cfg.getSearchMode() == Config.SearchMode.FULL)
            dp.addCrossTreeTransitions(clusters,
                cfg.getComputeMode() == Config.ComputeMode.GPU && GPUDPBuilder.tryLoad());
        WeightTable weights = new WeightTable(dp, partitions, clusters,
            candidates, originalTrees);
        Inference inference = new Inference();
        String tree = inference.run(dp, weights, clusters, candidates, registry, hasher);
        Logging.info("S5 accumulated DP: %d greedy trees + 1 UPGMA guide; "
            + "%d original genes supply triplet scores only; score %s",
            greedyCount, originalTrees.size(), inference.getLastTripletScore());
        return new Result(tree, inference.getLastTripletScore());
    }

    /** Compare rooted topology, ignoring Newick child order and branch lengths. */
    public static boolean sameRootedTopology(String first, String second,
                                             TaxonRegistry registry) {
        Tree a = TreeParser.parseSpeciesTreeNewick(first, registry);
        Tree b = TreeParser.parseSpeciesTreeNewick(second, registry);
        return signature(a.root).equals(signature(b.root));
    }

    private static String signature(TreeNode node) {
        if (node.isLeaf()) return Integer.toString(node.taxonId);
        String left = signature(node.left);
        String right = signature(node.right);
        return left.compareTo(right) <= 0
            ? "(" + left + "," + right + ")"
            : "(" + right + "," + left + ")";
    }
}
