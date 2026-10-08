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

/** Full STELAR triplet DP over S1 and accumulated S6 rooted proposals. */
public final class S6OverallSearch {
    private S6OverallSearch() {}

    public record Result(String newick, String tripletScore) {}

    /** Compare rooted topologies without depending on S5's search path. */
    public static boolean sameRootedTopology(String first, String second,
                                             TaxonRegistry registry) {
        Tree a = TreeParser.parseSpeciesTreeNewick(first, registry);
        Tree b = TreeParser.parseSpeciesTreeNewick(second, registry);
        return signature(a.root).equals(signature(b.root));
    }

    private static String signature(TreeNode node) {
        if (node.isLeaf()) return Integer.toString(node.taxonId);
        String left = signature(node.left), right = signature(node.right);
        return left.compareTo(right) <= 0
            ? "(" + left + "," + right + ")"
            : "(" + right + "," + left + ")";
    }

    public static Result infer(List<Tree> originalGenes, List<Tree> candidates,
                               TaxonRegistry registry) {
        Config cfg = Config.getInstance();
        TaxonHasher hasher = new TaxonHasher(registry.size(),
            cfg.getNumHashSeeds(), cfg.getBaseSeed());
        PrefixHashArrays candidatePrefix = new PrefixHashArrays(candidates, hasher);
        PrefixHashArrays genePrefix = new PrefixHashArrays(originalGenes, hasher);
        UniqueTaxonSubtreeHashes hashes =
            new UniqueTaxonSubtreeHashes(candidates, hasher);
        ClusterTable clusters = new ClusterTable(candidates, candidatePrefix,
            registry.size(), hashes);
        PartitionTable partitions = new PartitionTable(originalGenes, genePrefix, hashes);
        DPTable dp = new DPTable(candidates, candidatePrefix, clusters, hashes,
            Integer.MAX_VALUE);
        hashes.release();
        if (cfg.getSearchMode() == Config.SearchMode.FULL)
            dp.addCrossTreeTransitions(clusters,
                cfg.getComputeMode() == Config.ComputeMode.GPU && GPUDPBuilder.tryLoad());
        WeightTable weights = new WeightTable(dp, partitions, clusters,
            candidates, originalGenes);
        Inference inference = new Inference();
        String tree = inference.run(dp, weights, clusters, candidates, registry, hasher);
        Logging.info("S6 accumulated DP: %d S1 genes + %d rooted triplet-greedy "
            + "trees; score %s", originalGenes.size(),
            candidates.size() - originalGenes.size(), inference.getLastTripletScore());
        return new Result(tree, inference.getLastTripletScore());
    }
}
