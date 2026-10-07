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

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * In-process ASTRAL-Pro placement/NNI search for S5 candidate trees.
 * Its quartet objective proposes topologies only; original tagged gene trees
 * remain the sole observations scored by STELAR-Pro's triplet DP.
 */
public final class AProGreedySearch {
    private static boolean cpuLoadAttempted;
    private static boolean cpuLoaded;
    private static boolean gpuLoadAttempted;
    private static boolean gpuLoaded;

    private AProGreedySearch() {}

    public static synchronized boolean tryLoad() {
        if (!cpuLoadAttempted) {
            cpuLoadAttempted = true;
            try {
                System.loadLibrary("stelar_pro_s5_greedy");
                cpuLoaded = true;
            } catch (UnsatisfiedLinkError e) {
                cpuLoaded = false;
            }
        }
        return cpuLoaded;
    }

    public static synchronized boolean tryLoadGPU() {
        if (!gpuLoadAttempted) {
            gpuLoadAttempted = true;
            try {
                System.loadLibrary("stelar_pro_s5_greedy_gpu");
                gpuLoaded = true;
            } catch (UnsatisfiedLinkError e) {
                gpuLoaded = false;
            }
        }
        return gpuLoaded;
    }

    /** Build complete, rooted, single-copy candidate trees in native round order. */
    public static List<Tree> generate(String inputFile, TaxonRegistry registry,
                                      int firstTreeIndex, int threads, long seed,
                                      int rounds) {
        return generate(inputFile, registry, firstTreeIndex, threads, seed,
            rounds, false);
    }

    /** GPU mode uses the same placement/NNI search with exact device counters. */
    public static List<Tree> generate(String inputFile, TaxonRegistry registry,
                                      int firstTreeIndex, int threads, long seed,
                                      int rounds, boolean gpu) {
        return generate(inputFile, null, registry, firstTreeIndex,
            threads, seed, rounds, gpu);
    }

    /** The tagged input is needed only for the n>=100 triplet-DP backbone. */
    public static List<Tree> generate(String inputFile, String taggedInputFile,
                                      TaxonRegistry registry, int firstTreeIndex,
                                      int threads, long seed, int rounds,
                                      boolean gpu) {
        return generate(inputFile, taggedInputFile, registry, firstTreeIndex,
            threads, seed, rounds, gpu, null);
    }

    /** A prior STELAR-DP tree seeds each greedy round of an adaptive pass. */
    public static List<Tree> generate(String inputFile, String taggedInputFile,
                                      TaxonRegistry registry, int firstTreeIndex,
                                      int threads, long seed, int rounds,
                                      boolean gpu, String previousTree) {
        if (gpu ? !tryLoadGPU() : !tryLoad()) {
            throw new IllegalStateException(
                "S5 native greedy library is unavailable: " + (gpu
                    ? "libstelar_pro_s5_greedy_gpu.so"
                    : "libstelar_pro_s5_greedy.so"));
        }
        String[] newicks = gpu
            ? generateGPU(inputFile, taggedInputFile, previousTree,
                threads, seed, rounds)
            : generateCPU(inputFile, taggedInputFile, previousTree,
                threads, seed, rounds);
        List<Tree> trees = new ArrayList<>(newicks.length);
        for (int i = 0; i < newicks.length; i++) {
            Tree parsed = TreeParser.parseSpeciesTreeNewick(newicks[i], registry);
            markProposalNodes(parsed.root);
            trees.add(new Tree(firstTreeIndex + i, parsed.root,
                parsed.postorderArray, parsed.positionMap, parsed.leafCount,
                registry.size(), parsed.hasPolytomy));
        }
        return trees;
    }

    /** Called by native two-step search after building subset guide trees. */
    private static String inferTripletBackbone(String taggedInputFile,
                                               String[] guideNewicks,
                                               String[] selectedTaxa,
                                               boolean gpu) throws IOException {
        if (taggedInputFile == null || guideNewicks.length == 0
                || selectedTaxa.length < 4) {
            throw new IllegalArgumentException("Invalid S5 two-step backbone input");
        }
        TaxonRegistry subset = new TaxonRegistry();
        for (String name : selectedTaxa) subset.register(name);
        subset.lock();
        if (subset.size() != selectedTaxa.length) {
            throw new IllegalArgumentException("S5 backbone subset contains duplicate taxa");
        }

        List<Tree> guides = new ArrayList<>(guideNewicks.length);
        for (int i = 0; i < guideNewicks.length; i++) {
            Tree parsed = TreeParser.parseSpeciesTreeNewick(guideNewicks[i], subset);
            markProposalNodes(parsed.root);
            guides.add(new Tree(i, parsed.root, parsed.postorderArray,
                parsed.positionMap, parsed.leafCount, subset.size(), false));
        }

        // Restrict only the ORIGINAL tagged multicopy trees to the backbone taxa.
        // Greedy guides contribute candidates; quartet scores do not enter this DP.
        List<Tree> genes = TreeParser.parseRestrictedGeneTrees(
            taggedInputFile, subset, false).trees();
        TaxonHasher hasher = new TaxonHasher(subset.size(),
            Config.getInstance().getNumHashSeeds(), Config.getInstance().getBaseSeed());
        PrefixHashArrays guidePrefix = new PrefixHashArrays(guides, hasher);
        PrefixHashArrays genePrefix = new PrefixHashArrays(genes, hasher);
        UniqueTaxonSubtreeHashes guideHashes =
            new UniqueTaxonSubtreeHashes(guides, hasher);
        UniqueTaxonSubtreeHashes geneHashes =
            new UniqueTaxonSubtreeHashes(genes, hasher);
        ClusterTable clusters = new ClusterTable(guides, guidePrefix,
            subset.size(), guideHashes, 0);
        PartitionTable partitions = new PartitionTable(genes, genePrefix, geneHashes);
        DPTable dp = new DPTable(guides, guidePrefix, clusters, guideHashes, 0);
        guideHashes.release();
        geneHashes.release();
        dp.addCrossTreeTransitions(clusters, gpu && GPUDPBuilder.tryLoad());
        WeightTable weights = new WeightTable(dp, partitions, clusters, guides, genes);
        Inference inference = new Inference();
        String rooted = inference.run(dp, weights, clusters, guides, subset, hasher);
        Logging.info("S5 two-step backbone: %d taxa, %d quartet-greedy guides, "
            + "STELAR triplet score %s", subset.size(), guides.size(),
            inference.getLastTripletScore());

        // A-Pro places onto an unrooted backbone anchored at taxon 0. Retain
        // the triplet-DP-selected topology while changing only its Newick root.
        Tree optimal = TreeParser.parseSpeciesTreeNewick(rooted, subset);
        return anchorAtLeaf(optimal, subset.getId(selectedTaxa[0]), subset);
    }

    private static String anchorAtLeaf(Tree tree, int taxonId,
                                       TaxonRegistry registry) {
        TreeNode leaf = findLeaf(tree.root, taxonId);
        if (leaf == null) throw new IllegalStateException("S5 backbone anchor is missing");
        List<TreeNode> adjacent = neighbors(leaf, tree.root);
        if (adjacent.size() != 1) {
            throw new IllegalStateException("S5 backbone anchor is not a leaf edge");
        }
        return "(" + registry.getName(taxonId) + ","
            + renderAwayFrom(adjacent.get(0), leaf, tree.root, registry) + ");";
    }

    /** Induce the previous STELAR tree on an A-Pro subset, then anchor it for placement. */
    private static String inducePreviousBackbone(String previousTree,
                                                 String[] allTaxa,
                                                 String[] selectedTaxa) {
        if (previousTree == null || selectedTaxa.length < 4)
            throw new IllegalArgumentException("Invalid S5 prior-backbone input");
        TaxonRegistry all = new TaxonRegistry();
        for (String name : allTaxa) all.register(name);
        all.lock();
        if (all.size() != allTaxa.length)
            throw new IllegalArgumentException("Duplicate taxon in S5 prior tree");
        BitSet keep = new BitSet(all.size());
        TaxonRegistry subset = new TaxonRegistry();
        for (String name : selectedTaxa) {
            keep.set(all.getId(name));
            subset.register(name);
        }
        subset.lock();
        if (subset.size() != selectedTaxa.length)
            throw new IllegalArgumentException("Duplicate taxon in S5 subset");
        Tree full = TreeParser.parseSpeciesTreeNewick(previousTree, all);
        String induced = renderSelected(full.root, keep, all);
        if (induced == null) throw new IllegalStateException("S5 subset is absent from prior tree");
        Tree restricted = TreeParser.parseSpeciesTreeNewick(induced + ";", subset);
        return anchorAtLeaf(restricted, subset.getId(selectedTaxa[0]), subset);
    }

    /** Contract unary paths introduced by removing taxa outside the subset. */
    private static String renderSelected(TreeNode node, BitSet keep,
                                         TaxonRegistry all) {
        if (node.isLeaf())
            return keep.get(node.taxonId) ? all.getName(node.taxonId) : null;
        String left = renderSelected(node.left, keep, all);
        String right = renderSelected(node.right, keep, all);
        if (left == null) return right;
        if (right == null) return left;
        return "(" + left + "," + right + ")";
    }

    private static TreeNode findLeaf(TreeNode node, int taxonId) {
        if (node.isLeaf()) return node.taxonId == taxonId ? node : null;
        TreeNode left = findLeaf(node.left, taxonId);
        return left != null ? left : findLeaf(node.right, taxonId);
    }

    /** Treat the original degree-two Newick root as a suppressed unrooted edge. */
    private static List<TreeNode> neighbors(TreeNode node, TreeNode root) {
        List<TreeNode> result = new ArrayList<>(3);
        if (node.parent == root) {
            result.add(root.left == node ? root.right : root.left);
        } else if (node.parent != null) {
            result.add(node.parent);
        }
        if (!node.isLeaf()) {
            result.add(node.left);
            result.add(node.right);
        }
        return result;
    }

    private static String renderAwayFrom(TreeNode node, TreeNode from,
                                         TreeNode root, TaxonRegistry registry) {
        if (node.isLeaf()) return registry.getName(node.taxonId);
        List<TreeNode> children = new ArrayList<>(2);
        for (TreeNode next : neighbors(node, root)) {
            if (next != from) children.add(next);
        }
        if (children.size() != 2) {
            throw new IllegalStateException("S5 backbone is not an unrooted binary tree");
        }
        return "(" + renderAwayFrom(children.get(0), node, root, registry)
            + "," + renderAwayFrom(children.get(1), node, root, registry) + ")";
    }

    /** Every greedy-tree node is a proposal, not a gene-duplication observation. */
    private static void markProposalNodes(TreeNode node) {
        if (node.isLeaf()) return;
        node.isDuplicationNode = false;
        node.isSpeciationNode = true;
        markProposalNodes(node.left);
        markProposalNodes(node.right);
    }

    private static native String[] generateCPU(String inputFile, String taggedInputFile,
                                                String previousTree,
                                                int threads,
                                                long seed, int rounds);
    private static native String[] generateGPU(String inputFile, String taggedInputFile,
                                                String previousTree,
                                                int threads,
                                                long seed, int rounds);
}
