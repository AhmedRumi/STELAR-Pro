package stelarx.pro;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;

import java.util.ArrayList;
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
        if (gpu ? !tryLoadGPU() : !tryLoad()) {
            throw new IllegalStateException(
                "S5 native greedy library is unavailable: " + (gpu
                    ? "libstelar_pro_s5_greedy_gpu.so"
                    : "libstelar_pro_s5_greedy.so"));
        }
        String[] newicks = gpu
            ? generateGPU(inputFile, threads, seed, rounds)
            : generateCPU(inputFile, threads, seed, rounds);
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

    /** Every greedy-tree node is a proposal, not a gene-duplication observation. */
    private static void markProposalNodes(TreeNode node) {
        if (node.isLeaf()) return;
        node.isDuplicationNode = false;
        node.isSpeciationNode = true;
        markProposalNodes(node.left);
        markProposalNodes(node.right);
    }

    private static native String[] generateCPU(String inputFile, int threads,
                                                long seed, int rounds);
    private static native String[] generateGPU(String inputFile, int threads,
                                                long seed, int rounds);
}
