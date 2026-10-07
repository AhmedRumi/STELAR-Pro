import stelarx.pro.AProGreedySearch;
import stelarx.pro.GeneTreeRooterTagger;
import stelarx.pro.S5OverallSearch;
import stelarx.Config;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/** In-process CPU/GPU smoke tests for S5's greedy proposal generator. */
public final class AProGreedySearchTest {
    public static void main(String[] args) throws Exception {
        if (!AProGreedySearch.tryLoad()) {
            throw new AssertionError("S5 greedy native library was not loaded");
        }
        if (args.length > 0 && args[0].equals("--large-smoke")) {
            largeTwoStepSmoke();
            return;
        }
        if (args.length > 0) {
            TaxonRegistry registry = new TaxonRegistry();
            Path tagged = Files.createTempFile("s5-greedy-tagged-", ".tre");
            try {
                GeneTreeRooterTagger.run(args[0], tagged.toString(), null, null);
                TreeParser.parseGeneTrees(tagged.toString(), registry);
                int rounds = args.length > 1 ? Integer.parseInt(args[1]) : 4;
                int threads = args.length > 2 ? Integer.parseInt(args[2]) : 2;
                long seed = args.length > 3 ? Long.parseLong(args[3]) : 233L;
                Threading.start(threads);
                try {
                    Config.getInstance().setComputeMode(Config.ComputeMode.CPU);
                    List<Tree> cpu = AProGreedySearch.generate(
                        args[0], tagged.toString(), registry, 0, threads, seed, rounds, false);
                    List<Tree> gpu = null;
                    if ((args.length < 5 || !args[4].equalsIgnoreCase("cpu-only"))
                            && AProGreedySearch.tryLoadGPU()) {
                        Config.getInstance().setComputeMode(Config.ComputeMode.GPU);
                        gpu = AProGreedySearch.generate(args[0], tagged.toString(),
                            registry, 0, threads, seed, rounds, true);
                    }
                    for (int i = 0; i < rounds; i++) {
                        String proposal = cpu.get(i).toNewick(registry);
                        if (gpu != null && !proposal.equals(gpu.get(i).toNewick(registry))) {
                            throw new AssertionError("CPU/GPU greedy-tree mismatch at round " + i);
                        }
                        System.out.println(proposal);
                    }
                    System.err.println("S5 " + rounds + "-round CPU"
                        + (gpu == null ? "" : "/GPU") + " comparison passed");
                } finally {
                    Threading.shutdown();
                }
            } finally {
                Files.deleteIfExists(tagged);
            }
            return;
        }
        Path input = Files.createTempFile("s5-greedy-test-", ".tre");
        try {
            Files.writeString(input, "(((A,A),B),(C,D));\n((A,C),(B,D));\n");
            TaxonRegistry registry = new TaxonRegistry();
            for (String name : List.of("A", "B", "C", "D")) registry.register(name);
            registry.lock();
            List<Tree> candidates = AProGreedySearch.generate(
                input.toString(), registry, 2, 2, 1234L, 2);
            if (candidates.size() != 2) throw new AssertionError("Expected two rounds");
            for (int i = 0; i < candidates.size(); i++) {
                Tree tree = candidates.get(i);
                if (tree.treeIndex != i + 2 || tree.leafCount != 4 || !tree.isComplete) {
                    throw new AssertionError("Invalid greedy tree " + i);
                }
                if (!tree.root.isSpeciation()) {
                    throw new AssertionError("Synthetic root is not a candidate node");
                }
            }
            if (AProGreedySearch.tryLoadGPU()) {
                List<Tree> gpu = AProGreedySearch.generate(
                    input.toString(), registry, 2, 2, 1234L, 2, true);
                for (int i = 0; i < candidates.size(); i++) {
                    String cpuTree = candidates.get(i).toNewick(registry);
                    String gpuTree = gpu.get(i).toNewick(registry);
                    if (!cpuTree.equals(gpuTree)) {
                        throw new AssertionError("CPU/GPU greedy-tree mismatch: "
                            + cpuTree + " versus " + gpuTree);
                    }
                }
                System.out.println("S5 CPU/GPU greedy smoke test passed");
            } else {
                System.out.println("S5 CPU greedy smoke test passed; CUDA unavailable");
            }
        } finally {
            Files.deleteIfExists(input);
        }
    }

    /** Exercise the >=100-taxon subset backbone rather than only direct placement. */
    private static void largeTwoStepSmoke() throws Exception {
        Path input = Files.createTempFile("s5-two-step-test-", ".tre");
        try {
            String gene = balancedTree(1, 101) + ";\n";
            Files.writeString(input, gene + gene);
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> genes = TreeParser.parseGeneTrees(input.toString(), registry);
            Threading.start(4);
            try {
                Config.getInstance().setComputeMode(Config.ComputeMode.CPU);
                List<Tree> cpu = AProGreedySearch.generate(
                    input.toString(), input.toString(), registry, 2, 4, 233L, 2, false);
                List<Tree> gpu = null;
                if (AProGreedySearch.tryLoadGPU()) {
                    Config.getInstance().setComputeMode(Config.ComputeMode.GPU);
                    gpu = AProGreedySearch.generate(
                        input.toString(), input.toString(), registry, 2, 4, 233L, 2, true);
                }
                if (cpu.size() != 2) throw new AssertionError("Expected two two-step rounds");
                for (int i = 0; i < cpu.size(); i++) {
                    Tree tree = cpu.get(i);
                    if (tree.treeIndex != i + 2 || tree.leafCount != 100 || !tree.isComplete) {
                        throw new AssertionError("Invalid large-taxon greedy tree " + i);
                    }
                    if (gpu != null && !tree.toNewick(registry)
                            .equals(gpu.get(i).toNewick(registry))) {
                        throw new AssertionError("CPU/GPU two-step mismatch at round " + i);
                    }
                }
                List<Tree> accumulated = new ArrayList<>(genes);
                accumulated.addAll(cpu);
                Config.getInstance().setComputeMode(Config.ComputeMode.CPU);
                String previous = S5OverallSearch.infer(
                    genes, accumulated, registry).newick();
                List<Tree> guidedCPU = AProGreedySearch.generate(
                    input.toString(), input.toString(), registry, accumulated.size(),
                    4, 234L, 2, false, previous);
                if (guidedCPU.size() != 2 || guidedCPU.get(0).leafCount != 100)
                    throw new AssertionError("Invalid guided S5 adaptive round");
                if (gpu != null) {
                    Config.getInstance().setComputeMode(Config.ComputeMode.GPU);
                    List<Tree> guidedGPU = AProGreedySearch.generate(
                        input.toString(), input.toString(), registry, accumulated.size(),
                        4, 234L, 2, true, previous);
                    for (int i = 0; i < guidedCPU.size(); i++) {
                        if (!guidedCPU.get(i).toNewick(registry)
                                .equals(guidedGPU.get(i).toNewick(registry)))
                            throw new AssertionError("CPU/GPU adaptive round mismatch " + i);
                    }
                }
                System.out.println("S5 large-taxon two-step CPU"
                    + (gpu == null ? "" : "/GPU") + " adaptive smoke test passed");
            } finally {
                Threading.shutdown();
            }
        } finally {
            Files.deleteIfExists(input);
        }
    }

    private static String balancedTree(int begin, int end) {
        if (end - begin == 1) return "S" + begin;
        int mid = (begin + end) / 2;
        return "(" + balancedTree(begin, mid) + "," + balancedTree(mid, end) + ")";
    }
}
