import stelarx.pro.AProGreedySearch;
import stelarx.pro.GeneTreeRooterTagger;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Small in-process smoke test for S5's CPU quartet-greedy tree generator. */
public final class AProGreedySearchTest {
    public static void main(String[] args) throws Exception {
        if (!AProGreedySearch.tryLoad()) {
            throw new AssertionError("S5 greedy native library was not loaded");
        }
        if (args.length > 0) {
            TaxonRegistry registry = new TaxonRegistry();
            Path tagged = Files.createTempFile("s5-greedy-tagged-", ".tre");
            try {
                GeneTreeRooterTagger.run(args[0], tagged.toString(), null, null);
                TreeParser.parseGeneTrees(tagged.toString(), registry);
            } finally {
                Files.deleteIfExists(tagged);
            }
            int rounds = args.length > 1 ? Integer.parseInt(args[1]) : 4;
            int threads = args.length > 2 ? Integer.parseInt(args[2]) : 2;
            long seed = args.length > 3 ? Long.parseLong(args[3]) : 233L;
            List<Tree> cpu = AProGreedySearch.generate(
                args[0], registry, 0, threads, seed, rounds);
            List<Tree> gpu = AProGreedySearch.tryLoadGPU()
                ? AProGreedySearch.generate(args[0], registry, 0, threads, seed, rounds, true)
                : null;
            for (int i = 0; i < rounds; i++) {
                String proposal = cpu.get(i).toNewick(registry);
                if (gpu != null && !proposal.equals(gpu.get(i).toNewick(registry))) {
                    throw new AssertionError("CPU/GPU greedy-tree mismatch at round " + i);
                }
                System.out.println(proposal);
            }
            System.err.println("S5 " + rounds + "-round CPU"
                + (gpu == null ? "" : "/GPU") + " comparison passed");
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
}
