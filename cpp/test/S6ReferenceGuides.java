import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import stelarx.Config;
import stelarx.pro.S6S2BBGuideData;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;
import stelarx.cluster.*;
import stelarx.dp.*;
import stelarx.hash.*;
import stelarx.pro.UniqueTaxonSubtreeHashes;

/** Diagnostic helper: emit the first subset's proposals using retained Java. */
public final class S6ReferenceGuides {
    public static void main(String[] args) throws Exception {
        Config cfg = Config.getInstance();
        cfg.setComputeMode(Config.ComputeMode.CPU);
        cfg.setThreadCount(1);
        cfg.setSearchSpace(Config.SearchSpace.S6);
        Threading.start(1);
        try {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> genes = TreeParser.parseGeneTrees(args[0], registry);
            int n = registry.size();
            List<Integer> order = new ArrayList<>();
            for (int t = 0; t < n; ++t) order.add(t);
            Collections.shuffle(order, new Random(cfg.getS6Seed()));
            int m = Math.max(4, Math.min(n, (int)(Math.sqrt(n) * Math.log(n) / Math.log(2.0) / 4.0)));
            TaxonRegistry subset = new TaxonRegistry();
            for (int t : order.subList(0,m)) subset.register(registry.getName(t));
            subset.lock();
            List<Tree> guides = S6S2BBGuideData.prepare(genes, registry).restrict(subset);
            StringBuilder out = new StringBuilder();
            for (Tree tree : guides) out.append(tree.toNewick(subset)).append('\n');
            Files.writeString(Path.of(args[1]), out);
            TaxonHasher hasher = new TaxonHasher(m, cfg.getNumHashSeeds(), cfg.getBaseSeed());
            PrefixHashArrays pref = new PrefixHashArrays(guides, hasher);
            UniqueTaxonSubtreeHashes hashes = new UniqueTaxonSubtreeHashes(guides, hasher);
            ClusterTable clusters = new ClusterTable(guides, pref, m, hashes);
            DPTable dp = new DPTable(guides, pref, clusters, hashes, Integer.MAX_VALUE);
            Set<ClusterHash> reachable = dp.reachableClusters();
            List<String> lines = new ArrayList<>();
            for (var entry : dp.entries()) if (reachable.contains(entry.getKey()))
                for (BipartitionSplit split : entry.getValue()) {
                    String a = names(split.lo, clusters, guides, subset);
                    String b = names(split.hi, clusters, guides, subset);
                    lines.add(a.compareTo(b) < 0 ? a + " | " + b : b + " | " + a);
                }
            Collections.sort(lines);
            Files.writeString(Path.of(args[1] + ".splits"), String.join("\n", lines) + "\n");
        } finally { Threading.shutdown(); }
    }
    static String names(ClusterHash hash, ClusterTable clusters, List<Tree> trees, TaxonRegistry registry) {
        Cluster c = clusters.get(hash).exemplar;
        Tree tree = trees.get(c.treeIndex);
        List<String> names = new ArrayList<>();
        for (int i = c.left; i < c.right; ++i) names.add(registry.getName(tree.postorderArray[i]));
        Collections.sort(names);
        return String.join(",", names);
    }
}
