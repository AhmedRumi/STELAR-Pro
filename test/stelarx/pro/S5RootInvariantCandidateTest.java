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
import stelarx.util.Threading;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** S5's greedy candidates must not depend on A-Pro's arbitrary printed root. */
public final class S5RootInvariantCandidateTest {
    private static void tagProposal(TreeNode node) {
        if (node.isLeaf()) return;
        node.isSpeciationNode = true;
        tagProposal(node.left);
        tagProposal(node.right);
    }

    private static Map<ClusterHash, Set<BipartitionSplit>> candidates(
            String newick, TaxonRegistry registry) {
        Tree tree = TreeParser.parseSpeciesTreeNewick(newick, registry);
        tagProposal(tree.root);
        List<Tree> trees = List.of(tree);
        TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 233L);
        PrefixHashArrays prefix = new PrefixHashArrays(trees, hasher);
        UniqueTaxonSubtreeHashes unique = new UniqueTaxonSubtreeHashes(trees, hasher);
        ClusterTable clusters = new ClusterTable(
            trees, prefix, registry.size(), unique, 0);
        DPTable dp = new DPTable(trees, prefix, clusters, unique, 0);
        unique.release();
        dp.addCrossTreeTransitions(clusters, false);
        Map<ClusterHash, Set<BipartitionSplit>> result = new HashMap<>();
        for (var entry : dp.entries()) {
            result.put(entry.getKey(), new HashSet<>(entry.getValue()));
        }
        return result;
    }

    public static void main(String[] args) {
        Threading.start(2);
        try {
            TaxonRegistry registry = new TaxonRegistry();
            for (String taxon : List.of("A", "B", "C", "D")) registry.register(taxon);
            registry.lock();
            var first = candidates("((A,B),(C,D));", registry);
            var rerooted = candidates("(A,(B,(C,D)));", registry);
            if (!first.equals(rerooted)) {
                throw new AssertionError("S5 candidate DP depends on greedy Newick root");
            }
            TaxonRegistry fiveTaxa = new TaxonRegistry();
            for (String taxon : List.of("A", "B", "C", "D", "E")) {
                fiveTaxa.register(taxon);
            }
            fiveTaxa.lock();
            var fiveFirst = candidates("(((A,B),C),(D,E));", fiveTaxa);
            var fiveRerooted = candidates("(A,(B,(C,(D,E))));", fiveTaxa);
            if (!fiveFirst.equals(fiveRerooted)) {
                throw new AssertionError("S5 internal-edge rotations depend on Newick root");
            }
            System.out.println("S5 greedy-root invariance passed (" + first.size()
                + " four-taxon, " + fiveFirst.size() + " five-taxon DP parent clusters)");
        } finally {
            Threading.shutdown();
        }
    }
}
