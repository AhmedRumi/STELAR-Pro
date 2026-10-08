package stelarx.pro;

import stelarx.Config;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.BitSet;
import java.util.List;

/** Exhaustive small-set parity with STELAR's binary split-weight equation. */
public final class S6TripletSupportTest {
    private static long manual(TreeNode node, Tree tree, BitSet a, BitSet b) {
        if (node.isLeaf()) return 0;
        long score = manual(node.left, tree, a, b)
            + manual(node.right, tree, a, b);
        if (!node.isSpeciation()) return score;
        BitSet m0 = members(node.left, tree), m1 = members(node.right, tree);
        int a0 = overlap(a, m0), a1 = overlap(a, m1);
        int b0 = overlap(b, m0), b1 = overlap(b, m1);
        return score + ((long) a0 * b1 * (a0 + b1 - 2L)
            + (long) a1 * b0 * (a1 + b0 - 2L)) / 2L;
    }

    private static BitSet members(TreeNode node, Tree tree) {
        BitSet result = new BitSet();
        for (int pos = node.rangeStart; pos < node.rangeEnd; pos++)
            result.set(tree.postorderArray[pos]);
        return result;
    }

    private static int overlap(BitSet a, BitSet b) {
        BitSet both = (BitSet) a.clone();
        both.and(b);
        return both.cardinality();
    }

    private static int[] values(BitSet set) {
        int[] result = new int[set.cardinality()];
        int i = 0;
        for (int bit = set.nextSetBit(0); bit >= 0; bit = set.nextSetBit(bit + 1))
            result[i++] = bit;
        return result;
    }

    public static void main(String[] args) throws Exception {
        Path input = Files.createTempFile("s6-triplet-support-test-", ".tre");
        try {
            Files.writeString(input,
                "(((A,A)D,B),(C,D));\n"
              + "((A,B),(C,D));\n"
              + "(((A,C),B),D);\n");
            Threading.start(2);
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> genes = TreeParser.parseGeneTrees(input.toString(), registry);
            S6TripletSupport support = S6TripletSupport.build(genes,
                registry.size(), 2);
            S6GreedySearch.verifyInsertionDeltas(support);
            List<Tree> proposals = S6GreedySearch.generate(support, genes,
                input.toString(), registry, genes.size(), 233L, 2, null);
            for (Tree proposal : proposals) {
                if (proposal.leafCount != registry.size())
                    throw new AssertionError("S6 proposal omitted a species");
                BitSet present = members(proposal.root, proposal);
                if (present.cardinality() != registry.size())
                    throw new AssertionError("S6 proposal repeated a species");
                if (!proposal.root.isSpeciation())
                    throw new AssertionError("S6 proposal is not speciation rooted");
            }
            if (args.length == 1 && args[0].equals("--gpu")) {
                Config.getInstance().setComputeMode(Config.ComputeMode.GPU);
                Config.getInstance().setGpuStrict(true);
                S6TripletSupport device = S6TripletSupport.build(genes,
                    registry.size(), 2);
                for (int c = 0; c < registry.size(); c++)
                    for (int a = 0; a < registry.size(); a++)
                        for (int b = a + 1; b < registry.size(); b++)
                            if (device.support(a, b, c) != support.support(a, b, c))
                                throw new AssertionError("S6 CPU/GPU tensor mismatch");
                Config.getInstance().setComputeMode(Config.ComputeMode.CPU);
            }
            int tested = 0;
            for (int code = 0; code < 81; code++) {
                int state = code;
                BitSet a = new BitSet(), b = new BitSet();
                for (int taxon = 0; taxon < 4; taxon++) {
                    int color = state % 3;
                    state /= 3;
                    if (color == 1) a.set(taxon);
                    if (color == 2) b.set(taxon);
                }
                if (a.isEmpty() || b.isEmpty()) continue;
                long expected = 0;
                for (Tree gene : genes) expected += manual(gene.root, gene, a, b);
                long actual = support.splitWeight(values(a), values(b));
                if (actual != expected) {
                    throw new AssertionError("S6 support mismatch: " + a + " | " + b
                        + " expected " + expected + " observed " + actual);
                }
                tested++;
            }
            System.out.println("S6 triplet support parity passed (" + tested
                + " disjoint splits; duplicate and speciation roots)");
        } finally {
            Threading.shutdown();
            Files.deleteIfExists(input);
        }
    }
}
