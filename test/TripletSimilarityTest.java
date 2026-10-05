import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SimilarityMatrixBuilder;
import stelarx.gpu.GPUSimilarityMatrix;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Hand-calculated oracle and packed-storage parity for S4 guide similarity. */
public final class TripletSimilarityTest {
    public static void main(String[] args) throws Exception {
        Path input = Files.createTempFile("stelar-pro-triplet-sim-", ".tre");
        Files.writeString(input, "((A,B),(C,D));\n((A,C),B);\n");
        Threading.start(Math.min(4, Runtime.getRuntime().availableProcessors()));
        try {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> trees = TreeParser.parseGeneTrees(input.toString(), registry);
            SimilarityMatrix dense = SimilarityMatrixBuilder.buildTripletCPU(trees, 4);
            check(dense, registry.getId("A"), registry.getId("B"), 2.0 / 3.0);
            check(dense, registry.getId("A"), registry.getId("C"), 1.0 / 3.0);
            check(dense, registry.getId("C"), registry.getId("D"), 1.0);
            check(dense, registry.getId("A"), registry.getId("D"), 0.0);

            System.setProperty("stelarpro.similarity.forcePacked", "true");
            SimilarityMatrix packed;
            try {
                packed = SimilarityMatrixBuilder.buildTripletCPU(trees, 4);
            } finally {
                System.clearProperty("stelarpro.similarity.forcePacked");
            }
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    check(packed, a, b, dense.getSim(a, b));
                }
            }
            if (args.length > 0 && args[0].equals("--gpu")) {
                if (!GPUSimilarityMatrix.tryLoad()) {
                    throw new AssertionError("CUDA similarity backend did not load");
                }
                SimilarityMatrix gpu = SimilarityMatrixBuilder.buildTripletGPU(trees, 4);
                compare(dense, gpu);
                System.setProperty("stelarpro.similarity.forceWide", "true");
                System.setProperty("stelarpro.similarity.forcePacked", "true");
                try {
                    compare(dense, SimilarityMatrixBuilder.buildTripletGPU(trees, 4));
                } finally {
                    System.clearProperty("stelarpro.similarity.forceWide");
                    System.clearProperty("stelarpro.similarity.forcePacked");
                }
            }
        } finally {
            Threading.shutdown();
            Files.deleteIfExists(input);
        }
        System.out.println("Rooted-triplet similarity oracle/parity: PASS");
    }

    private static void check(SimilarityMatrix matrix, int a, int b, double expected) {
        double actual = matrix.getSim(a, b);
        if (Math.abs(actual - expected) > 1e-12) {
            throw new AssertionError("pair " + a + "," + b
                + ": expected=" + expected + " actual=" + actual);
        }
    }

    private static void compare(SimilarityMatrix expected, SimilarityMatrix actual) {
        for (int a = 0; a < expected.n; a++) {
            for (int b = 0; b < expected.n; b++) {
                check(actual, a, b, expected.getSim(a, b));
            }
        }
    }
}
