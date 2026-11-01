import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SimilarityMatrixBuilder;
import stelarx.gpu.GPUSimilarityMatrix;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.util.ArrayList;
import java.util.List;

/** Checks that bounded host batches preserve the quartet guide matrix. */
public final class CompactSimilarityHostBatchTest {
    private static void compare(SimilarityMatrix expected, SimilarityMatrix actual) {
        for (int a = 0; a < expected.n; a++) {
            for (int b = 0; b < expected.n; b++) {
                if (Math.abs(expected.getSim(a, b) - actual.getSim(a, b)) > 1e-12) {
                    throw new AssertionError("similarity mismatch at " + a + "," + b);
                }
            }
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass a gene-tree file");
        if (!GPUSimilarityMatrix.tryLoad()) {
            throw new IllegalStateException("CUDA similarity backend is unavailable");
        }
        Threading.start(Math.min(8, Runtime.getRuntime().availableProcessors()));
        try {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> input = TreeParser.parseGeneTrees(args[0], registry);
            List<Tree> repeated = new ArrayList<>(input.size() * 8);
            for (int i = 0; i < 8; i++) repeated.addAll(input);

            // One MiB forces several host batches on the supplied 40-taxon fixture.
            System.setProperty("stelarpro.similarity.compactHostBatchMiB", "1");
            SimilarityMatrix batched = SimilarityMatrixBuilder.buildGPU(repeated, registry.size());
            System.setProperty("stelarpro.similarity.compactHostBatchMiB", "1024");
            SimilarityMatrix oneBatch = SimilarityMatrixBuilder.buildGPU(repeated, registry.size());
            compare(oneBatch, batched);
            compare(SimilarityMatrixBuilder.buildCPU(repeated, registry.size()), batched);

            System.setProperty("stelarpro.similarity.forcePacked", "true");
            System.setProperty("stelarpro.similarity.compactHostBatchMiB", "1");
            SimilarityMatrix packed = SimilarityMatrixBuilder.buildGPU(repeated, registry.size());
            compare(oneBatch, packed);
            System.out.println("Compact GPU similarity host batching: PASS (dense and packed)");
        } finally {
            System.clearProperty("stelarpro.similarity.compactHostBatchMiB");
            System.clearProperty("stelarpro.similarity.forcePacked");
            Threading.shutdown();
        }
    }
}
