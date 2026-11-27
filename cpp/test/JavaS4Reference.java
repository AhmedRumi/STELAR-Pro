import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SimilarityMatrixBuilder;
import stelarx.completion.SortedRows;
import stelarx.completion.SortedRowsBuilder;
import stelarx.completion.TreeCompleter;
import stelarx.completion.UPGMAClusterer;
import stelarx.pro.DiscoDecomposer;
import stelarx.pro.ScopedTreeCompleter;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.io.PrintWriter;
import java.util.List;

/** Machine-readable CPU reference; uses the unchanged production Java classes. */
public final class JavaS4Reference {
    public static void main(String[] args) throws Exception {
        Threading.start(Integer.parseInt(args[2]));
        try (PrintWriter out = new PrintWriter(args[1])) {
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> genes = TreeParser.parseGeneTrees(args[0], registry);
            out.println("N " + registry.size());
            for (int taxon = 0; taxon < registry.size(); ++taxon)
                out.println("T " + taxon + " " + registry.getName(taxon));
            for (DiscoDecomposer.CompletionScope scope : DiscoDecomposer.CompletionScope.values()) {
                DiscoDecomposer.Result d = DiscoDecomposer.decomposeAll(
                    genes, registry.size(), 0, 4, scope);
                out.println("S " + scope + " " + d.duplicationCuts() + " "
                    + d.discardedSmallTrees() + " " + d.trees().size());
                if (d.trees().isEmpty()) continue;
                for (Tree tree : d.trees()) {
                    out.print("R " + tree.toNewick(registry) + " M");
                    if (tree instanceof DiscoDecomposer.ScopedTree scoped) {
                        long[] words = scoped.missingWords();
                        for (int taxon = 0; taxon < registry.size(); ++taxon)
                            if (taxon / 64 < words.length && (words[taxon / 64] & (1L << taxon)) != 0)
                                out.print(" " + taxon);
                    }
                    out.println();
                }
                SimilarityMatrix matrix = SimilarityMatrixBuilder.buildTripletCPU(d.trees(), registry.size());
                SortedRows rows = SortedRowsBuilder.build(matrix);
                for (int a = 0; a < registry.size(); ++a) {
                    for (int b = a + 1; b < registry.size(); ++b)
                        out.println("P " + a + " " + b + " " + Double.toHexString(matrix.getSim(a,b)));
                    out.print("A " + a);
                    for (int rank = 0; rank < registry.size(); ++rank) out.print(" " + rows.get(a,rank));
                    out.println();
                }
                out.println("U " + UPGMAClusterer.build(matrix,0).toNewick(registry));
                List<Tree> completed = scope == DiscoDecomposer.CompletionScope.GLOBAL
                    ? TreeCompleter.completeAll(d.trees(), rows, registry.size())
                    : ScopedTreeCompleter.prepare(rows,d.trees(),registry.size()).completeBatch(d.trees());
                for (Tree tree : completed) out.println("C " + tree.toNewick(registry));
            }
        } finally { Threading.shutdown(); }
    }
}
