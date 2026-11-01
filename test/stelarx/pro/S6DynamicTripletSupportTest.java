package stelarx.pro;

import stelarx.Config;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.util.Threading;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

/** Exhaustive score and proposal parity against the old dense test oracle. */
public final class S6DynamicTripletSupportTest {
    private static int[] group(int code, int n, int colour) {
        int[] result = new int[n];
        int count = 0;
        for (int taxon = 0; taxon < n; taxon++) {
            if (code % 4 == colour) result[count++] = taxon;
            code /= 4;
        }
        return Arrays.copyOf(result, count);
    }

    private static void check(S6TripletSupport dense, S6DynamicTripletSupport dynamic,
                              int n) {
        try (S6ScoreSource session = dynamic.fork()) {
            for (int x = 0; x < n; x++) {
                S6ScoreSource.Focus focus = session.focus(x);
                for (int a = 0; a < n; a++) for (int b = 0; b < n; b++)
                    for (int c = 0; c < n; c++) {
                        if (a != x && b != x && c != x) continue;
                        long expected = dense.support(a, b, c);
                        long found = focus.support(a, b, c);
                        if (found != expected)
                            throw new AssertionError("Focused triplet " + a + "," + b
                                + "|" + c + ": " + found + " != " + expected);
                    }
            }
            int limit = 1 << (2 * n);
            for (int code = 0; code < limit; code++) {
                int[] a = group(code, n, 0), b = group(code, n, 1);
                int[] c = group(code, n, 2);
                if (!Arrays.equals(session.groupScores(a, b, c),
                                   dense.groupScores(a, b, c)))
                    throw new AssertionError("Colour-score mismatch at " + code);
                if (session.splitWeight(a, b) != dense.splitWeight(a, b))
                    throw new AssertionError("Split-weight mismatch at " + code);
            }
        }
    }

    /** Restricting gene families must preserve selected-species triplets. */
    private static void checkRestricted(Path input, TaxonRegistry full,
                                        S6TripletSupport dense) throws Exception {
        TaxonRegistry subset = new TaxonRegistry();
        for (String name : List.of("A", "C", "D")) subset.register(name);
        subset.lock();
        List<Tree> genes = TreeParser.parseRestrictedGeneTrees(
            input.toString(), subset, false).trees();
        try (S6DynamicTripletSupport dynamic = S6DynamicTripletSupport.build(
                 genes, subset.size(), 2);
             S6ScoreSource scorer = dynamic.fork()) {
            for (int x = 0; x < 3; x++) {
                S6ScoreSource.Focus focused = scorer.focus(x);
                for (int a = 0; a < 3; a++)
                    for (int b = a + 1; b < 3; b++)
                        for (int c = 0; c < 3; c++) {
                            if (c == a || c == b ||
                                    (x != a && x != b && x != c)) continue;
                            int fa = full.getId(subset.getName(a));
                            int fb = full.getId(subset.getName(b));
                            int fc = full.getId(subset.getName(c));
                            if (focused.support(a, b, c) != dense.support(fa, fb, fc))
                                throw new AssertionError("Restricted triplet changed");
                        }
            }
        }
    }

    public static void main(String[] args) throws Exception {
        Path input = Files.createTempFile("s6-dynamic-parity-", ".tre");
        try {
            // A appears on both sides of a speciation root; D is duplication.
            Files.writeString(input, "((A,B),(A,C));\n"
                + "(((A,A)D,B),(C,D));\n"
                + "((A,B,C),D);\n");
            TaxonRegistry registry = new TaxonRegistry();
            List<Tree> genes = TreeParser.parseGeneTrees(input.toString(), registry, true);
            Config.getInstance().setComputeMode(Config.ComputeMode.CPU);
            S6TripletSupport dense = S6TripletSupport.build(genes, registry.size(), 2);
            checkRestricted(input, registry, dense);
            try (S6DynamicTripletSupport dynamic = S6DynamicTripletSupport.build(
                    genes, registry.size(), 2)) {
                check(dense, dynamic, registry.size());
                Threading.start(2);
                List<Tree> expected = S6GreedySearch.generate(dense, genes,
                    input.toString(), registry, genes.size(), 91L, 2, null);
                List<Tree> actual = S6GreedySearch.generate(dynamic, genes,
                    input.toString(), registry, genes.size(), 91L, 2, null);
                for (int i = 0; i < expected.size(); i++)
                    if (!expected.get(i).toNewick(registry)
                            .equals(actual.get(i).toNewick(registry)))
                        throw new AssertionError("Greedy proposal " + i + " changed");
            }
            if (args.length == 1 && args[0].equals("--gpu")) {
                Config.getInstance().setComputeMode(Config.ComputeMode.GPU);
                Config.getInstance().setGpuStrict(true);
                checkRestricted(input, registry, dense);
                try (S6DynamicTripletSupport dynamic = S6DynamicTripletSupport.build(
                        genes, registry.size(), 2)) {
                    check(dense, dynamic, registry.size());
                }
            }
            System.out.println("S6 dynamic triplet score and proposal parity: PASS");
        } finally {
            Threading.shutdown();
            Files.deleteIfExists(input);
        }
    }
}
