package stelarx.completion;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Self-tests for the guide-tree distance machinery. Run with:
 *
 * <pre>   java -cp build stelarx.completion.CompletionSelfTest   </pre>
 *
 * These check properties that are easy to get silently wrong and impossible to
 * notice downstream: a distance matrix that is merely plausible still produces a
 * guide tree, and a wrong guide tree only shows up as slightly worse RF weeks later.
 * The strongest check here is {@link #njRecoversTheGeneratingTree()} - NJ is provably
 * consistent on an additive matrix, so if the matrix is right NJ must return the exact
 * topology it was derived from, and if it is wrong that almost never holds.
 */
public final class CompletionSelfTest {

    private static int passed = 0, failed = 0;

    public static void main(String[] args) {
        internodeDistancesOnBalancedQuartet();
        internodeDistancesOnCaterpillar();
        rootCorrectionIsNotAConstantOffset();
        fourPointConditionHolds();
        njRecoversTheGeneratingTree();
        averagingAcrossTrees();
        closureFillsGapsWithShortestTwoLeg();
        pauplinEqualsBranchLengthSumOnAdditiveInput();
        standardErrorFollowsGeneTreeVariance();
        ensembleDiversifiesAndStaysSane();
        specMatchersAcceptOnlyWhatTheyShould();

        System.out.printf("%n%d passed, %d failed%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    // ── Internode distance ───────────────────────────────────────────────────

    /** ((A,B),(C,D)): siblings 2 apart, everything across the root 3. */
    private static void internodeDistancesOnBalancedQuartet() {
        var t = tree("((A,B),(C,D));", "A", "B", "C", "D");
        double[] d = InternodeDistanceMatrix.singleTreeForTest(t.tree.root, 4);
        check("quartet d(A,B)", 2.0, at(d, t, "A", "B", 4));
        check("quartet d(C,D)", 2.0, at(d, t, "C", "D", 4));
        check("quartet d(A,C)", 3.0, at(d, t, "A", "C", 4));
        check("quartet d(A,D)", 3.0, at(d, t, "A", "D", 4));
        check("quartet d(B,C)", 3.0, at(d, t, "B", "C", 4));
        check("quartet d(B,D)", 3.0, at(d, t, "B", "D", 4));
    }

    /** ((((A,B),C),D),E) unrooted is a 5-leaf caterpillar; check every pair by hand. */
    private static void internodeDistancesOnCaterpillar() {
        var t = tree("((((A,B),C),D),E);", "A", "B", "C", "D", "E");
        double[] d = InternodeDistanceMatrix.singleTreeForTest(t.tree.root, 5);
        check("cat d(A,B)", 2.0, at(d, t, "A", "B", 5));
        check("cat d(A,C)", 3.0, at(d, t, "A", "C", 5));
        check("cat d(A,D)", 4.0, at(d, t, "A", "D", 5));
        check("cat d(A,E)", 4.0, at(d, t, "A", "E", 5));   // root edge removed
        check("cat d(C,D)", 3.0, at(d, t, "C", "D", 5));
        check("cat d(C,E)", 3.0, at(d, t, "C", "E", 5));
        check("cat d(D,E)", 2.0, at(d, t, "D", "E", 5));
    }

    /**
     * The degree-2 root correction applies to a topology-dependent subset of pairs, so
     * forgetting it is not a constant offset that cancels. Two rootings of the SAME
     * unrooted tree must give the identical matrix; they only do if the correction is
     * applied where the LCA is the root.
     */
    private static void rootCorrectionIsNotAConstantOffset() {
        String[] taxa = {"A", "B", "C", "D", "E"};
        var t1 = tree("((((A,B),C),D),E);", taxa);
        var t2 = tree("(((A,B),C),(D,E));", taxa);
        double[] d1 = InternodeDistanceMatrix.singleTreeForTest(t1.tree.root, 5);
        double[] d2 = InternodeDistanceMatrix.singleTreeForTest(t2.tree.root, 5);
        // Both Newicks describe the same unrooted topology: A,B | C | D,E on a path.
        boolean same = true;
        for (String x : taxa) {
            for (String y : taxa) {
                if (x.equals(y)) continue;
                if (Math.abs(at(d1, t1, x, y, 5) - at(d2, t2, x, y, 5)) > 1e-9) same = false;
            }
        }
        record("two rootings of one unrooted topology agree", same,
            same ? "" : "matrices differ, so the root edge is being counted");
    }

    /**
     * A tree metric satisfies the four-point condition: for any four taxa, two of the
     * three pair-sums are equal and the third is no larger. This is the defining
     * property of additivity, which is what NJ's consistency assumes.
     */
    private static void fourPointConditionHolds() {
        String[] taxa = {"A", "B", "C", "D", "E", "F", "G", "H"};
        var t = tree("(((A,B),(C,D)),((E,F),(G,H)));", taxa);
        double[] d = InternodeDistanceMatrix.singleTreeForTest(t.tree.root, 8);
        int n = 8, violations = 0;
        for (int a = 0; a < n; a++)
          for (int b = a + 1; b < n; b++)
            for (int c = b + 1; c < n; c++)
              for (int e = c + 1; e < n; e++) {
                double s1 = d[a * n + b] + d[c * n + e];
                double s2 = d[a * n + c] + d[b * n + e];
                double s3 = d[a * n + e] + d[b * n + c];
                double[] s = {s1, s2, s3};
                Arrays.sort(s);
                if (Math.abs(s[1] - s[2]) > 1e-9) violations++;   // top two must tie
              }
        record("four-point condition on all 70 quadruples", violations == 0,
            violations + " violation(s): the matrix is not additive");
    }

    /**
     * The end-to-end check. NJ is consistent on an additive matrix, so running it on a
     * single tree's own internode distances must return that exact topology.
     */
    private static void njRecoversTheGeneratingTree() {
        String[] taxa = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};
        String nwk = "((((A,B),C),((D,E),F)),((G,H),(I,J)));";
        var t = tree(nwk, taxa);
        int n = taxa.length;
        var m = InternodeDistanceMatrix.build(List.of(t.tree), n);
        if (m == null) { record("NJ recovers the generating tree", false, "matrix was null"); return; }
        Tree nj = NeighborJoiningClusterer.buildFromDistances(
            m.distances(), n, 0, NeighborJoiningClusterer.Variant.NJ,
            NeighborJoiningClusterer.Transform.ONE_MINUS_SIM);
        var want = clades(t.tree, n);
        var got = clades(nj, n);
        // Compare unrooted: a clade and its complement are the same split.
        var wantSplits = splits(want, n);
        var gotSplits = splits(got, n);
        boolean ok = wantSplits.equals(gotSplits);
        record("NJ recovers the generating tree from its own internode matrix", ok,
            ok ? "" : wantSplits.size() + " wanted, " + gotSplits.size() + " got, "
                + diff(wantSplits, gotSplits) + " missing, "
                + diff(gotSplits, wantSplits) + " extra");
    }

    /** Two components, one pair in both: the entry must be the arithmetic mean. */
    private static void averagingAcrossTrees() {
        String[] taxa = {"A", "B", "C", "D"};
        var t1 = tree("((A,B),(C,D));", taxa);   // d(A,C) = 3
        var t2 = tree("((A,C),(B,D));", taxa);   // d(A,C) = 2
        var m = InternodeDistanceMatrix.build(List.of(t1.tree, t2.tree), 4);
        int ia = t1.reg.getId("A"), ic = t1.reg.getId("C");
        check("mean of 3 and 2", 2.5, m.distances()[ia * 4 + ic]);
        check("both components counted", 2.0, (double) m.stats().treesUsed());
        check("all 6 pairs observed", 6.0, (double) m.stats().pairsObserved());
        check("no pair missing", 0.0, (double) m.stats().pairsMissing());
    }

    /** A gap must be filled by the cheapest two-leg detour, not by any detour. */
    private static void closureFillsGapsWithShortestTwoLeg() {
        int n = 4;
        double[] d = new double[n * n];
        // Known: A-B 1, A-C 10, B-C 2, C-D 1, A-D 4. Unknown: B-D.
        set(d, n, 0, 1, 1); set(d, n, 0, 2, 10); set(d, n, 1, 2, 2);
        set(d, n, 2, 3, 1); set(d, n, 0, 3, 4);
        d[1 * n + 3] = Double.NaN; d[3 * n + 1] = Double.NaN;
        long filled = InternodeDistanceMatrix.closurePassForTest(d, n);
        check("closure filled one pair", 1.0, (double) filled);
        // via A: 1+4 = 5; via C: 2+1 = 3. Must take 3.
        check("closure takes the shortest detour", 3.0, d[1 * n + 3]);
        check("closure writes both triangles", 3.0, d[3 * n + 1]);
    }

    /** The dormant BME check, now actually run: Pauplin on an additive matrix. */
    private static void pauplinEqualsBranchLengthSumOnAdditiveInput() {
        String[] taxa = {"A", "B", "C", "D", "E", "F"};
        var t = tree("(((A,B),C),((D,E),F));", taxa);
        int n = taxa.length;
        var m = InternodeDistanceMatrix.build(List.of(t.tree), n);
        double got = BalancedMEOptimizer.lengthOfRootedForTest(m.distances(), n, t.tree);
        // Unrooted topology on 6 leaves: 2n-3 = 9 edges, each of unit length.
        check("Pauplin length equals the unit branch-length sum (9 edges)", 9.0, got);
    }

    /**
     * The standard error must come from variation across gene trees, not from a count
     * of derived quantities. A pair with the same distance in every tree has SE 0; a
     * pair that disagrees has SE = sd/sqrt(trees).
     */
    private static void standardErrorFollowsGeneTreeVariance() {
        String[] taxa = {"A", "B", "C", "D"};
        var t1 = tree("((A,B),(C,D));", taxa);   // d(A,B)=2, d(A,C)=3
        var t2 = tree("((A,B),(C,D));", taxa);   // identical
        var t3 = tree("((A,C),(B,D));", taxa);   // d(A,B)=3, d(A,C)=2
        int n = 4;

        var same = InternodeDistanceMatrix.build(List.of(t1.tree, t2.tree), n);
        int ia = t1.reg.getId("A"), ib = t1.reg.getId("B");
        check("SE is 0 when every tree agrees", 0.0, same.standardErrors()[ia * n + ib]);
        check("count is the number of trees, not quartets", 2.0,
            (double) same.counts()[ia * n + ib]);

        var mixed = InternodeDistanceMatrix.build(List.of(t1.tree, t3.tree), n);
        // distances 2 and 3: sample sd = 0.7071, SE = sd/sqrt(2) = 0.5
        check("SE equals sd/sqrt(trees) when they disagree", 0.5,
            mixed.standardErrors()[ia * n + ib]);
        check("mean of 2 and 3", 2.5, mixed.distances()[ia * n + ib]);
        record("median count is the gene-tree count, so no 100x fudge is needed",
            Math.abs(mixed.stats().medianCount() - 2.0) < 1e-9,
            "got " + mixed.stats().medianCount());
    }

    /**
     * An ensemble must produce more than one topology on data with real variation, must
     * never emit a negative distance, and must agree with the single-tree builder when
     * the perturbation is switched off.
     */
    private static void ensembleDiversifiesAndStaysSane() {
        String[] taxa = {"A", "B", "C", "D", "E", "F", "G", "H"};
        // Conflicting components, so the per-pair SE is nonzero.
        var a = tree("(((A,B),(C,D)),((E,F),(G,H)));", taxa);
        var b = tree("(((A,C),(B,D)),((E,G),(F,H)));", taxa);
        var c = tree("(((A,E),(B,F)),((C,G),(D,H)));", taxa);
        int n = taxa.length;
        var m = InternodeDistanceMatrix.build(List.of(a.tree, b.tree, c.tree), n);

        var off = m.ensemble("astrid-ensemble", 0, 4, 7L, 0.0, 10);
        record("scale 0 collapses to exactly one topology", off.size() == 1,
            "got " + off.size() + " member(s)");
        var njOnly = m.guide("astrid", 0, 10);
        record("scale 0 reproduces the unperturbed guide tree",
            ResampledGuideEnsembleAccess.sameClades(off.get(0).tree(), njOnly, n),
            "the zero-perturbation member differs from the plain guide");

        var on = m.ensemble("astrid-ensemble", 0, 12, 7L, 1.0, 10);
        record("scale 1 on conflicting input gives more than one topology", on.size() > 1,
            "got " + on.size() + " member(s); the SE may be degenerate");
        record("every member is a usable tree", on.stream().allMatch(
            x -> x.tree() != null && x.tree().root != null), "a member was null");
    }

    /** The CLI gatekeepers: a typo must be rejected, not silently ignored. */
    private static void specMatchersAcceptOnlyWhatTheyShould() {
        String[] good = {"astrid", "astrid-bionj", "astrid-fastme"};
        String[] bad = {"astrid-linear", "astrid-corrected", "astrid-nj", "astri", "nj"};
        boolean ok = true;
        for (String s : good) ok &= InternodeDistanceMatrix.isAstridSpec(s);
        for (String s : bad) ok &= !InternodeDistanceMatrix.isAstridSpec(s);
        record("isAstridSpec accepts the three specs and nothing else", ok, "mismatch");

        String[] goodE = {"astrid-ensemble:20", "astrid-ensemble:20:42",
                          "astrid-ensemble:20:42:1.5", "astrid-ensemble-bionj:5",
                          "astrid-ensemble-fastme:5:1:2"};
        String[] badE = {"astrid-ensemble", "astrid-ensemble:", "astrid-ensemble:abc",
                         "astrid:20", "astrid-ensemble:20:42:1.5:9"};
        boolean okE = true;
        for (String s : goodE) okE &= InternodeDistanceMatrix.isAstridEnsembleSpec(s);
        for (String s : badE) okE &= !InternodeDistanceMatrix.isAstridEnsembleSpec(s);
        record("isAstridEnsembleSpec accepts k[:seed[:scale]] and rejects malformed",
            okE, "mismatch");
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private record Built(Tree tree, TaxonRegistry reg) {}

    private static Built tree(String newick, String... taxa) {
        var reg = new TaxonRegistry();
        for (String s : taxa) reg.register(s);
        reg.lock();
        return new Built(TreeParser.parseSpeciesTreeNewick(newick, reg), reg);
    }

    private static double at(double[] d, Built b, String x, String y, int n) {
        return d[b.reg.getId(x) * n + b.reg.getId(y)];
    }

    private static void set(double[] d, int n, int a, int b, double v) {
        d[a * n + b] = v; d[b * n + a] = v;
    }

    private static List<java.util.BitSet> clades(Tree t, int n) {
        var out = new ArrayList<java.util.BitSet>();
        collect(t.root, n, out);
        return out;
    }

    private static java.util.BitSet collect(stelarx.tree.TreeNode node, int n,
                                            List<java.util.BitSet> out) {
        var bs = new java.util.BitSet(n);
        if (node.isLeaf()) { bs.set(node.taxonId); return bs; }
        bs.or(collect(node.left, n, out));
        bs.or(collect(node.right, n, out));
        if (bs.cardinality() > 1 && bs.cardinality() < n) out.add(bs);
        return bs;
    }

    /**
     * Canonicalize each clade to the side holding taxon 0, so the two sides of one edge
     * collapse to one split and rooting does not matter.
     *
     * <p>Splits with a singleton on either side are dropped. They are leaf edges, present
     * in every tree on the same taxa, and carry no topological claim. NJ returns a tree
     * rooted at its last join, which can be a leaf, and that rooting manufactures exactly
     * one such clade - counting it made this comparison fail against a correct matrix.
     */
    private static java.util.Set<java.util.BitSet> splits(List<java.util.BitSet> cs, int n) {
        var out = new java.util.HashSet<java.util.BitSet>();
        for (var c : cs) {
            if (c.cardinality() < 2 || n - c.cardinality() < 2) continue;
            var comp = new java.util.BitSet(n);
            comp.set(0, n);
            comp.andNot(c);
            out.add(c.get(0) ? c : comp);
        }
        return out;
    }

    private static int diff(java.util.Set<java.util.BitSet> a, java.util.Set<java.util.BitSet> b) {
        var copy = new java.util.HashSet<>(a);
        copy.removeAll(b);
        return copy.size();
    }

    private static void check(String name, double want, double got) {
        record(name, Math.abs(want - got) < 1e-9, "want " + want + ", got " + got);
    }

    private static void record(String name, boolean ok, String detail) {
        if (ok) { passed++; System.out.println("  ok    " + name); }
        else { failed++; System.out.println("  FAIL  " + name + "  -- " + detail); }
    }
}
