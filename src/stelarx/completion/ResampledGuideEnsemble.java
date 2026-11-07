package stelarx.completion;

import stelarx.Logging;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.SplittableRandom;

/**
 * An ensemble of NJ guide trees built from standard-error-perturbed resamples of the
 * quartet-similarity matrix, as many candidate sources instead of one.
 *
 * <h3>Why many trees rather than one more method</h3>
 * Measured on 200-taxon duplication data, a single global guide tree holds 199 rooted
 * clades; stacking five different distance methods reaches only 242.1, i.e. 1.22x.
 * Each method after the second adds ~11 clades and declining, so the single-tree axis
 * saturates. Resampling attacks that directly: k trees drawn from the matrix's own
 * uncertainty are not constrained to agree with each other, and they concentrate
 * disagreement exactly where the data is ambiguous.
 *
 * <h3>The perturbation is calibrated, not arbitrary</h3>
 * {@code sim} is a ratio of accumulated counts, {@code p = numSum / denSum}, so each
 * entry carries a binomial standard error that is already implied by the data:
 *
 * <pre>   se(a,b) = sqrt( p(1-p) / den(a,b) )   </pre>
 *
 * A member perturbs each pair by {@code z * se} with {@code z ~ N(0,1)}, one draw per
 * unordered pair, then applies the chosen similarity-to-distance transform. Pairs with
 * {@code den = 0} (never co-occur) and pairs at {@code p = 0} or {@code p = 1} have
 * {@code se = 0} and are left alone, so the perturbation cannot invent signal where
 * there is none. {@code scale} multiplies the standard error for sensitivity analysis;
 * 1.0 is the calibrated setting.
 *
 * <p>Note the known limitation: the Wald form understates uncertainty near {@code p = 0}
 * and {@code p = 1}, and at high gene counts {@code den} is large enough that {@code se}
 * may be too small to flip any NJ join. Both are why this class logs the {@code den} and
 * {@code se} distributions and the realized topological diversity of its own output -
 * an ensemble of k identical trees is a silent failure otherwise.
 *
 * <h3>What this is really probing</h3>
 * New clusters are no longer the interesting quantity: enrichment of the candidate set
 * was measured to raise the DP's objective without improving accuracy. k trees donate
 * k*(n-1) <em>splits</em>, and missing splits are what makes 40-55% of the clusters in
 * X unreachable and therefore discarded unscored. The reachability figure in
 * {@code WeightTable}'s prune line is the number to watch.
 */
public final class ResampledGuideEnsemble {

    private ResampledGuideEnsemble() {}

    /** One ensemble member: its name for per-source attribution, and its tree. */
    public record Member(String name, Tree tree) {}

    /**
     * Build k deduplicated NJ guide trees from perturbed resamples.
     *
     * @param firstTreeIndex tree index for the first member; members take consecutive
     *                       indices from here, so the caller must reserve k of them.
     * @return the members, possibly fewer than k after deduplication, or an empty list
     *         when the matrix is packed or the accumulators were discarded.
     */
    public static List<Member> build(SimilarityMatrix sim, int firstTreeIndex,
                                     NeighborJoiningClusterer.Variant variant,
                                     NeighborJoiningClusterer.Transform transform,
                                     int k, long seed, double scale, String label) {
        if (sim.isPacked() || sim.sim == null || sim.denSum == null) {
            Logging.warn("%s ensemble skipped: needs the dense similarity accumulators "
                + "(packed matrix at n=%d discards the denominator)", label, sim.n);
            return List.of();
        }
        int n = sim.n;
        if (n < 4 || k < 1) return List.of();

        int pairs = n * (n - 1) / 2;
        double[] p = new double[pairs], se = new double[pairs], den = new double[pairs];
        int[] pi = new int[pairs], pj = new int[pairs];
        int zeroDen = 0, zeroSe = 0, idx = 0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++, idx++) {
                double pv = sim.sim[i * n + j];
                double dv = sim.denSum[i * n + j];
                p[idx] = pv; den[idx] = dv; pi[idx] = i; pj[idx] = j;
                if (dv <= 0.0) { zeroDen++; se[idx] = 0.0; }
                else se[idx] = Math.sqrt(Math.max(0.0, pv * (1.0 - pv)) / dv);
                if (se[idx] <= 0.0) zeroSe++;
            }
        }
        Logging.info("%s ensemble: k=%d members, seed=%d, se scale=%.2f", label, k, seed, scale);
        Logging.info("%s ensemble   den per pair: %s", label, quantiles(den));
        Logging.info("%s ensemble    se per pair: %s  (%d pairs have se=0, %d have den=0)",
            label, quantiles(se), zeroSe, zeroDen);

        List<Member> members = new ArrayList<>(k);
        List<Set<BitSet>> cladeSets = new ArrayList<>(k);
        Set<Set<BitSet>> seen = new HashSet<>();
        double[] d = new double[n * n];
        int duplicates = 0, treeIndex = firstTreeIndex;

        for (int r = 0; r < k; r++) {
            SplittableRandom rng = new SplittableRandom(seed + r);
            for (int t = 0; t < pairs; t++) {
                double pv = p[t];
                if (se[t] > 0.0) {
                    pv = clamp01(pv + scale * se[t] * rng.nextGaussian());
                }
                double dv = NeighborJoiningClusterer.distance(pv, transform);
                d[pi[t] * n + pj[t]] = dv;
                d[pj[t] * n + pi[t]] = dv;
            }
            Tree member = NeighborJoiningClusterer.buildFromDistances(
                d, n, treeIndex, variant, transform);
            if (member == null) continue;
            Set<BitSet> clades = clades(member, n);
            if (!seen.add(clades)) { duplicates++; continue; }
            cladeSets.add(clades);
            members.add(new Member(String.format("%s~e%02d", label, r + 1), member));
            treeIndex++;
        }

        reportDiversity(label, k, duplicates, members, cladeSets, n);
        return members;
    }

    // ── Diagnostics ──────────────────────────────────────────────────────────

    /**
     * Report what the ensemble actually achieved. An ensemble whose members are all
     * the same tree is the expected failure mode when {@code se} is small relative to
     * the gaps between NJ's candidate joins, so this is not optional output.
     */
    private static void reportDiversity(String label, int k, int duplicates,
                                        List<Member> members,
                                        List<Set<BitSet>> cladeSets, int n) {
        if (members.isEmpty()) {
            Logging.warn("%s ensemble produced no usable member", label);
            return;
        }
        Set<BitSet> union = new LinkedHashSet<>();
        for (Set<BitSet> s : cladeSets) union.addAll(s);
        double jaccard = 0.0;
        int comparisons = 0;
        for (int a = 0; a < cladeSets.size(); a++) {
            for (int b = a + 1; b < cladeSets.size(); b++) {
                Set<BitSet> x = cladeSets.get(a), y = cladeSets.get(b);
                int inter = 0;
                for (BitSet c : x) if (y.contains(c)) inter++;
                int uni = x.size() + y.size() - inter;
                if (uni > 0) { jaccard += (double) inter / uni; comparisons++; }
            }
        }
        int perMember = cladeSets.get(0).size();
        Logging.info("%s ensemble: %d distinct topologies of %d drawn (%d duplicate(s) "
            + "dropped), %d clades per member, %d distinct clades in the union "
            + "(%.2fx one member), mean pairwise clade Jaccard %.3f",
            label, members.size(), k, duplicates, perMember, union.size(),
            perMember == 0 ? 0.0 : (double) union.size() / perMember,
            comparisons == 0 ? 1.0 : jaccard / comparisons);
        if (members.size() == 1 && k > 1) {
            Logging.warn("%s ensemble collapsed to a single topology: the standard error "
                + "is too small to move any NJ join. Raise the se scale or resample "
                + "genes instead of perturbing the matrix.", label);
        }
    }

    /** Rooted clades of size 2..n-1, as taxon-id bitsets. */
    private static Set<BitSet> clades(Tree tree, int n) {
        Set<BitSet> out = new HashSet<>();
        collect(tree.root, tree, n, out, true);
        return out;
    }

    private static void collect(TreeNode node, Tree tree, int n,
                                Set<BitSet> out, boolean isRoot) {
        if (node == null || node.isLeaf()) return;
        int size = node.rangeEnd - node.rangeStart;
        if (!isRoot && size >= 2 && size <= n - 1) {
            BitSet b = new BitSet(n);
            for (int pos = node.rangeStart; pos < node.rangeEnd; pos++) {
                b.set(tree.postorderArray[pos]);
            }
            out.add(b);
        }
        collect(node.left, tree, n, out, false);
        collect(node.right, tree, n, out, false);
    }

    private static String quantiles(double[] values) {
        double[] s = values.clone();
        Arrays.sort(s);
        return String.format("min=%.4g p25=%.4g median=%.4g p75=%.4g max=%.4g",
            s[0], at(s, 0.25), at(s, 0.50), at(s, 0.75), s[s.length - 1]);
    }

    private static double at(double[] sorted, double q) {
        return sorted[Math.min(sorted.length - 1, (int) (q * sorted.length))];
    }

    private static double clamp01(double v) {
        return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
    }
}
