package stelarx.completion;

import stelarx.Logging;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.List;

/**
 * ASTRID-style average internode distance matrix, and guide trees built from it.
 *
 * <h3>Why this and not another transform of the quartet matrix</h3>
 * Every guide method before this one reads the same {@link SimilarityMatrix}: the
 * transforms change <em>how</em> it is read, NJ/BIONJ change the agglomeration rule,
 * BME changes the objective, and the resampling ensemble changes its noise. All of
 * them are bounded by what quartet similarity records. That matrix is largely
 * saturated - on a 200-taxon run, 14,631 of 20,100 pairs sit at or below the 1/3
 * random-quartet level and have their distance floored - so the top two thirds of it
 * carries no usable signal, and no amount of re-reading recovers what is not there.
 *
 * This measures something different. For each pair of taxa it averages the
 * <em>topological internode distance</em> - the number of edges on the path between
 * them - over the gene trees in which both appear. That is the matrix ASTRID uses.
 *
 * <h3>Additivity, and why there is no transform</h3>
 * Quartet similarity needed {@code -log s} to become additive, since NJ's consistency
 * argument assumes an additive matrix. An internode count is already additive: on a
 * tree, the path length between two leaves is exactly the sum of the unit branch
 * lengths along it. So these distances are fed to NJ unchanged, and the
 * {@link NeighborJoiningClusterer.Transform} enum is deliberately bypassed.
 *
 * <h3>Rooted input, unrooted distances</h3>
 * ASTRID's matrix is defined on unrooted gene trees, but DISCO components arrive
 * rooted, and a rooted binary tree carries a degree-2 root that is not an edge of the
 * unrooted topology. Every path crossing the root is therefore one edge too long. The
 * accumulator subtracts that edge for exactly the pairs whose LCA is the root, which
 * reproduces the unrooted path length. Without this correction the matrix is not
 * additive on any tree, and the error is not a constant offset - it applies to a
 * topology-dependent subset of pairs.
 *
 * <h3>Pairs never seen together</h3>
 * A pair of taxa can appear in no common component, in which case it has no observed
 * distance. Those entries are completed by one or more single-intermediate passes,
 * {@code D(a,b) = min_c D(a,c) + D(c,b)} over taxa c where both legs are known, which
 * is the metric closure restricted to one hop and costs O(missing * n) per pass.
 * Anything still unknown is filled with the largest observed distance - a deliberately
 * conservative "as far apart as anything observed" choice. Both counts are logged,
 * because the fill policy only matters if it ever binds, and on the 200-taxon
 * conditions measured so far it does not.
 */
public final class InternodeDistanceMatrix {

    /** Passes of single-intermediate closure before falling back to a constant fill. */
    private static final int MAX_CLOSURE_PASSES = 3;

    private final int n;
    private final double[] d;
    /** Standard error of each pair's mean, NaN where the pair was never observed. */
    private final double[] se;
    /** Gene trees in which each pair co-occurred; the real independent sample size. */
    private final int[] count;
    private final Stats stats;

    /** What the build actually saw; logged, and worth asserting on in tests. */
    public record Stats(int n, int treesUsed, long pairsObserved, long pairsMissing,
                        double meanObserved, double minObserved, double maxObserved,
                        int closurePasses, long filledByClosure, long filledByConstant,
                        long buildMs, double medianCount, double medianSe) {
        public long pairsTotal() { return (long) n * (n - 1) / 2; }
    }

    private InternodeDistanceMatrix(int n, double[] d, double[] se, int[] count, Stats stats) {
        this.n = n;
        this.d = d;
        this.se = se;
        this.count = count;
        this.stats = stats;
    }

    /** Dense n*n symmetric matrix, zero diagonal. Not copied: callers must not mutate. */
    public double[] distances() { return d; }

    public Stats stats() { return stats; }

    public int size() { return n; }

    /**
     * Standard error of each pair's mean internode distance,
     * {@code sd / sqrt(count)}, with {@code count} the number of gene trees the pair
     * co-occurred in. NaN for a pair never observed, 0 for a pair observed once or
     * with no variation.
     *
     * <p><b>Why this matters.</b> Tier 2's resampling ensemble perturbed quartet
     * similarity by a binomial standard error {@code sqrt(p(1-p)/den)} where
     * {@code den} counts quartets - about 1e7 per pair, all derived from 1000 gene
     * trees. Treating them as independent understated the error by roughly
     * {@code sqrt(1e7/1e3)} ~ 100x, which is why that ensemble needed a hand-tuned
     * scale of 100 to produce any diversity at all. Here the denominator is the gene
     * tree count itself, so the error needs no fudge factor: this is an empirical
     * standard error over independent observations, and scale 1.0 is the honest setting.
     */
    public double[] standardErrors() { return se; }

    /** Gene trees each pair co-occurred in; the sample size behind its mean. */
    public int[] counts() { return count; }

    // ── Construction ─────────────────────────────────────────────────────────

    /**
     * Accumulate the average internode distance over {@code trees}, which are the
     * DISCO single-copy components: each taxon appears at most once per tree, and
     * most trees hold a small subset of the taxa.
     */
    public static InternodeDistanceMatrix build(List<Tree> trees, int n) {
        long t0 = System.currentTimeMillis();
        if (n < 2) return null;
        long cells = (long) n * n;
        if (cells > Integer.MAX_VALUE - 8) {
            Logging.info("ASTRID matrix: n=%d needs %d cells, beyond a Java array; skipped",
                n, cells);
            return null;
        }

        double[] sum = new double[(int) cells];
        double[] sq = new double[(int) cells];
        int[] cnt = new int[(int) cells];
        int used = 0;
        for (Tree t : trees) {
            if (t == null || t.root == null) continue;
            if (accumulate(t.root, n, sum, sq, cnt)) used++;
        }

        // Normalize observed pairs; collect the statistics the fill policy needs,
        // and the per-pair standard error of the mean from the same accumulators.
        double[] seOut = new double[(int) cells];
        long observed = 0, missing = 0;
        double total = 0.0, min = Double.MAX_VALUE, max = 0.0;
        double[] obsCounts = new double[(int) ((long) n * (n - 1) / 2)];
        double[] obsSe = new double[obsCounts.length];
        int obsIdx = 0;
        for (int a = 0; a < n; a++) {
            for (int b = a + 1; b < n; b++) {
                int ab = a * n + b, ba = b * n + a;
                int k = cnt[ab];
                if (k > 0) {
                    double v = sum[ab] / k;
                    sum[ab] = v;
                    sum[ba] = v;
                    // Population variance of the per-tree distances, then SE of the mean.
                    double var = Math.max(0.0, sq[ab] / k - v * v);
                    double s = k > 1 ? Math.sqrt(var * k / (k - 1.0)) / Math.sqrt(k) : 0.0;
                    seOut[ab] = s;
                    seOut[ba] = s;
                    obsCounts[obsIdx] = k;
                    obsSe[obsIdx] = s;
                    obsIdx++;
                    observed++;
                    total += v;
                    if (v < min) min = v;
                    if (v > max) max = v;
                } else {
                    sum[ab] = Double.NaN;
                    sum[ba] = Double.NaN;
                    seOut[ab] = Double.NaN;
                    seOut[ba] = Double.NaN;
                    missing++;
                }
            }
            sum[a * n + a] = 0.0;
            seOut[a * n + a] = 0.0;
        }
        if (observed == 0) {
            Logging.info("ASTRID matrix: no pair of taxa co-occurs in any component; skipped");
            return null;
        }
        double mean = total / observed;
        if (min == Double.MAX_VALUE) min = 0.0;

        long byClosure = 0;
        int passes = 0;
        for (; missing > 0 && passes < MAX_CLOSURE_PASSES; passes++) {
            long filled = closurePass(sum, n);
            if (filled == 0) break;
            byClosure += filled;
            missing -= filled;
        }
        long byConstant = 0;
        if (missing > 0) {
            for (int i = 0; i < cells; i++) {
                if (Double.isNaN(sum[i])) { sum[i] = max; byConstant++; }
            }
            byConstant /= 2;      // counted both triangles
        }

        double medCount = median(obsCounts, obsIdx);
        double medSe = median(obsSe, obsIdx);
        Stats s = new Stats(n, used, observed, missing, mean, min, max,
            passes, byClosure, byConstant, System.currentTimeMillis() - t0,
            medCount, medSe);
        log(s);
        return new InternodeDistanceMatrix(n, sum, seOut, cnt, s);
    }

    private static double median(double[] v, int len) {
        if (len == 0) return 0.0;
        double[] c = java.util.Arrays.copyOf(v, len);
        java.util.Arrays.sort(c);
        return len % 2 == 1 ? c[len / 2] : 0.5 * (c[len / 2 - 1] + c[len / 2]);
    }

    private static void log(Stats s) {
        Logging.info("ASTRID matrix: %d component(s), %d/%d pairs observed (%.2f%%), "
                + "distance min=%.2f mean=%.2f max=%.2f, built in %d ms",
            s.treesUsed(), s.pairsObserved(), s.pairsTotal(),
            100.0 * s.pairsObserved() / s.pairsTotal(),
            s.minObserved(), s.meanObserved(), s.maxObserved(), s.buildMs());
        Logging.info("ASTRID matrix: median %.0f gene tree(s) per pair, median SE of the "
                + "mean %.4f (%.2f%% of the mean distance) - an empirical SE over "
                + "independent trees, so no scale factor is needed",
            s.medianCount(), s.medianSe(),
            s.meanObserved() == 0 ? 0.0 : 100.0 * s.medianSe() / s.meanObserved());
        if (s.filledByClosure() > 0 || s.filledByConstant() > 0) {
            Logging.info("ASTRID matrix: %d pair(s) completed by %d closure pass(es), "
                    + "%d filled with the maximum observed distance",
                s.filledByClosure(), s.closurePasses(), s.filledByConstant());
        }
        if (s.filledByConstant() > 0) {
            Logging.info("ASTRID matrix: WARNING the constant fill bound on %d pair(s); "
                + "those entries are not estimates", s.filledByConstant());
        }
    }

    /**
     * One single-intermediate metric-closure pass. Reads only entries known before the
     * pass began, so the result does not depend on iteration order.
     */
    private static long closurePass(double[] d, int n) {
        long filled = 0;
        int[] known = new int[n];
        double[] patch = null;
        int[] patchAt = null;
        int pending = 0;
        for (int a = 0; a < n; a++) {
            int kc = 0;
            for (int c = 0; c < n; c++) {
                if (c != a && !Double.isNaN(d[a * n + c])) known[kc++] = c;
            }
            if (kc == 0) continue;
            for (int b = a + 1; b < n; b++) {
                if (!Double.isNaN(d[a * n + b])) continue;
                double best = Double.NaN;
                for (int i = 0; i < kc; i++) {
                    int c = known[i];
                    if (c == b) continue;
                    double cb = d[c * n + b];
                    if (Double.isNaN(cb)) continue;
                    double v = d[a * n + c] + cb;
                    if (Double.isNaN(best) || v < best) best = v;
                }
                if (!Double.isNaN(best)) {
                    if (patch == null) { patch = new double[16]; patchAt = new int[16]; }
                    if (pending == patch.length) {
                        patch = java.util.Arrays.copyOf(patch, pending * 2);
                        patchAt = java.util.Arrays.copyOf(patchAt, pending * 2);
                    }
                    patch[pending] = best;
                    patchAt[pending] = a * n + b;
                    pending++;
                    filled++;
                }
            }
        }
        for (int i = 0; i < pending; i++) {
            int idx = patchAt[i];
            int a = idx / n, b = idx % n;
            d[a * n + b] = patch[i];
            d[b * n + a] = patch[i];
        }
        return filled;
    }

    // ── Per-tree all-pairs path lengths ──────────────────────────────────────

    /**
     * Add every leaf pair's path length in one component. Each pair is visited exactly
     * once, at its LCA, so the work is the number of pairs in the component rather than
     * the number of nodes times leaves.
     *
     * @return false when the component has fewer than two leaves and contributes nothing
     */
    private static boolean accumulate(TreeNode root, int n, double[] sum, double[] sq,
                                      int[] cnt) {
        int leaves = countLeaves(root);
        if (leaves < 2) return false;
        int[] buf = new int[leaves];        // taxon ids, grouped by subtree
        int[] depth = new int[leaves];      // edges from the component root
        int[] cursor = {0};
        gather(root, 0, buf, depth, cursor);
        cross(root, 0, n, sum, sq, cnt, buf, depth, new int[]{0}, true);
        return true;
    }

    private static int countLeaves(TreeNode node) {
        if (node.isLeaf()) return 1;
        if (node.isPolytomous()) {
            int c = 0;
            for (TreeNode ch : node.children) c += countLeaves(ch);
            return c;
        }
        return countLeaves(node.left) + countLeaves(node.right);
    }

    /** Fill {@code buf}/{@code depth} in leaf order, so every subtree owns a range. */
    private static void gather(TreeNode node, int dep, int[] buf, int[] depth, int[] cursor) {
        if (node.isLeaf()) {
            buf[cursor[0]] = node.taxonId;
            depth[cursor[0]] = dep;
            cursor[0]++;
            return;
        }
        if (node.isPolytomous()) {
            for (TreeNode ch : node.children) gather(ch, dep + 1, buf, depth, cursor);
        } else {
            gather(node.left, dep + 1, buf, depth, cursor);
            gather(node.right, dep + 1, buf, depth, cursor);
        }
    }

    /**
     * Walk the same order as {@link #gather}, and at each internal node add the pairs
     * whose LCA it is. {@code dep} is the node's depth, so a pair (u,v) meeting here is
     * {@code depth[u] + depth[v] - 2*dep} edges apart.
     *
     * <p>{@code atRoot} subtracts the degree-2 root's spurious edge: in the unrooted
     * topology the root's two children are joined directly, so a path crossing it is
     * one edge shorter than the rooted representation says.
     *
     * @return the number of leaves in this subtree, which advances the caller's cursor
     */
    private static int cross(TreeNode node, int dep, int n, double[] sum, double[] sq,
                             int[] cnt, int[] buf, int[] depth, int[] cursor,
                             boolean atRoot) {
        if (node.isLeaf()) { cursor[0]++; return 1; }

        int start = cursor[0];
        int[] childStart, childEnd;
        int kids;
        if (node.isPolytomous()) {
            kids = node.children.length;
            childStart = new int[kids];
            childEnd = new int[kids];
            for (int i = 0; i < kids; i++) {
                childStart[i] = cursor[0];
                cross(node.children[i], dep + 1, n, sum, sq, cnt, buf, depth, cursor, false);
                childEnd[i] = cursor[0];
            }
        } else {
            kids = 2;
            childStart = new int[2];
            childEnd = new int[2];
            childStart[0] = cursor[0];
            cross(node.left, dep + 1, n, sum, sq, cnt, buf, depth, cursor, false);
            childEnd[0] = cursor[0];
            childStart[1] = cursor[0];
            cross(node.right, dep + 1, n, sum, sq, cnt, buf, depth, cursor, false);
            childEnd[1] = cursor[0];
        }

        // A degree-2 root contributes no edge to the unrooted topology.
        int correction = (atRoot && kids == 2) ? 1 : 0;
        for (int i = 0; i < kids; i++) {
            for (int j = i + 1; j < kids; j++) {
                for (int p = childStart[i]; p < childEnd[i]; p++) {
                    int u = buf[p], du = depth[p];
                    for (int q = childStart[j]; q < childEnd[j]; q++) {
                        int v = buf[q];
                        int len = du + depth[q] - 2 * dep - correction;
                        int ab = u * n + v, ba = v * n + u;
                        double l2 = (double) len * len;
                        sum[ab] += len; sq[ab] += l2; cnt[ab]++;
                        sum[ba] += len; sq[ba] += l2; cnt[ba]++;
                    }
                }
            }
        }
        return cursor[0] - start;
    }

    // ── Guide trees ──────────────────────────────────────────────────────────

    /**
     * Build one guide tree from this matrix. Specs: {@code astrid} (NJ),
     * {@code astrid-bionj}, {@code astrid-fastme} (NJ seed, then BME/NNI). No
     * similarity transform is applied - see the class note on additivity.
     */
    public Tree guide(String spec, int treeIndex, int bmeMaxPasses) {
        var variant = spec.endsWith("-bionj")
            ? NeighborJoiningClusterer.Variant.BIONJ
            : NeighborJoiningClusterer.Variant.NJ;
        // Carried only for log labels; the distances are already additive.
        var label = NeighborJoiningClusterer.Transform.ONE_MINUS_SIM;
        Tree seed = NeighborJoiningClusterer.buildFromDistances(d, n, treeIndex, variant, label);
        if (seed == null || !spec.endsWith("-fastme")) return seed;
        return BalancedMEOptimizer.optimize(d, n, seed, treeIndex, label, bmeMaxPasses);
    }

    /** True for the three specs {@link #guide} accepts. */
    public static boolean isAstridSpec(String spec) {
        return spec.equals("astrid") || spec.equals("astrid-bionj")
            || spec.equals("astrid-fastme");
    }

    // ── Resampling ensemble ──────────────────────────────────────────────────

    /**
     * k guide trees from standard-error-perturbed copies of this matrix.
     *
     * <p>This is Tier 2's idea applied to a matrix that supports it properly. There,
     * the standard error came from a binomial over ~1e7 quartets per pair that were
     * really 1000 gene trees in disguise, so it was about 100x too small and the
     * ensemble collapsed to 12 distinct topologies until the scale was hand-set to
     * 100. Here {@link #standardErrors()} divides by the number of gene trees the pair
     * was actually seen in, so the perturbation {@code d' = d + scale*z*se} is on the
     * scale the data supports and {@code scale = 1.0} is the principled default.
     *
     * <p>Perturbed distances are clamped at 0: a negative distance is meaningless and
     * would make NJ's reduction formula produce nonsense rather than merely a different
     * tree. Symmetry is preserved by drawing once per unordered pair.
     *
     * <p>Members are deduplicated by exact rooted clade-set equality, so k is an upper
     * bound on the count returned. A collapse to one topology means the perturbation is
     * too small to matter and is reported as such.
     */
    public List<ResampledGuideEnsemble.Member> ensemble(String spec, int firstTreeIndex,
                                                        int k, long seed, double scale,
                                                        int bmeMaxPasses) {
        if (k < 1) return List.of();
        var variant = spec.contains("bionj")
            ? NeighborJoiningClusterer.Variant.BIONJ
            : NeighborJoiningClusterer.Variant.NJ;
        boolean bme = spec.contains("fastme");
        var label = NeighborJoiningClusterer.Transform.ONE_MINUS_SIM;

        var out = new java.util.ArrayList<ResampledGuideEnsemble.Member>(k);
        var seen = new java.util.HashSet<java.util.Set<java.util.BitSet>>();
        double[] perturbed = new double[n * n];
        long zeroSe = 0;
        for (int i = 0; i < n * n; i++) if (se[i] == 0.0) zeroSe++;

        for (int r = 0; r < k; r++) {
            var rng = new java.util.SplittableRandom(seed + r);
            for (int a = 0; a < n; a++) perturbed[a * n + a] = 0.0;
            for (int a = 0; a < n; a++) {
                for (int b = a + 1; b < n; b++) {
                    int ab = a * n + b;
                    double s = se[ab];
                    double v = d[ab];
                    if (!Double.isNaN(s) && s > 0.0) {
                        v = Math.max(0.0, v + scale * s * gaussian(rng));
                    }
                    perturbed[ab] = v;
                    perturbed[b * n + a] = v;
                }
            }
            Tree tree = NeighborJoiningClusterer.buildFromDistances(
                perturbed, n, firstTreeIndex + out.size(), variant, label);
            if (tree == null) continue;
            if (bme) {
                tree = BalancedMEOptimizer.optimize(
                    perturbed, n, tree, firstTreeIndex + out.size(), label, bmeMaxPasses);
            }
            if (!seen.add(ResampledGuideEnsemble.cladeSet(tree, n))) continue;
            out.add(new ResampledGuideEnsemble.Member(
                String.format("%s~e%02d", spec, r + 1), tree));
        }

        // Union and mean Jaccard, reported in the same terms as the Tier 2 ensemble so
        // the two are directly comparable: diversity is the only reason to build k trees.
        var cladeSets = new java.util.ArrayList<java.util.Set<java.util.BitSet>>(out.size());
        for (var m : out) cladeSets.add(ResampledGuideEnsemble.cladeSet(m.tree(), n));
        var union = new java.util.HashSet<java.util.BitSet>();
        for (var s : cladeSets) union.addAll(s);
        double jac = 0.0;
        int pairsCompared = 0;
        for (int a = 0; a < cladeSets.size(); a++) {
            for (int b = a + 1; b < cladeSets.size(); b++) {
                var x = cladeSets.get(a);
                var y = cladeSets.get(b);
                var inter = new java.util.HashSet<>(x);
                inter.retainAll(y);
                var uni = new java.util.HashSet<>(x);
                uni.addAll(y);
                if (!uni.isEmpty()) { jac += (double) inter.size() / uni.size(); pairsCompared++; }
            }
        }
        int perMember = cladeSets.isEmpty() ? 0 : cladeSets.get(0).size();
        Logging.info("ASTRID ensemble (%s): k=%d requested, %d distinct topologies, "
                + "seed=%d, se scale=%.2f, %d pair(s) had se=0 and were not perturbed",
            spec, k, out.size(), seed, scale, zeroSe / 2);
        Logging.info("ASTRID ensemble (%s): %d clades per member, %d distinct in the "
                + "union (%.2fx one member), mean pairwise clade Jaccard %.3f",
            spec, perMember, union.size(),
            perMember == 0 ? 0.0 : (double) union.size() / perMember,
            pairsCompared == 0 ? 1.0 : jac / pairsCompared);
        if (out.size() == 1 && k > 1) {
            Logging.info("ASTRID ensemble (%s): WARNING collapsed to one topology - the "
                + "standard error is too small at this scale to change any join", spec);
        }
        return out;
    }

    /** Box-Muller, since SplittableRandom has no Gaussian of its own. */
    private static double gaussian(java.util.SplittableRandom rng) {
        double u1 = Math.max(rng.nextDouble(), 1e-300);
        double u2 = rng.nextDouble();
        return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    }

    /** True for an {@code astrid-ensemble[-bionj|-fastme]:k[:seed[:scale]]} spec. */
    public static boolean isAstridEnsembleSpec(String spec) {
        return spec.matches("astrid-ensemble(-bionj|-fastme)?:\\d{1,4}(:\\d{1,19}(:\\d+(\\.\\d+)?)?)?");
    }

    // ── Self-test hooks ──────────────────────────────────────────────────────

    /** Exposed for the self-test: all-pairs path lengths of a single rooted tree. */
    static double[] singleTreeForTest(TreeNode root, int n) {
        double[] sum = new double[n * n];
        double[] sq = new double[n * n];
        int[] cnt = new int[n * n];
        accumulate(root, n, sum, sq, cnt);
        for (int i = 0; i < n * n; i++) if (cnt[i] > 0) sum[i] /= cnt[i];
        return sum;
    }

    /** Exposed for the self-test: one closure pass over a matrix holding NaN gaps. */
    static long closurePassForTest(double[] d, int n) { return closurePass(d, n); }
}
