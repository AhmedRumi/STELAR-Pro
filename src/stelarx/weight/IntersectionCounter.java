package stelarx.weight;

import stelarx.tree.Tree;

/**
 * Duplicate-invariant intersections between postorder ranges.
 *
 * For complete gene trees (Lg = L), complement intersections reduce to:
 *   |comp(P) ∩ Q| = |Q| - |P ∩ Q|
 *   |P ∩ comp(Q)| = |P| - |P ∩ Q|
 *   |comp(P) ∩ comp(Q)| = n - |P| - |Q| + |P ∩ Q|
 *
 * So we only ever need the one "core" non-complement count.
 */
public final class IntersectionCounter {

    private IntersectionCounter() {}

    /**
     * |range_in_treeA ∩ range_in_treeB| -- both ranges non-complement.
     * Iterates over the smaller occurrence range. A taxon is visited only at its
     * first copy in that range, then its sorted position vector in the other tree
     * is binary-searched for membership.
     */
    public static int coreIntersect(Tree tA, int loA, int hiA,
                                     Tree tB, int loB, int hiB) {
        if (hiA - loA <= hiB - loB) {
            return scanUniqueTaxa(tA, loA, hiA, tB, loB, hiB);
        }
        return scanUniqueTaxa(tB, loB, hiB, tA, loA, hiA);
    }

    private static int scanUniqueTaxa(Tree source, int sourceLo, int sourceHi,
                                      Tree target, int targetLo, int targetHi) {
        int count = 0;
        for (int pos = sourceLo; pos < sourceHi; pos++) {
            int taxon = source.postorderArray[pos];
            if (!source.taxonPositions.isFirstInRange(taxon, pos, sourceLo)) continue;
            if (target.taxonPositions.containsInRange(taxon, targetLo, targetHi)) count++;
        }
        return count;
    }

    /**
     * |M_range ∩ cluster| where M_range is non-complement (a gene-tree subtree
     * range) and cluster may be complement (a species-tree candidate part).
     *
     * For complete trees: |comp(cluster_range) ∩ M_range| = |M_range| - |cluster_range ∩ M_range|
     *
     * @param tGT       gene-tree tree
     * @param loGT, hiGT  range in gene-tree postorder array (non-complement)
     * @param tC        exemplar tree of the candidate cluster
     * @param loC, hiC  range in candidate cluster's exemplar tree
     * @param cComp     whether the candidate cluster is complement w.r.t. its tree
     * @param sizeGTRange number of distinct taxa in the gene-tree range
     */
    public static int intersect(Tree tGT, int loGT, int hiGT,
                                 Tree tC, int loC, int hiC, boolean cComp,
                                 int sizeGTRange) {
        int core = coreIntersect(tGT, loGT, hiGT, tC, loC, hiC);
        return cComp ? (sizeGTRange - core) : core;
    }

    /**
     * |A ∩ Lg_GT| — how many of the gene-tree's leaves fall inside cluster A.
     *
     * For a super-complement cluster A = S\sub(u):
     *   |A ∩ Lg_GT| = L_GT - |sub(u) ∩ Lg_GT|  (since M ⊆ Lg_GT ⊆ S)
     * For a plain cluster A = sub(u):
     *   |A ∩ Lg_GT| = |sub(u) ∩ Lg_GT|
     *
     * Used to compute the correct row sum for the A-row of the intersection matrix
     * when tGT is an incomplete gene tree (L_GT < n).
     * For complete gene trees (L_GT == n) this returns sizeA directly.
     *
     * @param tGT   gene tree
     * @param tC    exemplar tree for cluster A
     * @param loC, hiC  range of A in tC's postorder array (the un-complemented range)
     * @param cComp true if A is a complement cluster (A = S\[loC,hiC))
     */
    public static int intersectWithFullTree(Tree tGT, Tree tC, int loC, int hiC, boolean cComp) {
        int L_GT = tGT.distinctTaxonCount;
        int core = coreIntersect(tGT, 0, tGT.leafCount, tC, loC, hiC);
        return cComp ? (L_GT - core) : core;
    }

    // -------------------------------------------------------------------------
    // Multi-range cluster variants (DOCS/multi-range-cluster-design.md §5.1).
    //
    // A multi-range cluster's positive part is a union of PAIRWISE-DISJOINT
    // ranges {[los[j],his[j])} in its exemplar tree tC.  Because the ranges are
    // disjoint, |M ∩ (⋃ ranges)| = Σ_j |M ∩ [los[j],his[j])| — so each core count
    // is a sum of per-range coreIntersect()s (each still walks the smaller side).
    // The complement subtract trick (cComp ? size−core : core) is unchanged:
    // |comp ∩ M| = |M| − |(⋃ ranges) ∩ M|, valid since M ⊆ Lg ⊆ S.
    //
    // These mirror the single-range methods exactly when los.length == 1, so a
    // caller may always dispatch on Cluster.isMultiRange() with identical results.
    // -------------------------------------------------------------------------

    /** |[loGT,hiGT) ∩ union(ranges)|, counting each taxon once across all ranges. */
    public static int coreIntersectMulti(Tree tGT, int loGT, int hiGT,
                                          Tree tC, int[] los, int[] his) {
        int count = 0;
        for (int pos = loGT; pos < hiGT; pos++) {
            int taxon = tGT.postorderArray[pos];
            if (!tGT.taxonPositions.isFirstInRange(taxon, pos, loGT)) continue;
            for (int j = 0; j < los.length; j++) {
                if (tC.taxonPositions.containsInRange(taxon, los[j], his[j])) {
                    count++;
                    break;
                }
            }
        }
        return count;
    }

    /** Multi-range analogue of {@link #intersect}. */
    public static int intersectMulti(Tree tGT, int loGT, int hiGT,
                                     Tree tC, int[] los, int[] his, boolean cComp,
                                     int sizeGTRange) {
        int core = coreIntersectMulti(tGT, loGT, hiGT, tC, los, his);
        return cComp ? (sizeGTRange - core) : core;
    }

    /** Multi-range analogue of {@link #intersectWithFullTree} (row sum vs full gene tree). */
    public static int intersectWithFullTreeMulti(Tree tGT, Tree tC,
                                                 int[] los, int[] his, boolean cComp) {
        int L_GT = tGT.distinctTaxonCount;
        int core = coreIntersectMulti(tGT, 0, tGT.leafCount, tC, los, his);
        return cComp ? (L_GT - core) : core;
    }
}
