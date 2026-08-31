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

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Reports occurrence-sensitive and duplicate-invariant S1 candidate counts. */
public final class CandidateCountReport {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected tagged gene trees");

        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(args[0], registry, false);
        TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 0xDEADBEEFCAFEL);
        PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);

        ClusterTable previousClusters = new ClusterTable(trees, pref, registry.size());
        DPTable previousDP = new DPTable(trees, pref, previousClusters);

        UniqueTaxonSubtreeHashes unique = new UniqueTaxonSubtreeHashes(trees, hasher);
        ClusterTable currentClusters = new ClusterTable(trees, pref, registry.size(), unique);
        DPTable currentDP = new DPTable(trees, pref, currentClusters, unique);

        long leaves = 0;
        long speciation = 0;
        long duplication = 0;
        Set<BipartitionSplit> allOccurrenceSplits = new HashSet<>();
        Set<BipartitionSplit> speciationOccurrenceSplits = new HashSet<>();
        Set<BipartitionSplit> speciationSpeciesSplits = new HashSet<>();
        Set<BipartitionSplit> legacyAllSplits = new HashSet<>();
        Set<BipartitionSplit> legacySpeciationSplits = new HashSet<>();
        for (Tree tree : trees) {
            leaves += tree.leafCount;
            int[] firstPosition = new int[registry.size()];
            Arrays.fill(firstPosition, -1);
            for (int position = 0; position < tree.leafCount; position++) {
                int taxon = tree.postorderArray[position];
                if (firstPosition[taxon] < 0) firstPosition[taxon] = position;
            }
            Map<TreeNode, int[]> legacyRanges = new IdentityHashMap<>();
            buildLegacyRanges(tree.root, firstPosition, legacyRanges);

            ArrayDeque<TreeNode> nodes = new ArrayDeque<>();
            nodes.push(tree.root);
            while (!nodes.isEmpty()) {
                TreeNode node = nodes.pop();
                if (node.isLeaf()) continue;
                if (node.isDuplication()) duplication++;
                else if (node.isSpeciation()) speciation++;
                if (node.isPolytomous()) {
                    for (TreeNode child : node.children) nodes.push(child);
                } else {
                    BipartitionSplit occurrenceSplit = new BipartitionSplit(
                        occurrenceHash(pref, tree.treeIndex, node.left),
                        occurrenceHash(pref, tree.treeIndex, node.right));
                    allOccurrenceSplits.add(occurrenceSplit);
                    int[] legacyLeft = legacyRanges.get(node.left);
                    int[] legacyRight = legacyRanges.get(node.right);
                    BipartitionSplit legacySplit = new BipartitionSplit(
                        occurrenceHash(pref, tree.treeIndex, legacyLeft[0], legacyLeft[1]),
                        occurrenceHash(pref, tree.treeIndex, legacyRight[0], legacyRight[1]));
                    legacyAllSplits.add(legacySplit);
                    if (node.isSpeciation()) {
                        speciationOccurrenceSplits.add(occurrenceSplit);
                        legacySpeciationSplits.add(legacySplit);
                        speciationSpeciesSplits.add(new BipartitionSplit(
                            unique.get(tree.treeIndex, node.left),
                            unique.get(tree.treeIndex, node.right)));
                    }
                    nodes.push(node.left);
                    nodes.push(node.right);
                }
            }
        }

        System.out.printf("trees=%d taxa=%d leaf_occurrences=%d%n",
            trees.size(), registry.size(), leaves);
        System.out.printf("speciation_nodes=%d duplication_nodes=%d%n",
            speciation, duplication);
        System.out.printf("all_internal_occurrence_sensitive_splits=%d%n",
            allOccurrenceSplits.size());
        System.out.printf("speciation_occurrence_sensitive_splits=%d%n",
            speciationOccurrenceSplits.size());
        System.out.printf("speciation_duplicate_collapsed_splits=%d%n",
            speciationSpeciesSplits.size());
        System.out.printf("historical_first_copy_all_splits=%d%n", legacyAllSplits.size());
        System.out.printf("historical_first_copy_speciation_splits=%d%n",
            legacySpeciationSplits.size());
        System.out.printf("previous_unique_subtree_bipartitions=%d%n", previousClusters.size());
        System.out.printf("current_unique_subtree_bipartitions=%d%n", currentClusters.size());
        System.out.printf("previous_nontrivial_subtree_bipartitions=%d%n",
            previousClusters.entries().stream().filter(entry -> entry.hash.size > 1).count());
        System.out.printf("current_nontrivial_subtree_bipartitions=%d%n",
            currentClusters.entries().stream().filter(entry -> entry.hash.size > 1).count());
        System.out.printf("previous_unique_dp_splits=%d%n", previousDP.numUniqueSplits());
        System.out.printf("current_unique_dp_splits=%d%n", currentDP.numUniqueSplits());
        System.out.printf("current_emitted_speciation_splits=%d%n", currentDP.numEmitted());
        System.out.printf("overlapping_tagged_speciations_skipped=%d%n",
            currentDP.numOverlappingSpeciationNodesSkipped());
    }

    private static ClusterHash occurrenceHash(PrefixHashArrays pref, int treeIndex,
                                               TreeNode node) {
        return occurrenceHash(pref, treeIndex, node.rangeStart, node.rangeEnd);
    }

    private static ClusterHash occurrenceHash(PrefixHashArrays pref, int treeIndex,
                                               int start, int end) {
        int m = pref.numSeeds();
        long[] sums = new long[m];
        long[] xors = new long[m];
        for (int seed = 0; seed < m; seed++) {
            sums[seed] = pref.rangeSum(treeIndex, seed, start, end);
            xors[seed] = pref.rangeXor(treeIndex, seed, start, end);
        }
        return new ClusterHash(sums, xors, end - start, m);
    }

    /** Reproduce the former range builder, which mapped every copy to its first position. */
    private static int[] buildLegacyRanges(TreeNode node, int[] firstPosition,
                                           Map<TreeNode, int[]> ranges) {
        int[] range;
        if (node.isLeaf()) {
            int position = firstPosition[node.taxonId];
            range = new int[]{position, position + 1};
        } else {
            int start = Integer.MAX_VALUE;
            int end = Integer.MIN_VALUE;
            if (node.isPolytomous()) {
                for (TreeNode child : node.children) {
                    int[] childRange = buildLegacyRanges(child, firstPosition, ranges);
                    start = Math.min(start, childRange[0]);
                    end = Math.max(end, childRange[1]);
                }
            } else {
                int[] left = buildLegacyRanges(node.left, firstPosition, ranges);
                int[] right = buildLegacyRanges(node.right, firstPosition, ranges);
                start = Math.min(left[0], right[0]);
                end = Math.max(left[1], right[1]);
            }
            range = new int[]{start, end};
        }
        ranges.put(node, range);
        return range;
    }
}
