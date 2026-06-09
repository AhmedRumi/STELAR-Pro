package core;

import java.util.*;
import tree.Tree;
import tree.TreeNode;

/**
 * Manages inverse permutation arrays for efficient range intersection calculations.
 * 
 * This class builds and maintains:
 * 1. Gene tree orderings: [treeIndex][position] = taxonId
 * 2. Position index: [treeIndex][taxonId] = sorted positions of that taxon
 * 
 * Enables O(min(|A|, |B|)) intersection counting instead of O(n) BitSet operations.
 * 
 * The intersection algorithm works by:
 * 1. Choosing the smaller of two ranges for iteration
 * 2. For each unique taxon in the smaller range, binary-searching that taxon's
 *    position array in the other tree to determine whether it occurs in the
 *    target range
 */
public class InverseIndexManager {
    
    // Legacy single-position view retained for the existing GPU/JNA path.
    // STELAR-Pro CPU scoring uses taxonPositionsByTree instead.
    private final int[][] inverseIndex; // [treeIndex][taxonId] = first position, or -1
    private final int[][][] taxonPositionsByTree; // [treeIndex][taxonId] = sorted positions
    private final int[][] geneTreeOrderings; // [treeIndex][position] = taxonId
    private final int numTrees;
    private final int numTaxa;
    
    // Statistics for logging
    private long totalIntersectionCalls = 0;
    private long totalElementsProcessed = 0;
    private long maxRangeSize = 0;
    private long minRangeSize = Long.MAX_VALUE;
    
    public InverseIndexManager(List<Tree> geneTrees, int realTaxaCount) {
        System.out.println("==== INITIALIZING INVERSE INDEX MANAGER ====");
        System.out.println("Number of gene trees: " + geneTrees.size());
        System.out.println("Number of taxa: " + realTaxaCount);
        
        this.numTrees = geneTrees.size();
        this.numTaxa = realTaxaCount;
        this.inverseIndex = new int[numTrees][numTaxa];
        this.taxonPositionsByTree = new int[numTrees][numTaxa][];
        this.geneTreeOrderings = new int[numTrees][];
        
        long startTime = System.currentTimeMillis();
        buildInverseIndex(geneTrees);
        long endTime = System.currentTimeMillis();
        
        System.out.println("Inverse index construction completed in " + (endTime - startTime) + " ms");
        System.out.println("Memory allocated: " + 
                         estimateIndexMemoryBytes() / (1024 * 1024) + " MB for index arrays");
        System.out.println("==== INVERSE INDEX MANAGER READY ====");
    }
    
    /**
     * Build the per-tree position index from left-to-right leaf orderings.
     *
     * STELAR-X stored one position per taxon because every taxon appeared at
     * most once in a gene tree. STELAR-Pro gene family trees can contain several
     * copies from the same species, so one position is not enough. For each tree
     * and each taxon/species, we store every occurrence position in a sorted int
     * array:
     *
     *   taxonPositionsByTree[treeIdx][taxonId] = {p1, p2, p3, ...}
     *
     * This lets range intersections ask: "does this species occur anywhere in
     * [start, end) of the other tree?" via binary search.
     */
    private void buildInverseIndex(List<Tree> geneTrees) {
        System.out.println("Building inverse index mappings...");
        
        System.out.println("Initializing legacy inverse index with sentinel values (-1 for non-existent taxa)...");
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            java.util.Arrays.fill(inverseIndex[treeIdx], -1);
            for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
                taxonPositionsByTree[treeIdx][taxonId] = new int[0];
            }
        }
        
        int processedTrees = 0;
        int totalLeaves = 0;
        int totalAbsentTaxa = 0;
        
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            Tree tree = geneTrees.get(treeIdx);
            
            // Get left-to-right ordering (same as MemoryEfficientBipartitionManager)
            List<Integer> ordering = new ArrayList<>();
            collectLeavesInOrder(tree.root, ordering);
            
            geneTreeOrderings[treeIdx] = ordering.stream().mapToInt(Integer::intValue).toArray();
            totalLeaves += ordering.size();

            @SuppressWarnings("unchecked")
            List<Integer>[] positionsByTaxon = new ArrayList[numTaxa];
            for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
                positionsByTaxon[taxonId] = new ArrayList<>();
            }

            for (int pos = 0; pos < geneTreeOrderings[treeIdx].length; pos++) {
                int taxonId = geneTreeOrderings[treeIdx][pos];
                if (taxonId >= 0 && taxonId < numTaxa) {
                    positionsByTaxon[taxonId].add(pos);

                    // Legacy single-position view: keep the first occurrence.
                    // Correct STELAR-Pro intersections use taxonPositionsByTree.
                    if (inverseIndex[treeIdx][taxonId] == -1) {
                        inverseIndex[treeIdx][taxonId] = pos;
                    }
                } else {
                    System.err.println("WARNING: Invalid taxon ID " + taxonId + 
                                     " in tree " + treeIdx + " at position " + pos);
                }
            }

            int absentTaxaInTree = 0;
            for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
                List<Integer> positions = positionsByTaxon[taxonId];
                int[] compactPositions = new int[positions.size()];
                for (int i = 0; i < positions.size(); i++) {
                    compactPositions[i] = positions.get(i);
                }
                taxonPositionsByTree[treeIdx][taxonId] = compactPositions;

                if (compactPositions.length == 0) {
                    absentTaxaInTree++;
                }
            }
            totalAbsentTaxa += absentTaxaInTree;
            
            processedTrees++;
            
            // Log progress for large datasets (uncomment if needed for debugging)
            if (processedTrees % 100 == 0 || processedTrees == numTrees) {
                System.out.println("Processed " + processedTrees + "/" + numTrees + 
                                 " trees, average leaves per tree: " + 
                                 (totalLeaves / (double) processedTrees));
            }
        }
        
        System.out.println("Inverse index built successfully");
        System.out.println("Total leaves processed: " + totalLeaves);
        System.out.println("Average leaves per tree: " + (totalLeaves / (double) numTrees));
        System.out.println("Total absent taxon entries: " + totalAbsentTaxa);
        System.out.println("Average missing taxa per tree: " + (totalAbsentTaxa / (double) numTrees));
        
        validateInverseIndex();
    }
    
    /**
     * Collect leaves in left-to-right order (inorder traversal).
     * Same implementation as MemoryEfficientBipartitionManager for consistency.
     */
    private void collectLeavesInOrder(TreeNode node, List<Integer> ordering) {
        if (node.isLeaf()) {
            ordering.add(node.taxon.id);
            return;
        }
        
        // For binary trees, visit left child, then right child
        if (node.childs != null && node.childs.size() >= 2) {
            collectLeavesInOrder(node.childs.get(0), ordering);
            collectLeavesInOrder(node.childs.get(1), ordering);
        }
        
        // Handle any additional children (though binary trees should only have 2)
        for (int i = 2; i < (node.childs != null ? node.childs.size() : 0); i++) {
            collectLeavesInOrder(node.childs.get(i), ordering);
        }
    }
    
    /**
     * Validate that every taxon position vector matches the gene tree ordering.
     *
     * Duplicates are now valid. For example, if taxon 20 occurs at positions
     * 4, 7, and 10, taxonPositionsByTree[tree][20] must be exactly {4, 7, 10}.
     * Taxa absent from the tree must have an empty vector and a legacy sentinel
     * of -1.
     */
    private void validateInverseIndex() {
        System.out.println("Validating inverse index consistency (position vectors and sentinel values)...");
        
        int validationErrors = 0;
        int sentinelValidationErrors = 0;
        
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            int[] ordering = geneTreeOrderings[treeIdx];
            
            for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
                List<Integer> expectedPositions = new ArrayList<>();
                for (int pos = 0; pos < ordering.length; pos++) {
                    if (ordering[pos] == taxonId) {
                        expectedPositions.add(pos);
                    }
                }

                int[] actualPositions = taxonPositionsByTree[treeIdx][taxonId];
                if (actualPositions.length != expectedPositions.size()) {
                    System.err.println("Position-vector validation error: tree " + treeIdx
                            + ", taxon " + taxonId
                            + ", expected " + expectedPositions.size()
                            + " positions, got " + actualPositions.length);
                    validationErrors++;
                } else {
                    for (int i = 0; i < actualPositions.length; i++) {
                        if (actualPositions[i] != expectedPositions.get(i)) {
                            System.err.println("Position-vector validation error: tree " + treeIdx
                                    + ", taxon " + taxonId
                                    + ", index " + i
                                    + ", expected pos " + expectedPositions.get(i)
                                    + ", got pos " + actualPositions[i]);
                            validationErrors++;
                            break;
                        }
                    }
                }

                if (expectedPositions.isEmpty() && inverseIndex[treeIdx][taxonId] != -1) {
                    System.err.println("Sentinel validation error: tree " + treeIdx
                            + ", taxon " + taxonId
                            + " absent but legacy inverse index is " + inverseIndex[treeIdx][taxonId]);
                    sentinelValidationErrors++;
                }

                if (validationErrors >= 10 || sentinelValidationErrors >= 10) {
                    System.err.println("Too many inverse-index validation errors, stopping validation");
                    break;
                }
            }
            
            if (validationErrors >= 10 || sentinelValidationErrors >= 10) break;
        }
        
        if (validationErrors == 0 && sentinelValidationErrors == 0) {
            System.out.println("Inverse index validation passed (position vectors and sentinel values)");
        } else {
            System.err.println("Inverse index validation failed: " + 
                             validationErrors + " position-vector errors, " + 
                             sentinelValidationErrors + " sentinel errors");
        }
    }
    
    /**
     * Calculate intersection size between two ranges using inverse index.
     * 
     * STELAR-Pro behavior with duplicated taxa:
     * - Each species is counted at most once in the intersection.
     * - The smaller range is scanned by occurrence positions.
     * - A small-range HashSet suppresses repeated copies of the same species.
     * - Membership in the other range is checked by binary search over that
     *   species' sorted position vector in the other tree.
     * 
     * Complexity: O(min(|range1|, |range2|)) instead of O(n) for BitSet operations.
     * 
     * @param tree1 First tree index
     * @param start1 Start position in first tree (inclusive)
     * @param end1 End position in first tree (exclusive)
     * @param tree2 Second tree index
     * @param start2 Start position in second tree (inclusive)
     * @param end2 End position in second tree (exclusive)
     * @return Number of taxa in the intersection of the two ranges
     */
    public int getRangeIntersectionSize(int tree1, int start1, int end1, 
                                       int tree2, int start2, int end2) {
        // Validate input parameters
        if (tree1 < 0 || tree1 >= numTrees || tree2 < 0 || tree2 >= numTrees) {
            System.err.println("Invalid tree indices: tree1=" + tree1 + ", tree2=" + tree2);
            return 0;
        }
        
        if (start1 < 0 || end1 < start1 || start2 < 0 || end2 < start2) {
            System.err.println("Invalid range parameters: [" + start1 + "," + end1 + ") and [" + start2 + "," + end2 + ")");
            return 0;
        }
        
        // Additional validation for tree orderings
        if (geneTreeOrderings[tree1] == null || geneTreeOrderings[tree2] == null) {
            System.err.println("Null gene tree orderings for tree1=" + tree1 + " or tree2=" + tree2);
            return 0;
        }
        
        // Update statistics
        totalIntersectionCalls++;
        int size1 = end1 - start1;
        int size2 = end2 - start2;
        maxRangeSize = Math.max(maxRangeSize, Math.max(size1, size2));
        minRangeSize = Math.min(minRangeSize, Math.min(size1, size2));
        
        // Choose smaller range for iteration (complexity optimization)
        // This is especially important when trees have different taxa counts
        if (size1 <= size2) {
            totalElementsProcessed += size1;
            return countIntersection(tree1, start1, end1, tree2, start2, end2);
        } else {
            totalElementsProcessed += size2;
            return countIntersection(tree2, start2, end2, tree1, start1, end1);
        }
    }
    
    /**
     * Count an intersection by scanning the smaller occurrence range.
     *
     * The input ranges are still ranges of gene-tree leaf occurrences. With
     * paralogs, a range can contain the same species multiple times. The triplet
     * scoring model works on species sets, so we count each taxon/species once:
     *
     *   small range: [20, 20, 21, 20]
     *   target range contains 20 and 21
     *   intersection count = 2, not 4
     *
     * For each first-seen taxon in the smaller range, the position vector in the
     * larger tree is binary-searched. If the first occurrence >= largeStart is
     * still < largeEnd, that species intersects the target range.
     * 
     * @param smallTree Tree index for the smaller range
     * @param smallStart Start position in smaller range
     * @param smallEnd End position in smaller range
     * @param largeTree Tree index for the larger range
     * @param largeStart Start position in larger range
     * @param largeEnd End position in larger range
     * @return Intersection count
     */
    private int countIntersection(int smallTree, int smallStart, int smallEnd,
                                 int largeTree, int largeStart, int largeEnd) {
        int count = 0;
        int[] smallOrdering = geneTreeOrderings[smallTree];
        
        // Validate array bounds
        if (smallOrdering == null) {
            System.err.println("Null ordering for tree " + smallTree);
            return 0;
        }
        
        int maxPos = Math.min(smallEnd, smallOrdering.length);
        Set<Integer> seenInSmallRange = new HashSet<>();
        
        for (int pos = smallStart; pos < maxPos; pos++) {
            int taxonId = smallOrdering[pos];
            
            // Validate taxon ID
            if (taxonId < 0 || taxonId >= numTaxa) {
                System.err.println("Invalid taxon ID " + taxonId + " at position " + pos + " in tree " + smallTree);
                continue;
            }

            if (!seenInSmallRange.add(taxonId)) {
                continue;
            }
            
            if (hasPositionInRange(largeTree, taxonId, largeStart, largeEnd)) {
                count++;
            }
        }

        return count;
    }

    /**
     * Return true if a taxon has at least one occurrence in [start, end).
     *
     * The position array is sorted because it is collected from the left-to-right
     * gene tree ordering. Arrays.binarySearch returns either the matching index
     * or the insertion point for start; in both cases, checking the candidate
     * index tells us whether an occurrence lies inside the target range.
     */
    private boolean hasPositionInRange(int treeIdx, int taxonId, int start, int end) {
        if (treeIdx < 0 || treeIdx >= numTrees || taxonId < 0 || taxonId >= numTaxa) {
            return false;
        }

        int[] positions = taxonPositionsByTree[treeIdx][taxonId];
        if (positions == null || positions.length == 0) {
            return false;
        }

        int idx = Arrays.binarySearch(positions, start);
        if (idx < 0) {
            idx = -idx - 1;
        }
        return idx < positions.length && positions[idx] < end;
    }
    
    /**
     * Get statistics about intersection calculations for performance monitoring.
     * 
     * ENHANCED: Includes information about sentinel value handling.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Inverse Index Manager Statistics:\n");
        sb.append("  Trees: ").append(numTrees).append("\n");
        sb.append("  Taxa: ").append(numTaxa).append("\n");
        sb.append("  Total intersection calls: ").append(totalIntersectionCalls).append("\n");
        sb.append("  Total elements processed: ").append(totalElementsProcessed).append("\n");
        
        if (totalIntersectionCalls > 0) {
            sb.append("  Average elements per call: ").append(totalElementsProcessed / (double) totalIntersectionCalls).append("\n");
            sb.append("  Min range size: ").append(minRangeSize == Long.MAX_VALUE ? 0 : minRangeSize).append("\n");
            sb.append("  Max range size: ").append(maxRangeSize).append("\n");
        }
        
        int absentTaxonEntries = 0;
        int totalStoredPositions = 0;
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
                int[] positions = taxonPositionsByTree[treeIdx][taxonId];
                if (positions == null || positions.length == 0) {
                    absentTaxonEntries++;
                } else {
                    totalStoredPositions += positions.length;
                }
            }
        }
        
        sb.append("  Stored taxon occurrence positions: ").append(totalStoredPositions).append("\n");
        sb.append("  Absent taxon entries: ").append(absentTaxonEntries).append("\n");
        sb.append("  Average missing taxa per tree: ").append(absentTaxonEntries / (double) numTrees).append("\n");
        sb.append("  Taxa coverage: ").append(String.format("%.2f%%", 
                 100.0 * (numTrees * numTaxa - absentTaxonEntries) / (double)(numTrees * numTaxa))).append("\n");
        
        return sb.toString();
    }
    
    /**
     * Reset statistics counters.
     */
    public void resetStatistics() {
        totalIntersectionCalls = 0;
        totalElementsProcessed = 0;
        maxRangeSize = 0;
        minRangeSize = Long.MAX_VALUE;
    }
    
    // Getters for external access (e.g., GPU implementation).
    // NOTE: getInverseIndex() is a legacy single-position view. It is not
    // sufficient for STELAR-Pro duplicated taxa; CPU scoring uses
    // taxonPositionsByTree through getRangeIntersectionSize().
    public int[][] getInverseIndex() { 
        return inverseIndex; 
    }
    
    public int[][] getGeneTreeOrderings() { 
        return geneTreeOrderings; 
    }
    
    public int getNumTrees() { 
        return numTrees; 
    }
    
    public int getNumTaxa() { 
        return numTaxa; 
    }

    private long estimateIndexMemoryBytes() {
        long bytes = (long) numTrees * numTaxa * 4L; // legacy inverseIndex
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
                int[] positions = taxonPositionsByTree[treeIdx][taxonId];
                if (positions != null) {
                    bytes += (long) positions.length * 4L;
                }
            }
        }
        return bytes;
    }
}
