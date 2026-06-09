package core;

import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicLong;

import com.sun.jna.Memory;

import preprocessing.GeneTrees;
import tree.RangeBipartition;
import tree.Tree;
import tree.TreeNode;
import utils.Config;
import utils.Threading;
import core.WeightCalculator.WeightCalcLib;
import core.WeightCalculator.WeightCalcLib.CompactBipartition;

/**
 * Efficient triplet score calculator for species tree vs gene trees.
 * 
 * Uses the same range-based bipartition representation as gene trees.
 * The species tree is treated as an additional tree in the inverse index,
 * allowing O(min(|A|, |B|)) intersection calculations.
 * 
 * Supports CPU_SINGLE, CPU_PARALLEL, and GPU_PARALLEL modes.
 * No BitSets used - everything is range-based for efficiency.
 */
public class SpeciesTreeScorer {
    
    private final GeneTrees geneTrees;
    private final int numGeneTrees;
    private final int numTaxa;
    
    // Extended inverse index: includes gene trees + species tree.
    // Species tree is at index = numGeneTrees.
    //
    // extendedInverseIndex is a legacy first-position view used only for simple
    // species-tree leaf range extraction. STELAR-Pro scoring uses
    // extendedTaxonPositions so duplicated gene-tree taxa are not collapsed to
    // one arbitrary occurrence.
    private int[][] extendedInverseIndex;
    private int[][][] extendedTaxonPositions;
    private int[][] extendedOrderings;
    
    // Species tree bipartitions as RangeBipartitions
    private List<RangeBipartition> speciesBipartitions;
    
    public SpeciesTreeScorer(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
        this.numGeneTrees = geneTrees.geneTrees.size();
        this.numTaxa = geneTrees.realTaxaCount;
    }
    
    /**
     * Calculate triplet score between species tree and gene trees.
     * 
     * @param speciesTree The species tree to score
     * @return Total triplet score
     */
    public double calculateScore(Tree speciesTree) {
        System.out.println("\n=== Species Tree Triplet Score Calculation ===");
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Gene trees: " + numGeneTrees);
        System.out.println("Gene tree bipartitions: " + geneTrees.rangeBipartitions.size());
        System.out.println("Taxa: " + numTaxa);
        
        long startTime = System.currentTimeMillis();
        
        // Step 1: Build species tree ordering and add to extended inverse index
        buildExtendedInverseIndex(speciesTree);
        
        // Step 2: Extract bipartitions from species tree as RangeBipartitions
        extractSpeciesBipartitions(speciesTree);
        System.out.println("Species tree bipartitions: " + speciesBipartitions.size());
        
        // Step 3: Calculate total score based on computation mode
        double totalScore;
        
        switch (Config.COMPUTATION_MODE) {
            case GPU_PARALLEL:
                System.out.println("Using GPU_PARALLEL mode for score calculation");
                totalScore = calculateScoreGPU();
                break;
            case CPU_PARALLEL:
                System.out.println("Using CPU_PARALLEL mode for score calculation");
                totalScore = calculateScoreCPUParallel();
                break;
            case CPU_SINGLE:
            default:
                System.out.println("Using CPU_SINGLE mode for score calculation");
                totalScore = calculateScoreCPUSingle();
                break;
        }
        
        long endTime = System.currentTimeMillis();
        
        System.out.println("Time: " + (endTime - startTime) + " ms");
        System.out.println("===============================================\n");
        
        return totalScore;
    }
    
    /**
     * CPU single-threaded score calculation.
     */
    private double calculateScoreCPUSingle() {
        double totalScore = 0.0;
        
        for (RangeBipartition speciesBip : speciesBipartitions) {
            double bipScore = calculateBipartitionWeight(speciesBip);
            totalScore += bipScore;
        }
        
        return totalScore;
    }
    
    /**
     * CPU parallel score calculation.
     */
    private double calculateScoreCPUParallel() {
        int numBips = speciesBipartitions.size();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        Threading.startThreading(numThreads);
        
        int chunkSize = Math.max(1, (numBips + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (numBips + chunkSize - 1) / chunkSize);
        
        CountDownLatch latch = new CountDownLatch(actualThreads);
        AtomicLong[] partialScores = new AtomicLong[actualThreads];
        
        for (int i = 0; i < actualThreads; i++) {
            partialScores[i] = new AtomicLong(0);
            final int threadId = i;
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, numBips);
            
            Threading.execute(() -> {
                try {
                    double localScore = 0.0;
                    for (int j = startIdx; j < endIdx; j++) {
                        localScore += calculateBipartitionWeight(speciesBipartitions.get(j));
                    }
                    partialScores[threadId].set(Double.doubleToLongBits(localScore));
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Score calculation interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        double totalScore = 0.0;
        for (int i = 0; i < actualThreads; i++) {
            totalScore += Double.longBitsToDouble(partialScores[i].get());
        }
        
        return totalScore;
    }
    
    /**
     * GPU parallel score calculation.
     * Reuses the existing launchCompactWeightCalculation kernel.
     */
    private double calculateScoreGPU() {
        System.out.println("Preparing data for GPU kernel...");
        
        int numSpeciesBips = speciesBipartitions.size();
        int numGeneTreeBips = geneTrees.rangeBipartitions.size();
        int totalTrees = numGeneTrees + 1;  // Gene trees + species tree
        
        System.out.println("GPU Parameters:");
        System.out.println("  Species bipartitions: " + numSpeciesBips);
        System.out.println("  Gene tree bipartitions: " + numGeneTreeBips);
        System.out.println("  Total trees (gene + species): " + totalTrees);
        System.out.println("  Taxa: " + numTaxa);
        
        try {
            // CRITICAL: Use Structure.toArray() for contiguous memory allocation
            // This is required by JNA for passing structure arrays to native code
            
            // Prepare species tree bipartitions using contiguous memory
            CompactBipartition speciesTemplate = new CompactBipartition();
            CompactBipartition[] speciesBipsArray = 
                (CompactBipartition[]) speciesTemplate.toArray(numSpeciesBips);
            
            for (int i = 0; i < numSpeciesBips; i++) {
                RangeBipartition rb = speciesBipartitions.get(i);
                speciesBipsArray[i].geneTreeIndex = rb.geneTreeIndex;  // = numGeneTrees for species tree
                speciesBipsArray[i].leftStart = rb.leftStart;
                speciesBipsArray[i].leftEnd = rb.leftEnd;
                speciesBipsArray[i].rightStart = rb.rightStart;
                speciesBipsArray[i].rightEnd = rb.rightEnd;
            }
            System.out.println("Species bipartitions array prepared (contiguous memory)");
            
            // Prepare gene tree bipartitions using contiguous memory
            CompactBipartition geneTreeTemplate = new CompactBipartition();
            CompactBipartition[] geneTreeBipsArray = 
                (CompactBipartition[]) geneTreeTemplate.toArray(numGeneTreeBips);
            int[] frequencies = new int[numGeneTreeBips];
            
            int idx = 0;
            for (Map.Entry<RangeBipartition, Integer> entry : geneTrees.rangeBipartitions.entrySet()) {
                RangeBipartition rb = entry.getKey();
                geneTreeBipsArray[idx].geneTreeIndex = rb.geneTreeIndex;
                geneTreeBipsArray[idx].leftStart = rb.leftStart;
                geneTreeBipsArray[idx].leftEnd = rb.leftEnd;
                geneTreeBipsArray[idx].rightStart = rb.rightStart;
                geneTreeBipsArray[idx].rightEnd = rb.rightEnd;
                frequencies[idx] = entry.getValue();
                idx++;
            }
            System.out.println("Gene tree bipartitions array prepared (contiguous memory)");
            
            // Flatten the occurrence-aware STELAR-Pro index for CUDA. The old
            // dense inverse index stored one position per taxon and therefore
            // could not represent paralogs correctly.
            GpuRangeIndexData gpuIndexData = buildGpuRangeIndexData();
            
            // Output array for weights
            double[] weights = new double[numSpeciesBips];
            
            System.out.println("Launching GPU kernel...");
            
            // Call the STELAR-Pro compact weight calculation kernel.
            try {
                WeightCalcLib.INSTANCE.launchCompactWeightCalculationPro(
                    speciesBipsArray,
                    geneTreeBipsArray,
                    frequencies,
                    weights,
                    gpuIndexData.positionOffsetsMemory,
                    gpuIndexData.positionsMemory,
                    gpuIndexData.orderingMemory,
                    gpuIndexData.orderingOffsetsMemory,
                    numSpeciesBips,
                    numGeneTreeBips,
                    totalTrees,
                    numTaxa
                );
            } catch (UnsatisfiedLinkError e) {
                System.err.println("STELAR-Pro GPU scorer kernel is not available: " + e.getMessage());
                System.err.println("Rebuild cuda/libweight_calc.so with nvcc to enable GPU scoring.");
                System.err.println("Falling back to CPU_PARALLEL mode...");
                return calculateScoreCPUParallel();
            }
            
            // Sum up weights to get total score
            double totalScore = 0.0;
            for (double w : weights) {
                totalScore += w;
            }
            
            System.out.println("GPU kernel completed successfully");
            return totalScore;
            
        } catch (Exception e) {
            System.err.println("GPU kernel failed: " + e.getMessage());
            e.printStackTrace();
            System.err.println("Falling back to CPU_PARALLEL mode...");
            return calculateScoreCPUParallel();
        }
    }
    
    /**
     * Build extended inverse index that includes the species tree.
     * Species tree gets index = numGeneTrees.
     */
    private void buildExtendedInverseIndex(Tree speciesTree) {
        int totalTrees = numGeneTrees + 1;
        extendedInverseIndex = new int[totalTrees][numTaxa];
        extendedTaxonPositions = new int[totalTrees][numTaxa][];
        extendedOrderings = new int[totalTrees][];
        
        // Copy gene tree inverse indices from existing InverseIndexManager
        // We'll build our own for gene trees too for simplicity
        for (int t = 0; t < numGeneTrees; t++) {
            Tree geneTree = geneTrees.geneTrees.get(t);
            List<Integer> ordering = new ArrayList<>();
            collectLeavesInOrder(geneTree.root, ordering);
            
            extendedOrderings[t] = ordering.stream().mapToInt(Integer::intValue).toArray();
            
            buildTreePositionIndex(t);
        }
        
        // Add species tree at index numGeneTrees
        int speciesIdx = numGeneTrees;
        List<Integer> speciesOrdering = new ArrayList<>();
        collectLeavesInOrder(speciesTree.root, speciesOrdering);
        
        extendedOrderings[speciesIdx] = speciesOrdering.stream().mapToInt(Integer::intValue).toArray();
        buildTreePositionIndex(speciesIdx);
        
        System.out.println("Extended inverse index built for " + totalTrees + " trees");
    }

    /**
     * Build both the legacy first-position view and the STELAR-Pro occurrence
     * vectors for one tree in the extended scorer index.
     */
    private void buildTreePositionIndex(int treeIdx) {
        Arrays.fill(extendedInverseIndex[treeIdx], -1);

        @SuppressWarnings("unchecked")
        List<Integer>[] positionsByTaxon = new ArrayList[numTaxa];
        for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
            positionsByTaxon[taxonId] = new ArrayList<>();
        }

        for (int pos = 0; pos < extendedOrderings[treeIdx].length; pos++) {
            int taxonId = extendedOrderings[treeIdx][pos];
            if (taxonId >= 0 && taxonId < numTaxa) {
                positionsByTaxon[taxonId].add(pos);
                if (extendedInverseIndex[treeIdx][taxonId] == -1) {
                    extendedInverseIndex[treeIdx][taxonId] = pos;
                }
            }
        }

        for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
            List<Integer> positions = positionsByTaxon[taxonId];
            int[] compact = new int[positions.size()];
            for (int i = 0; i < positions.size(); i++) {
                compact[i] = positions.get(i);
            }
            extendedTaxonPositions[treeIdx][taxonId] = compact;
        }
    }
    
    /**
     * Collect leaves in left-to-right order (inorder traversal).
     */
    private void collectLeavesInOrder(TreeNode node, List<Integer> ordering) {
        if (node == null) return;
        
        if (node.isLeaf()) {
            if (node.taxon != null) {
                ordering.add(node.taxon.id);
            }
            return;
        }
        
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                collectLeavesInOrder(child, ordering);
            }
        }
    }
    
    /**
     * Extract bipartitions from species tree as RangeBipartitions.
     * Species tree has treeIndex = numGeneTrees.
     */
    private void extractSpeciesBipartitions(Tree speciesTree) {
        speciesBipartitions = new ArrayList<>();
        if (speciesTree.root != null) {
            extractBipartitionsRecursive(speciesTree.root, numGeneTrees);
        }
    }
    
    /**
     * Recursively extract bipartitions.
     * Returns [minPos, maxPos] in the species tree ordering.
     */
    private int[] extractBipartitionsRecursive(TreeNode node, int treeIndex) {
        if (node.isLeaf()) {
            if (node.taxon != null) {
                int pos = extendedInverseIndex[treeIndex][node.taxon.id];
                return new int[]{pos, pos};
            }
            return null;
        }
        
        if (node.childs == null || node.childs.size() < 2) {
            // Not a valid binary internal node
            int minPos = Integer.MAX_VALUE;
            int maxPos = Integer.MIN_VALUE;
            if (node.childs != null) {
                for (TreeNode child : node.childs) {
                    int[] range = extractBipartitionsRecursive(child, treeIndex);
                    if (range != null) {
                        minPos = Math.min(minPos, range[0]);
                        maxPos = Math.max(maxPos, range[1]);
                    }
                }
            }
            return (minPos == Integer.MAX_VALUE) ? null : new int[]{minPos, maxPos};
        }
        
        // Binary node - get ranges for left and right subtrees
        int[] leftRange = extractBipartitionsRecursive(node.childs.get(0), treeIndex);
        int[] rightRange = extractBipartitionsRecursive(node.childs.get(1), treeIndex);
        
        // Handle additional children (multifurcating)
        for (int i = 2; i < node.childs.size(); i++) {
            int[] extraRange = extractBipartitionsRecursive(node.childs.get(i), treeIndex);
            if (extraRange != null && rightRange != null) {
                rightRange[0] = Math.min(rightRange[0], extraRange[0]);
                rightRange[1] = Math.max(rightRange[1], extraRange[1]);
            } else if (extraRange != null) {
                rightRange = extraRange;
            }
        }
        
        if (leftRange != null && rightRange != null) {
            // Create RangeBipartition for this internal node
            // Note: ranges are [inclusive, inclusive], convert to [inclusive, exclusive)
            RangeBipartition bip = new RangeBipartition(
                treeIndex,
                leftRange[0], leftRange[1] + 1,   // left range
                rightRange[0], rightRange[1] + 1  // right range
            );
            speciesBipartitions.add(bip);
            
            // Return combined range
            return new int[]{
                Math.min(leftRange[0], rightRange[0]),
                Math.max(leftRange[1], rightRange[1])
            };
        }
        
        return (leftRange != null) ? leftRange : rightRange;
    }
    
    /**
     * Calculate weight for a species tree bipartition.
     * Sum over all gene tree bipartitions.
     */
    private double calculateBipartitionWeight(RangeBipartition speciesBip) {
        double totalScore = 0.0;
        
        for (Map.Entry<RangeBipartition, Integer> entry : geneTrees.rangeBipartitions.entrySet()) {
            RangeBipartition geneBip = entry.getKey();
            int frequency = entry.getValue();
            
            double score = calculateRangeScore(speciesBip, geneBip);
            totalScore += score * frequency;
        }
        
        return totalScore;
    }
    
    /**
     * Calculate score between two RangeBipartitions using inverse index.
     * Same formula as MemoryOptimizedWeightCalculator.
     */
    private double calculateRangeScore(RangeBipartition bip1, RangeBipartition bip2) {
        // Calculate four intersection sizes
        int aa = getRangeIntersectionSize(
            bip1.geneTreeIndex, bip1.leftStart, bip1.leftEnd,
            bip2.geneTreeIndex, bip2.leftStart, bip2.leftEnd);
        
        int bb = getRangeIntersectionSize(
            bip1.geneTreeIndex, bip1.rightStart, bip1.rightEnd,
            bip2.geneTreeIndex, bip2.rightStart, bip2.rightEnd);
        
        int ab = getRangeIntersectionSize(
            bip1.geneTreeIndex, bip1.leftStart, bip1.leftEnd,
            bip2.geneTreeIndex, bip2.rightStart, bip2.rightEnd);
        
        int ba = getRangeIntersectionSize(
            bip1.geneTreeIndex, bip1.rightStart, bip1.rightEnd,
            bip2.geneTreeIndex, bip2.leftStart, bip2.leftEnd);
        
        // Triplet scoring formula
        double score1 = 0;
        if (aa + bb >= 2) {
            score1 = aa * bb * (aa + bb - 2) / 2.0;
        }
        
        double score2 = 0;
        if (ab + ba >= 2) {
            score2 = ab * ba * (ab + ba - 2) / 2.0;
        }
        
        return score1 + score2;
    }
    
    /**
     * Calculate intersection size between two ranges using inverse index.
     * Complexity: O(min(|range1|, |range2|))
     */
    private int getRangeIntersectionSize(int tree1, int start1, int end1,
                                         int tree2, int start2, int end2) {
        if (start1 < 0 || end1 <= start1 || start2 < 0 || end2 <= start2) {
            return 0;
        }
        
        int size1 = end1 - start1;
        int size2 = end2 - start2;
        
        // Iterate over smaller range
        if (size1 <= size2) {
            return countIntersection(tree1, start1, end1, tree2, start2, end2);
        } else {
            return countIntersection(tree2, start2, end2, tree1, start1, end1);
        }
    }
    
    /**
     * Count intersection by iterating over smaller range.
     */
    private int countIntersection(int smallTree, int smallStart, int smallEnd,
                                  int largeTree, int largeStart, int largeEnd) {
        int count = 0;
        int[] smallOrdering = extendedOrderings[smallTree];
        
        if (smallOrdering == null) return 0;
        
        int maxPos = Math.min(smallEnd, smallOrdering.length);
        Set<Integer> seenInSmallRange = new HashSet<>();
        
        for (int pos = smallStart; pos < maxPos; pos++) {
            int taxonId = smallOrdering[pos];
            
            if (taxonId < 0 || taxonId >= numTaxa) continue;

            // STELAR-Pro scores species-set intersections. Multiple copies of
            // the same taxon in the scanned range count once, and membership in
            // the other range is tested against all occurrence positions.
            if (!seenInSmallRange.add(taxonId)) {
                continue;
            }

            if (hasPositionInRange(largeTree, taxonId, largeStart, largeEnd)) {
                count++;
            }
        }
        
        return count;
    }

    private boolean hasPositionInRange(int treeIdx, int taxonId, int start, int end) {
        int[] positions = extendedTaxonPositions[treeIdx][taxonId];
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
     * JNA memory bundle for the occurrence-aware GPU scorer index.
     */
    private static class GpuRangeIndexData {
        Memory positionOffsetsMemory;
        Memory positionsMemory;
        Memory orderingMemory;
        Memory orderingOffsetsMemory;
    }

    /**
     * Flatten extendedTaxonPositions and extendedOrderings into the same Pro GPU
     * data format used by MemoryOptimizedWeightCalculator.
     */
    @SuppressWarnings("resource")
    private GpuRangeIndexData buildGpuRangeIndexData() {
        int totalTrees = numGeneTrees + 1;
        int taxonEntryCount = totalTrees * numTaxa;

        int totalPositionEntries = 0;
        for (int tree = 0; tree < totalTrees; tree++) {
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                totalPositionEntries += extendedTaxonPositions[tree][taxon].length;
            }
        }

        int totalOrderingEntries = 0;
        for (int tree = 0; tree < totalTrees; tree++) {
            totalOrderingEntries += extendedOrderings[tree].length;
        }

        GpuRangeIndexData data = new GpuRangeIndexData();
        data.positionOffsetsMemory = new Memory((long) (taxonEntryCount + 1) * Integer.BYTES);
        data.positionsMemory = new Memory((long) Math.max(1, totalPositionEntries) * Integer.BYTES);
        data.orderingOffsetsMemory = new Memory((long) (totalTrees + 1) * Integer.BYTES);
        data.orderingMemory = new Memory((long) Math.max(1, totalOrderingEntries) * Integer.BYTES);

        int positionCursor = 0;
        for (int tree = 0; tree < totalTrees; tree++) {
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                int entryIndex = tree * numTaxa + taxon;
                data.positionOffsetsMemory.setInt((long) entryIndex * Integer.BYTES, positionCursor);
                for (int position : extendedTaxonPositions[tree][taxon]) {
                    data.positionsMemory.setInt((long) positionCursor * Integer.BYTES, position);
                    positionCursor++;
                }
            }
        }
        data.positionOffsetsMemory.setInt((long) taxonEntryCount * Integer.BYTES, positionCursor);

        int orderingCursor = 0;
        for (int tree = 0; tree < totalTrees; tree++) {
            data.orderingOffsetsMemory.setInt((long) tree * Integer.BYTES, orderingCursor);
            for (int taxonId : extendedOrderings[tree]) {
                data.orderingMemory.setInt((long) orderingCursor * Integer.BYTES, taxonId);
                orderingCursor++;
            }
        }
        data.orderingOffsetsMemory.setInt((long) totalTrees * Integer.BYTES, orderingCursor);

        System.out.println("STELAR-Pro GPU scorer index flattened:");
        System.out.println("  Taxon occurrence positions: " + totalPositionEntries);
        System.out.println("  Leaf ordering entries: " + totalOrderingEntries);

        return data;
    }
}
