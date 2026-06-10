package tree;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import utils.HashUtils;
import utils.Threading;

/**
 * Memory-efficient bipartition manager that reduces memory usage from O(n²k) to O(nk).
 * 
 * Instead of immediately creating BitSet representations for all bipartitions,
 * this manager:
 * 1. Represents bipartitions as ranges during initial processing
 * 2. Uses hash-based equality checking to filter unique bipartitions
 * 3. Only converts unique bipartitions to BitSets at the end
 * 
 * This optimization is particularly effective when there are many duplicate
 * bipartitions across gene trees.
 */
public class MemoryEfficientBipartitionManager {
    
    // Configuration
    public static boolean ENABLE_DOUBLE_HASHING = true;
    public static boolean ENABLE_EXPENSIVE_EQUALITY_CHECKS = false; // Default: trust the hash functions
    public static RangeBipartition.HashFunction DEFAULT_HASH_FUNCTION = getDefaultHashFunction();
    
    private static RangeBipartition.HashFunction getDefaultHashFunction() {
        return ENABLE_DOUBLE_HASHING ? new RangeBipartition.DoubleHashFunction() : new RangeBipartition.SumHashFunction();
    }
    
    // Gene tree data structures
    private final List<Tree> geneTrees;
    private final int[][] geneTreeTaxaOrdering;  // [tree_index][taxa_position] = taxon_id
    private final long[][] prefixSums;           // Legacy occurrence-prefix sums used by older extension code
    private final long[][] prefixXORs;           // Legacy occurrence-prefix XORs used by older extension code
    
    // Processing results
    private final Map<Object, List<RangeBipartition>> hashToBipartitions;
    private final Map<RangeBipartition, Integer> uniqueRangeBipartitions;
    private final Map<RangeKey, RangeClusterInfo> rangeClusterInfoByRange;
    private final int realTaxaCount;

    /**
     * Immutable key for a contiguous occurrence range in one gene tree.
     *
     * The range is still expressed in leaf-occurrence coordinates because that
     * is how the original memory-efficient DP refers to subtrees. The value
     * associated with this key is duplicate-invariant: repeated species inside
     * the range have already been collapsed.
     */
    public static final class RangeKey {
        public final int geneTreeIndex;
        public final int start;
        public final int end;

        public RangeKey(int geneTreeIndex, int start, int end) {
            this.geneTreeIndex = geneTreeIndex;
            this.start = start;
            this.end = end;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof RangeKey)) return false;
            RangeKey other = (RangeKey) obj;
            return geneTreeIndex == other.geneTreeIndex && start == other.start && end == other.end;
        }

        @Override
        public int hashCode() {
            int result = geneTreeIndex;
            result = 31 * result + start;
            result = 31 * result + end;
            return result;
        }
    }

    /**
     * Compact, persistent cluster record for one subtree range.
     *
     * No taxon set is stored here. The temporary set used to build rawSum/rawXor
     * is released during traversal, leaving only O(1) data needed by the DP:
     * duplicate-collapsed size plus the final cluster hash.
     */
    public static final class RangeClusterInfo {
        public final int geneTreeIndex;
        public final int start;
        public final int end;
        public final int uniqueTaxonCount;
        public final long rawSum;
        public final long rawXor;
        public final ClusterHashPair hash;

        public RangeClusterInfo(int geneTreeIndex, int start, int end,
                                int uniqueTaxonCount, long rawSum, long rawXor,
                                ClusterHashPair hash) {
            this.geneTreeIndex = geneTreeIndex;
            this.start = start;
            this.end = end;
            this.uniqueTaxonCount = uniqueTaxonCount;
            this.rawSum = rawSum;
            this.rawXor = rawXor;
            this.hash = hash;
        }
    }

    /**
     * Mutable state used only while one gene tree is being traversed.
     *
     * The taxa set is the small-to-large merge payload. It is intentionally not
     * retained after its subtree has been merged upward, because the permanent
     * RangeClusterInfo cache already contains the compact hash/count summary.
     */
    private static final class ClusterBuildState {
        int geneTreeIndex;
        int start;
        int end;
        int uniqueTaxonCount;
        long rawSum;
        long rawXor;
        ClusterHashPair hash;
        Set<Integer> taxa;

        void releaseTaxa() {
            if (taxa != null) {
                taxa.clear();
                taxa = null;
            }
        }
    }
    
    public MemoryEfficientBipartitionManager(List<Tree> geneTrees, int realTaxaCount) {
        this.geneTrees = geneTrees;
        this.realTaxaCount = realTaxaCount;
        this.geneTreeTaxaOrdering = new int[geneTrees.size()][];
        this.prefixSums = new long[geneTrees.size()][];
        this.prefixXORs = new long[geneTrees.size()][];
        this.hashToBipartitions = new ConcurrentHashMap<>();
        this.uniqueRangeBipartitions = new ConcurrentHashMap<>();
        this.rangeClusterInfoByRange = new ConcurrentHashMap<>();
        
        initializeGeneTreeOrderings();
        calculatePrefixArrays();
    }
    
    /**
     * Initialize the left-to-right ordering for each gene tree.
     *
     * The ordering stores leaf occurrences, not only species IDs. This distinction
     * matters for STELAR-Pro because a rooted/tagged gene family tree can contain
     * several leaves with the same species label. The integer in the ordering is
     * still the species/taxon ID, but each leaf node also receives its occurrence
     * position through TreeNode.traversalIndex. Later range construction uses that
     * occurrence position so duplicate leaves remain in their true tree locations.
     */
    private void initializeGeneTreeOrderings() {
        System.out.println("Initializing gene tree taxa orderings...");
        
        for (int i = 0; i < geneTrees.size(); i++) {
            Tree tree = geneTrees.get(i);
            List<Integer> ordering = new ArrayList<>();
            collectLeavesInOrder(tree.root, ordering);
            
            geneTreeTaxaOrdering[i] = ordering.stream().mapToInt(Integer::intValue).toArray();
        }
        
        System.out.println("Initialized orderings for " + geneTrees.size() + " gene trees");
    }
    
    /**
     * Collect leaves in left-to-right order and remember each leaf occurrence.
     */
    private void collectLeavesInOrder(TreeNode node, List<Integer> ordering) {
        if (node.isLeaf()) {
            node.traversalIndex = ordering.size();
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
     * Calculate prefix sums and XORs over hashed taxon IDs for efficient range hash computation.
     * This provides much better hash distribution than using raw taxon IDs.
     */
    private void calculatePrefixArrays() {
        System.out.println("Calculating prefix sums and XORs over hashed taxon IDs for hash functions...");
        
        for (int i = 0; i < geneTrees.size(); i++) {
            int[] ordering = geneTreeTaxaOrdering[i];
            prefixSums[i] = new long[ordering.length];
            prefixXORs[i] = new long[ordering.length];
            
            if (ordering.length > 0) {
                // Hash the first taxon ID and store as initial values
                long hashedTaxon0 = hashSingleTaxon(ordering[0]);
                prefixSums[i][0] = hashedTaxon0;
                prefixXORs[i][0] = hashedTaxon0;
                
                for (int j = 1; j < ordering.length; j++) {
                    long hashedTaxon = hashSingleTaxon(ordering[j]);
                    prefixSums[i][j] = prefixSums[i][j - 1] + hashedTaxon;
                    prefixXORs[i][j] = prefixXORs[i][j - 1] ^ hashedTaxon;
                }
            }
        }
        
        System.out.println("Calculated prefix sums and XORs over hashed taxon IDs for efficient range hashing");
    }
    
    /**
     * Simple but effective hash function for individual taxon IDs.
     * Uses a combination of multiplications and XOR operations for good distribution.
     */
    private static long hashSingleTaxon(int taxonId) {
        // Use a strong hash mixing function for individual taxon IDs
        long x = taxonId;
        x ^= x >>> 16;
        x *= 0x85ebca6b;
        x ^= x >>> 13;
        x *= 0xc2b2ae35;
        x ^= x >>> 16;
        
        // Ensure we never return 0 to avoid issues with prefix operations
        return x == 0 ? 1 : x;
    }
    
    /**
     * Process all gene trees to extract range bipartitions in parallel.
     */
    public Map<RangeBipartition, Integer> processGeneTreesParallel() {
        System.out.println("Processing gene trees with memory-efficient range bipartitions...");

        hashToBipartitions.clear();
        uniqueRangeBipartitions.clear();
        rangeClusterInfoByRange.clear();
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        // Calculate optimal number of threads to avoid invalid ranges
        int chunkSize = Math.max(1, (geneTrees.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (geneTrees.size() + chunkSize - 1) / chunkSize);
        
        CountDownLatch latch = new CountDownLatch(actualThreads);
        AtomicInteger processedTrees = new AtomicInteger(0);
        
        System.out.println("Using " + actualThreads + " threads for parallel processing");
        
        // Process gene trees in parallel chunks
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, geneTrees.size());
            final int threadId = i;
            
            // Skip threads that would have invalid ranges
            if (startIdx >= geneTrees.size()) {
                continue;
            }
            
            Threading.execute(() -> {
                try {
                    Map<Object, List<RangeBipartition>> localHashMap = new HashMap<>();
                    
                    // Validate range before processing
                    if (startIdx >= endIdx || startIdx >= geneTrees.size()) {
                        System.out.println("Thread " + threadId + " skipped - invalid range [" + startIdx + ", " + endIdx + ")");
                        return;
                    }
                    
                    System.out.println("Thread " + threadId + " processing trees " + startIdx + " to " + (endIdx - 1));
                    
                    for (int treeIdx = startIdx; treeIdx < endIdx; treeIdx++) {
                        Tree tree = geneTrees.get(treeIdx);
                        if (tree.isRooted()) {
                            ClusterBuildState rootState = extractRangeBipartitions(tree.root, treeIdx, localHashMap);
                            if (rootState != null) {
                                rootState.releaseTaxa();
                            }
                        }
                        processedTrees.incrementAndGet();
                    }
                    
                    // Merge local results into global map
                    synchronized (hashToBipartitions) {
                        for (Map.Entry<Object, List<RangeBipartition>> entry : localHashMap.entrySet()) {
                            hashToBipartitions.computeIfAbsent(entry.getKey(), k -> new ArrayList<>())
                                             .addAll(entry.getValue());
                        }
                    }
                    
                    System.out.println("Thread " + threadId + " completed processing " + (endIdx - startIdx) + " trees");
                    
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Gene tree processing was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        System.out.println("Processed " + processedTrees.get() + " gene trees");
        System.out.println("Found " + hashToBipartitions.size() + " unique hash groups");
        
        // Convert range bipartitions to frequency map
        convertToFrequencyMap();
        
        return uniqueRangeBipartitions;
    }
    
    /**
     * Extract range bipartitions and subtree hashes in one postorder traversal.
     *
     * This is the efficient STELAR-Pro path. Each child returns a temporary set
     * of unique species in its subtree plus raw duplicate-invariant sum/XOR
     * accumulators. The parent first uses the child summaries to add a candidate
     * bipartition if the parent is a speciation node, then merges the child sets
     * with the small-to-large rule to build its own summary.
     *
     * Duplication nodes are deliberately not added as candidate bipartitions, but
     * they are still traversed and summarized because speciation descendants and
     * ancestor union clusters are needed by the DP.
     */
    private ClusterBuildState extractRangeBipartitions(TreeNode node, int treeIndex,
                                                       Map<Object, List<RangeBipartition>> localHashMap) {
        if (node == null) {
            return null;
        }

        if (node.isLeaf()) {
            return createLeafClusterState(node, treeIndex);
        }

        List<ClusterBuildState> childStates = new ArrayList<>();
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                ClusterBuildState childState = extractRangeBipartitions(child, treeIndex, localHashMap);
                if (childState != null) {
                    childStates.add(childState);
                }
            }
        }

        if (childStates.isEmpty()) {
            return null;
        }

        // STELAR-Pro only uses speciation-driven triplets. Therefore the
        // candidate split rooted at this node is admitted only when the tag says
        // this is not a duplication. The children have already been collapsed to
        // unique species summaries, so duplicate copies inside either side do
        // not change the candidate identity.
        if (node.childs != null && node.childs.size() == 2 && childStates.size() == 2 && !node.isDuplicationNode) {
            ClusterBuildState left = childStates.get(0);
            ClusterBuildState right = childStates.get(1);

            if (left.start <= left.end && right.start <= right.end) {
                RangeBipartition rangeBip = new RangeBipartition(treeIndex,
                    left.start, left.end,
                    right.start, right.end);

                Object hash = calculateSpeciesBipartitionHash(
                    left.hash, left.uniqueTaxonCount,
                    right.hash, right.uniqueTaxonCount);

                localHashMap.computeIfAbsent(hash, k -> new ArrayList<>()).add(rangeBip);
            }
        }

        ClusterBuildState merged = mergeChildClusterStates(childStates);
        storeRangeClusterInfo(merged);
        return merged;
    }

    /**
     * Build the base cluster for a leaf occurrence.
     *
     * A duplicated taxon may appear as several different leaf occurrences in a
     * gene tree, but each leaf cluster itself contains exactly one unique
     * species. Higher nodes collapse repeated species while merging these sets.
     */
    private ClusterBuildState createLeafClusterState(TreeNode node, int treeIndex) {
        if (node.traversalIndex < 0) {
            return null;
        }

        int taxonId = node.taxon.id;
        long hashedTaxon = HashUtils.hashSingleTaxon(taxonId);

        ClusterBuildState state = new ClusterBuildState();
        state.geneTreeIndex = treeIndex;
        state.start = node.traversalIndex;
        state.end = node.traversalIndex + 1;
        state.uniqueTaxonCount = 1;
        state.rawSum = hashedTaxon;
        state.rawXor = hashedTaxon;
        state.hash = HashUtils.computeClusterHashFromRaw(state.rawSum, state.rawXor, state.uniqueTaxonCount);
        state.taxa = new HashSet<>(1);
        state.taxa.add(taxonId);

        storeRangeClusterInfo(state);
        return state;
    }

    /**
     * Merge child cluster states with the small-to-large technique.
     *
     * The largest child set becomes the parent set. Every smaller set is scanned
     * once and immediately cleared after its species have been inserted into the
     * large set. Across a tree, each surviving set entry migrates from a
     * smaller set to a set at least twice as large. Entries for duplicated
     * species disappear once they collide with an existing copy, so the bound is
     * O(m log m) worst-case for m leaf occurrences and often lower on highly
     * duplicated gene trees.
     */
    private ClusterBuildState mergeChildClusterStates(List<ClusterBuildState> childStates) {
        ClusterBuildState base = childStates.get(0);
        int minStart = Integer.MAX_VALUE;
        int maxEnd = Integer.MIN_VALUE;

        for (ClusterBuildState child : childStates) {
            if (child.taxa != null && (base.taxa == null || child.taxa.size() > base.taxa.size())) {
                base = child;
            }
            minStart = Math.min(minStart, child.start);
            maxEnd = Math.max(maxEnd, child.end);
        }

        long rawSum = base.rawSum;
        long rawXor = base.rawXor;
        int uniqueTaxonCount = base.uniqueTaxonCount;

        for (ClusterBuildState child : childStates) {
            if (child == base || child.taxa == null) {
                continue;
            }

            for (int taxonId : child.taxa) {
                // Hash contribution is added only on the first insertion into
                // the surviving set. Extra copies of the same species therefore
                // do not affect rawSum/rawXor or the cluster size.
                if (base.taxa.add(taxonId)) {
                    long hashedTaxon = HashUtils.hashSingleTaxon(taxonId);
                    rawSum += hashedTaxon;
                    rawXor ^= hashedTaxon;
                    uniqueTaxonCount++;
                }
            }

            // The compact RangeClusterInfo for this child has already been
            // stored, so the temporary set can be released as soon as it has
            // been merged into the parent.
            child.releaseTaxa();
        }

        base.start = minStart;
        base.end = maxEnd;
        base.rawSum = rawSum;
        base.rawXor = rawXor;
        base.uniqueTaxonCount = uniqueTaxonCount;
        base.hash = HashUtils.computeClusterHashFromRaw(rawSum, rawXor, uniqueTaxonCount);
        return base;
    }

    /**
     * Persist the compact subtree summary used later by ClusterHashManager.
     */
    private void storeRangeClusterInfo(ClusterBuildState state) {
        if (state == null) {
            return;
        }

        RangeKey key = new RangeKey(state.geneTreeIndex, state.start, state.end);
        rangeClusterInfoByRange.put(key, new RangeClusterInfo(
            state.geneTreeIndex,
            state.start,
            state.end,
            state.uniqueTaxonCount,
            state.rawSum,
            state.rawXor,
            state.hash));
    }

    /**
     * Combine already-computed duplicate-invariant side hashes into an
     * orientation-free bipartition hash.
     *
     * This replaces the previous naive implementation that rescanned the left
     * and right occurrence ranges into Set<Integer> objects every time a
     * candidate was seen.
     */
    private Object calculateSpeciesBipartitionHash(ClusterHashPair leftHash, int leftSize,
                                                   ClusterHashPair rightHash, int rightSize) {
        ClusterHashPair first = leftHash;
        ClusterHashPair second = rightHash;
        int firstSize = leftSize;
        int secondSize = rightSize;

        if (shouldSwapCanonicalSides(leftHash, leftSize, rightHash, rightSize)) {
            first = rightHash;
            second = leftHash;
            firstSize = rightSize;
            secondSize = leftSize;
        }

        long sumComponent = first.sumHash * 0x9e3779b97f4a7c15L
            ^ second.sumHash * 0xc2b2ae3d27d4eb4fL
            ^ ((long) firstSize << 32)
            ^ secondSize;
        long xorComponent = first.xorHash
            ^ Long.rotateLeft(second.xorHash, 29)
            ^ ((long) (firstSize + secondSize) << 17);

        return new RangeBipartition.HashPair(
            Long.rotateLeft(sumComponent, 27) ^ (sumComponent >>> 33),
            splitMix64(xorComponent));
    }

    private boolean shouldSwapCanonicalSides(ClusterHashPair leftHash, int leftSize,
                                             ClusterHashPair rightHash, int rightSize) {
        if (leftSize != rightSize) {
            return leftSize > rightSize;
        }
        if (leftHash.sumHash != rightHash.sumHash) {
            return Long.compareUnsigned(leftHash.sumHash, rightHash.sumHash) > 0;
        }
        return Long.compareUnsigned(leftHash.xorHash, rightHash.xorHash) > 0;
    }

    private long splitMix64(long z) {
        z += 0x9e3779b97f4a7c15L;
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }
    
    /**
     * Convert range bipartitions to frequency map for unique ones.
     * Uses expensive equality checks if enabled, otherwise trusts hash function uniqueness.
     */
    private void convertToFrequencyMap() {
        System.out.println("Converting unique range bipartitions to frequency map...");
        System.out.println("Expensive equality checks: " + (ENABLE_EXPENSIVE_EQUALITY_CHECKS ? "ENABLED" : "DISABLED (trusting hash)"));
        
        int totalRanges = 0;
        int uniqueRanges = 0;
        
        if (ENABLE_EXPENSIVE_EQUALITY_CHECKS) {
            System.out.println("\n\nPerforming expensive equality checks...\n\n");
            // Original expensive approach with full equality checking
            for (Map.Entry<Object, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
                List<RangeBipartition> ranges = entry.getValue();
                totalRanges += ranges.size();
                
                if (ranges.isEmpty()) continue;
                
                // Group ranges by actual equality (not just hash)
                Map<RangeBipartition, Integer> actuallyUniqueRanges = new HashMap<>();
                
                for (RangeBipartition range : ranges) {
                    boolean found = false;
                    
                    for (RangeBipartition existing : actuallyUniqueRanges.keySet()) {
                        if (rangesAreEqual(range, existing)) {
                            actuallyUniqueRanges.put(existing, actuallyUniqueRanges.get(existing) + 1);
                            found = true;
                            break;
                        }
                    }
                    
                    if (!found) {
                        actuallyUniqueRanges.put(range, 1);
                    }
                }
                
                // Add unique ranges to frequency map
                for (Map.Entry<RangeBipartition, Integer> rangeEntry : actuallyUniqueRanges.entrySet()) {
                    RangeBipartition range = rangeEntry.getKey();
                    if (range != null) {
                        uniqueRangeBipartitions.merge(range, rangeEntry.getValue(), Integer::sum);
                        uniqueRanges++;
                    }
                }
            }
        } else {
            // Fast approach: trust the hash function, just take the first from each hash group
            for (Map.Entry<Object, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
                List<RangeBipartition> ranges = entry.getValue();
                totalRanges += ranges.size();
                
                if (ranges.isEmpty()) continue;
                
                // Just take the first range from each hash group and sum up all occurrences
                RangeBipartition representative = ranges.get(0);
                int totalCount = ranges.size(); // All ranges in this hash group are considered identical
                
                if (representative != null) {
                    uniqueRangeBipartitions.merge(representative, totalCount, Integer::sum);
                    uniqueRanges++;
                }
            }
        }
        
        System.out.println("Memory optimization results:");
        System.out.println("  Total range bipartitions processed: " + totalRanges);
        System.out.println("  Unique RangeBipartitions created: " + uniqueRanges);
        System.out.println("  Memory reduction factor: " + ((double) totalRanges / Math.max(1, uniqueRanges)));
        System.out.println("  Final unique RangeBipartitions: " + uniqueRangeBipartitions.size());
    }
    
    /**
     * Check equality through the cached duplicate-invariant side summaries.
     *
     * The old STELAR-X equality path compared occurrence ranges and, for
     * different trees, fell back to prefix-array hashes plus Set scans. That is
     * not correct for STELAR-Pro because two ranges with different copy counts
     * can represent the same species cluster. Here both orientations are tested
     * using the side cluster hashes/counts already produced by the traversal.
     */
    private boolean rangesAreEqual(RangeBipartition range1, RangeBipartition range2) {
        RangeClusterInfo left1 = getRangeClusterInfo(range1.geneTreeIndex, range1.leftStart, range1.leftEnd);
        RangeClusterInfo right1 = getRangeClusterInfo(range1.geneTreeIndex, range1.rightStart, range1.rightEnd);
        RangeClusterInfo left2 = getRangeClusterInfo(range2.geneTreeIndex, range2.leftStart, range2.leftEnd);
        RangeClusterInfo right2 = getRangeClusterInfo(range2.geneTreeIndex, range2.rightStart, range2.rightEnd);

        if (left1 == null || right1 == null || left2 == null || right2 == null) {
            return false;
        }

        boolean directMatch = clusterInfosEqual(left1, left2) && clusterInfosEqual(right1, right2);
        boolean symmetricMatch = clusterInfosEqual(left1, right2) && clusterInfosEqual(right1, left2);

        return directMatch || symmetricMatch;
    }

    private boolean clusterInfosEqual(RangeClusterInfo first, RangeClusterInfo second) {
        return first.uniqueTaxonCount == second.uniqueTaxonCount && first.hash.equals(second.hash);
    }

    private RangeClusterInfo getRangeClusterInfo(int geneTreeIndex, int start, int end) {
        return rangeClusterInfoByRange.get(new RangeKey(geneTreeIndex, start, end));
    }
    
    
    /**
     * Get processing statistics.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Memory-Efficient Bipartition Processing Statistics:\n");
        sb.append("  Gene trees processed: ").append(geneTrees.size()).append("\n");
        sb.append("  Hash function: ").append(DEFAULT_HASH_FUNCTION.getName()).append(" (using hashed taxon IDs)\n");
        sb.append("  Equality checking mode: ").append(ENABLE_EXPENSIVE_EQUALITY_CHECKS ? "EXPENSIVE" : "FAST (trust hash)").append("\n");
        sb.append("  Unique hash groups: ").append(hashToBipartitions.size()).append("\n");
        sb.append("  Final unique RangeBipartitions: ").append(uniqueRangeBipartitions.size()).append("\n");
        
        int totalRanges = hashToBipartitions.values().stream()
                         .mapToInt(List::size)
                         .sum();
        
        sb.append("  Total range bipartitions: ").append(totalRanges).append("\n");
        sb.append("  Memory reduction factor: ").append(String.format("%.2f", 
            (double) totalRanges / Math.max(1, uniqueRangeBipartitions.size()))).append("x\n");
        
        return sb.toString();
    }
    
    // Getters for accessing internal data structures (needed for memory-optimized weight calculation)
    
    /**
     * Get the hash-to-bipartitions mapping for memory-optimized processing.
     * This exposes the range bipartition data for efficient weight calculation.
     */
    public Map<Object, List<RangeBipartition>> getHashToBipartitions() {
        return hashToBipartitions;
    }

    /**
     * Get compact duplicate-invariant cluster summaries for subtree ranges.
     *
     * ClusterHashManager uses this map to reuse the hashes produced during
     * small-to-large traversal instead of rescanning occurrence ranges.
     */
    public Map<RangeKey, RangeClusterInfo> getRangeClusterInfoByRange() {
        return rangeClusterInfoByRange;
    }
    
    /**
     * Get the gene tree taxa orderings for inverse index construction.
     */
    public int[][] getGeneTreeTaxaOrdering() {
        return geneTreeTaxaOrdering;
    }
    
    /**
     * Get the prefix sums arrays for hash computation.
     */
    public long[][] getPrefixSums() {
        return prefixSums;
    }
    
    /**
     * Get the prefix XORs arrays for hash computation.
     */
    public long[][] getPrefixXORs() {
        return prefixXORs;
    }
}
