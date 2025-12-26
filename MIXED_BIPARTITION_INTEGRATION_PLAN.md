# Mixed Bipartition Integration Plan

## Overview

This document outlines the changes needed to integrate `MixedBipartition` (cross-tree bipartitions) into the weight computation and DP inference pipeline.

---

## Current State

### Current Representation: `RangeBipartition`
```java
RangeBipartition {
    int geneTreeIndex;    // Single tree for both sides
    int leftStart, leftEnd;
    int rightStart, rightEnd;
}
```
- Both sides (A and B) come from the **same gene tree**
- Ranges are contiguous in the gene tree's left-to-right ordering
- Intersection computation uses `InverseIndexManager.getRangeIntersectionSize()`

### New Representation: `MixedBipartition`
```java
MixedBipartition {
    int leftTreeIndex;     // Gene tree for left side
    int leftStart, leftEnd;
    int rightTreeIndex;    // Gene tree for right side (can be different!)
    int rightStart, rightEnd;
}
```
- Sides A and B can come from **different gene trees**
- Created via cross-tree recombination in `CandidateExtender`

---

## Key Insight: Intersection Computation is Already Cross-Tree!

### Current Method Already Handles Different Trees!
The `InverseIndexManager.getRangeIntersectionSize()` already supports ranges from **different trees**:

```java
// This ALREADY works for ranges from different trees!
int intersection = inverseIndexManager.getRangeIntersectionSize(
    tree1, start1, end1,  // Range in tree1
    tree2, start2, end2   // Range in tree2 (can be different from tree1!)
);
```

This means **no special set-based intersection is needed** for `MixedBipartition`!

### For MixedBipartition: Same Logic, Different Accessors

The only change is accessing the tree index separately for each side:

```java
// RangeBipartition (both sides same tree):
calculateRangeScore(candidate, geneTree) {
    int aa = getRangeIntersectionSize(
        candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd,   // same tree
        geneTree.geneTreeIndex, geneTree.leftStart, geneTree.leftEnd
    );
    int bb = getRangeIntersectionSize(
        candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd, // same tree
        geneTree.geneTreeIndex, geneTree.rightStart, geneTree.rightEnd
    );
    // ...
}

// MixedBipartition (sides may be from different trees):
calculateMixedScore(mixed, geneTree) {
    int aa = getRangeIntersectionSize(
        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,    // left's tree
        geneTree.geneTreeIndex, geneTree.leftStart, geneTree.leftEnd
    );
    int bb = getRangeIntersectionSize(
        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd, // right's tree (may differ!)
        geneTree.geneTreeIndex, geneTree.rightStart, geneTree.rightEnd
    );
    // ...
}
```

**This is a trivial change!** The core intersection logic stays the same.

---

## Changes Required

### 1. No Interface Needed!

Since the intersection computation is already cross-tree capable, we don't need a common interface.
Both `RangeBipartition` and `MixedBipartition` can be handled with their existing structures.

---

### 2. No Changes to InverseIndexManager Needed!

The existing `getRangeIntersectionSize(tree1, start1, end1, tree2, start2, end2)` already handles
ranges from different trees. No modifications required!

---

### 3. Add MixedBipartition Score Method to MemoryOptimizedWeightCalculator

Since the intersection logic is the same, we just add a new method to the existing calculator:

```java
// In MemoryOptimizedWeightCalculator.java - ADD this method:

/**
 * Calculate score between a MixedBipartition and a RangeBipartition.
 * Uses the SAME intersection logic - just accesses tree indices separately for each side.
 */
private double calculateMixedScore(MixedBipartition mixed, RangeBipartition geneTree) {
    // Same intersection logic - just use leftTreeIndex/rightTreeIndex separately
    int aa = inverseIndexManager.getRangeIntersectionSize(
        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
        geneTree.geneTreeIndex, geneTree.leftStart, geneTree.leftEnd);
        
    int bb = inverseIndexManager.getRangeIntersectionSize(
        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
        geneTree.geneTreeIndex, geneTree.rightStart, geneTree.rightEnd);
        
    int ab = inverseIndexManager.getRangeIntersectionSize(
        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
        geneTree.geneTreeIndex, geneTree.rightStart, geneTree.rightEnd);
        
    int ba = inverseIndexManager.getRangeIntersectionSize(
        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
        geneTree.geneTreeIndex, geneTree.leftStart, geneTree.leftEnd);
    
    // Same scoring formula
    double score1 = (aa + bb >= 2) ? aa * bb * (aa + bb - 2) / 2.0 : 0;
    double score2 = (ab + ba >= 2) ? ab * ba * (ab + ba - 2) / 2.0 : 0;
    
    return score1 + score2;
}

/**
 * Calculate total weight for a MixedBipartition.
 */
public double calculateMixedWeight(MixedBipartition mixed) {
    double totalScore = 0.0;
    
    for (Map.Entry<RangeBipartition, Integer> entry : geneTrees.rangeBipartitions.entrySet()) {
        double score = calculateMixedScore(mixed, entry.getKey());
        totalScore += score * entry.getValue();
    }
    
    return totalScore;
}
```

**NO new file needed!** Just add methods to existing `MemoryOptimizedWeightCalculator.java`

---

### 4. Update MemoryOptimizedInferenceDP

The DP needs to handle both `RangeBipartition` and `MixedBipartition` candidates:

**Option A: Unified Candidate List**
- Convert everything to a common representation
- Store weights in `Map<Object, Double>` where key can be either type

**Option B: Separate Handling**
- Keep `RangeBipartition` candidates with existing logic
- Add separate `MixedBipartition` candidates with new weight calculation
- Merge in the DP lookup

**Recommended: Option B** (less invasive, backward compatible)

```java
// Existing
private Map<ClusterHashPair, List<RangeBipartition>> clusterHashToRangeBips;
private Map<RangeBipartition, Double> rangeBipWeights;

// NEW: Add for mixed bipartitions
private Map<ClusterHashPair, List<MixedBipartition>> clusterHashToMixedBips;
private Map<MixedBipartition, Double> mixedBipWeights;
```

**File to modify:** `src/core/MemoryOptimizedInferenceDP.java`

---

### 5. Update DP Function

The `dp()` function needs to consider both types of candidates:

```java
private double dp(ClusterHashPair clusterHash) {
    // ... existing memoization check ...
    
    double maxScore = Double.NEGATIVE_INFINITY;
    Object bestChoice = null;  // Can be RangeBipartition or MixedBipartition
    
    // Check RangeBipartition candidates (existing)
    List<RangeBipartition> rangeCandidates = clusterHashToRangeBips.get(clusterHash);
    if (rangeCandidates != null) {
        for (RangeBipartition rangeBip : rangeCandidates) {
            // ... existing logic ...
        }
    }
    
    // NEW: Check MixedBipartition candidates
    List<MixedBipartition> mixedCandidates = clusterHashToMixedBips.get(clusterHash);
    if (mixedCandidates != null) {
        for (MixedBipartition mixedBip : mixedCandidates) {
            // Get child cluster hashes
            ClusterHashPair leftHash = getClusterHash(mixedBip.leftTreeIndex, 
                                                       mixedBip.leftStart, mixedBip.leftEnd);
            ClusterHashPair rightHash = getClusterHash(mixedBip.rightTreeIndex,
                                                        mixedBip.rightStart, mixedBip.rightEnd);
            
            double leftScore = dp(leftHash);
            double rightScore = dp(rightHash);
            double bipScore = mixedBipWeights.getOrDefault(mixedBip, 0.0);
            
            double totalScore = leftScore + rightScore + bipScore;
            
            if (totalScore > maxScore) {
                maxScore = totalScore;
                bestChoice = mixedBip;
            }
        }
    }
    
    // ... memoize and return ...
}
```

---

### 6. Update Tree Reconstruction

`reconstructTree()` needs to handle `MixedBipartition` choices:

```java
private TreeNode buildTreeNode(ClusterHashPair clusterHash, Tree tree) {
    Object choice = dpChoice.get(clusterHash);
    
    if (choice instanceof RangeBipartition) {
        RangeBipartition rangeBip = (RangeBipartition) choice;
        // ... existing logic ...
    } else if (choice instanceof MixedBipartition) {
        MixedBipartition mixedBip = (MixedBipartition) choice;
        ClusterHashPair leftHash = getClusterHash(mixedBip.leftTreeIndex,
                                                   mixedBip.leftStart, mixedBip.leftEnd);
        ClusterHashPair rightHash = getClusterHash(mixedBip.rightTreeIndex,
                                                    mixedBip.rightStart, mixedBip.rightEnd);
        
        TreeNode leftChild = buildTreeNode(leftHash, tree);
        TreeNode rightChild = buildTreeNode(rightHash, tree);
        // ... create internal node ...
    }
}
```

**File to modify:** `src/core/MemoryOptimizedInferenceDP.java`

---

### 7. Add Command-Line Flag

Add a flag to enable/disable mixed bipartition extension:

```java
// In Main.java
boolean useMixedBipartitions = false;

// Parse --use-mixed or --extend-candidates flag
if (args[i].equals("--use-mixed") || args[i].equals("--extend-candidates")) {
    useMixedBipartitions = true;
}

// Pass to inference
if (useMixedBipartitions) {
    List<MixedBipartition> mixedBips = geneTrees.getMixedBipartitions();
    // ... include in DP ...
}
```

**File to modify:** `src/Main.java`

---

### 8. CUDA Kernel Changes (GPU Mode)

The GPU kernel needs a small modification to handle `MixedBipartition`.

#### Current CUDA Struct (RangeBipartition)
```c
// In cuda/weight_calc.cu
struct CompactBipartition {
    int geneTreeIndex;   // Single tree index for both sides
    int leftStart;
    int leftEnd;
    int rightStart;
    int rightEnd;
};
```

#### New CUDA Struct (MixedBipartition)
```c
// Add new struct for mixed bipartitions
struct MixedCompactBipartition {
    int leftTreeIndex;   // Tree index for LEFT side
    int leftStart;
    int leftEnd;
    int rightTreeIndex;  // Tree index for RIGHT side (may differ!)
    int rightStart;
    int rightEnd;
};
```

#### Kernel Change (Minimal)
The intersection computation logic stays the same - just use different tree indices:

```c
// BEFORE (RangeBipartition - same tree for both sides):
__device__ double computeScore(CompactBipartition bip, CompactBipartition gt, ...) {
    int aa = rangeIntersection(bip.geneTreeIndex, bip.leftStart, bip.leftEnd,
                               gt.geneTreeIndex, gt.leftStart, gt.leftEnd, inverseIndex);
    int bb = rangeIntersection(bip.geneTreeIndex, bip.rightStart, bip.rightEnd,
                               gt.geneTreeIndex, gt.rightStart, gt.rightEnd, inverseIndex);
    // ...
}

// AFTER (MixedBipartition - separate tree indices):
__device__ double computeMixedScore(MixedCompactBipartition bip, CompactBipartition gt, ...) {
    int aa = rangeIntersection(bip.leftTreeIndex, bip.leftStart, bip.leftEnd,    // left's tree
                               gt.geneTreeIndex, gt.leftStart, gt.leftEnd, inverseIndex);
    int bb = rangeIntersection(bip.rightTreeIndex, bip.rightStart, bip.rightEnd, // right's tree
                               gt.geneTreeIndex, gt.rightStart, gt.rightEnd, inverseIndex);
    // Same scoring formula...
}
```

#### JNA Interface Update
```java
// In WeightCalcLib.java - add new struct
public static class MixedCompactBipartition extends Structure {
    public int leftTreeIndex;
    public int leftStart;
    public int leftEnd;
    public int rightTreeIndex;
    public int rightStart;
    public int rightEnd;
    
    // ... Structure boilerplate ...
}

// Add new kernel launch method
void launchMixedWeightCalculation(
    MixedCompactBipartition[] candidates,
    CompactBipartition[] geneTreeBips,
    int[] frequencies,
    double[] weights,
    Memory inverseIndex,
    Memory orderings,
    int numCandidates,
    int numGeneTreeBips,
    int numTrees,
    int numTaxa
);
```

**Files to modify:**
- `cuda/weight_calc.cu` - Add new struct and kernel
- `src/core/WeightCalculator.java` - Add JNA interface for new struct/kernel
- `src/core/MemoryOptimizedWeightCalculator.java` - Add GPU path for mixed bipartitions

**The core intersection logic in CUDA remains unchanged!** Only the struct and which tree index to use changes.

---

## Summary of Files to Modify

| File | Changes |
|------|---------|
| `src/core/MemoryOptimizedWeightCalculator.java` | Add `calculateMixedScore()`, `calculateMixedWeight()`, GPU path |
| `src/core/MemoryOptimizedInferenceDP.java` | Add mixed bipartition handling in DP |
| `src/core/InferenceDP.java` | Update to pass mixed bipartitions |
| `src/core/WeightCalculator.java` | Add JNA interface for `MixedCompactBipartition` |
| `src/Main.java` | Add `--use-mixed` flag |
| `cuda/weight_calc.cu` | Add `MixedCompactBipartition` struct and kernel |

**Key insight:** The intersection logic is unchanged! We just use separate tree indices for left/right sides.

---

## Backward Compatibility

- Default behavior: Use only `RangeBipartition` (gene tree bipartitions)
- With `--use-mixed` flag: Also include `MixedBipartition` candidates
- All existing code paths remain functional

---

## Performance Considerations

1. **Same intersection complexity** - Reuses existing `InverseIndexManager` with O(min(|A|, |B|)) complexity
2. **Memory**: `MixedBipartition` stores 6 ints vs 5 ints for `RangeBipartition` - negligible overhead
3. **GPU**: Minimal kernel change - same intersection logic, just different tree index accessors
4. **Candidate explosion**: Cross-tree recombination may generate many candidates
   - Already mitigated by hash-based deduplication in `CandidateExtender`

---

## Testing Plan

1. Run without `--use-mixed` → should produce same results as before
2. Run with `--use-mixed` on small dataset → verify mixed bipartitions are generated
3. Compare scores: with vs without mixed bipartitions
4. Verify tree reconstruction works with mixed bipartition choices

