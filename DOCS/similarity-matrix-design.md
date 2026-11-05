# Similarity Matrix — GPU-Parallel Design

## 1. Goal

Compute a taxon-similarity matrix `M[a,b]` from k gene trees over n global taxa:

```
M[a,b] = Σ_{t: a,b ∈ T_t}  num_t(a,b)
          ─────────────────────────────
          Σ_{t: a,b ∈ T_t}  den_t(a,b)
```

where for gene tree T_t with k_t leaves:

- `den_t(a,b) = C2(k_t − 2)`  (total resolved quartets in T_t containing a and b)
- `num_t(a,b)` = number of those quartets where a and b are on the **same side** (see §2)

After normalization: `M[a][a] = 1`, `M[a][b] = M[b][a]`.  
For tree completion we expose `dist[a][b] = 1 − M[a][b]`.

---

## 2. Exact Pairwise Formula

### 2.1 Node-level S[u]

For an internal node u in tree T_t, define its **incident components**:

- one component per child subtree `c₁, c₂, …, c_d`
- one **parent-side** component of size `k_t − subLeafCount[u]` (0 for the root)

Then:

```
S[u] = C2(k_t − subLeafCount[u])  +  Σ_c  C2(subLeafCount[c])
```

### 2.2 Per-pair contribution

For pair (a, b) with `u = LCA_t(a, b)`:

Let `c_a` = child of u on the path toward a,  
    `c_b` = child of u on the path toward b.

```
num_t(a,b) = S[u] − C2(subLeafCount[c_a]) − C2(subLeafCount[c_b])
```

This counts quartets {a,b,x,y} where x and y come from the **same** incident component
that is neither `c_a` nor `c_b` — exactly the quartets supporting a,b on the same side.

**This formula is correct for binary trees AND polytomies.**

For binary trees it simplifies: `num_t(a,b) = C2(k_t − subLeafCount[u])`,
but the implementation uses the general formula to support polytomies.

---

## 3. GPU Architectural Advantage

The original CPU code **scatters** node contributions into many matrix cells:

- many write conflicts
- requires O(n²) GPU accumulator
- heavy atomic writes

The pairwise reformulation turns this into **pair-owned independent evaluation**:

- thread (da, db) owns pair (a0+da, b0+db) for the current tile
- loops over a batch of Δ trees, reading only **its** pair's LCA data
- writes only to its own tile accumulator — **no atomics**
- GPU VRAM stays **O(B² + Δ·n·log n)**, never O(n²)

---

## 4. Tiling and Batching Strategy

Same architecture as the distance matrix kernel:

- **B×B output tile**: B = min(n, ceil(sqrt(n·k))). Tile VRAM = O(B²).
- **Δ-tree batching**: Δ chosen so that Δ·(per-tree GPU bytes) fits the configured
  tree-data cap (`--gpu-sim-vram-cap-mb`, default 512 MiB) and currently free VRAM.
  If the option was not set explicitly, the large-N packed path may raise this
  ceiling to 8 GiB (still clamped to free VRAM) to avoid repeating every taxon
  pair across many small tree batches; the established dense path is unchanged.
- **Upper-triangle tiling**: only tiles with a0 ≤ b0 are processed; each pair
  is written once into the packed symmetric upper triangle.
- **Bounded host flattening**: when the established wide blocked layout itself
  would exceed Java's single-array element limit, its trees are flattened and
  submitted in original-order host batches of at most 1 GiB. Each batch adds to
  the same unnormalized numerator and denominator, and the matrix is normalized
  once after the final batch. Inputs whose compact or wide arrays already fit
  retain their established one-shot paths. After the completion/consensus phases
  have consumed a streamed matrix, its reference is dropped and one full-GC hint
  is issued before the allocation-heavy tripartition scan; this cleanup is not
  applied to one-shot runs.

For each tile:
1. Zero `numTile[B×B]` and `denTile[B×B]` on GPU
2. For each tree batch of Δ trees: upload → launch kernel → (accumulate across Δ trees)
3. Download tile → CPU folds each tile row slice into the packed triangle
   (`Get/SetDoubleArrayRegion`; no whole-accumulator pin or copy)

---

## 5. Per-Tree Preprocessing (CPU, parallel)

For each gene tree T_t, build:

### 5.1 Standard Euler tour + RMQ (shared with distance matrix)

DFS produces:
- `eulerDepths[tourLen]` — depth at each tour position
- `firstOcc[n]`           — first Euler position of each leaf taxon; −1 if absent
- `leafDepth[n]`          — depth of each leaf; −1 if absent (presence test)
- `eulerLen`              — actual emitted tour length (3·k_t − 2 for a strictly binary rooted tree)
- `sparseMin[LOG][tourLen]` — standard min-depth sparse table (left-biased)

### 5.2 Bridge-formula Euler payload arrays

For each node `v`, preprocessing computes its descendant-leaf count `s(v)` and
the root-to-node prefix `F(v)` used by the bridge formula. At each Euler position:

- `eulerF[pos]` stores `F(node at pos)`.
- At an internal node's INTERMEDIATE position,
  `eulerLeftChildS/F` and `eulerRightChildS/F` store the two child payloads.

For a binary-tree leaf pair `(a,b)`, the leftmost minimum-depth Euler position
between their first occurrences is the unique INTERMEDIATE visit of their LCA.

### 5.3 Compact argmin sparse table

The current bridge-formula implementation stores one unsigned 16-bit Euler
position per sparse cell:

```
sparseArgmin[lvl][pos] = leftmost minimum-depth Euler position
                         in [pos, pos + 2^lvl)
```

Level 0 stores `pos`. Higher levels compare the two half-interval minima using
the same `leftDepth <= rightDepth` tie rule and propagate the winning position.
The child-size/F payloads are not replicated at every sparse level; after an
RMQ chooses the exact Euler position, the kernel fetches them from the base
Euler arrays. This preserves the selected position exactly while reducing a
sparse cell from 22 bytes to 2 bytes.

### 5.4 Wide-tour blocked RMQ

The compact table is used unchanged while every Euler position fits in an
unsigned 16-bit value (tour length at most 65,536) and the flattened Java arrays
fit safely. Larger inputs are dispatched to a separate exact blocked RMQ:

- Euler positions are divided into blocks of 256.
- Nine in-block sparse levels store unsigned-byte offsets within each block.
- A small 32-bit sparse table stores argmin positions over whole blocks.
- A query combines the left partial block, zero or more whole middle blocks, and
  the right partial block. Depth ties always keep the earlier position, so the
  result is the same left-biased argmin as the compact table.
- Euler depths and child-subtree sizes are 32-bit on this path, avoiding a
  second signed-16-bit limit above 32,767 taxa.

The wide arrays use the exact maximum tour length rather than power-of-two
padding. The normal compact kernel and its memory layout are not changed.

---

## 6. GPU Kernel Logic

Thread (da, db) handles pair `a = a0+da`, `b = b0+db` for the current B×B tile.

```
for each tree t in current Δ-batch:
    // Presence check
    if leafDepth[t][a] < 0 || leafDepth[t][b] < 0: continue

    fa = firstOcc[t][a];  fb = firstOcc[t][b]
    l  = min(fa, fb);     r  = max(fa, fb)

    // O(1) left-biased RMQ
    k_lvl = 31 - clz(r - l + 1)
    l2    = r - (1 << k_lvl) + 1

    p_l = sparseArgmin[t][k_lvl][l]
    p_r = sparseArgmin[t][k_lvl][l2]
    p   = (eulerDepth[p_l] <= eulerDepth[p_r]) ? p_l : p_r

    leftS,leftF,rightS,rightF = base Euler payloads at p
    aS,aF = (fa <= fb) ? (leftS,leftF)   : (rightS,rightF)
    bS,bF = (fa <= fb) ? (rightS,rightF) : (leftS,leftF)

    // Accumulate
    den  = C2(kt[t] - 2)
    Z = kt[t] - aS - bS
    twoQD = (eulerF[fa] - aF) + (eulerF[fb] - bF)
            + (aS - 1)*Z + (bS - 1)*Z
    num = den - twoQD/2
    numTile[da * bB + db] += num
    denTile[da * bB + db] += den
```

No atomics — each (da,db) owns a unique cell.

---

## 7. Flat GPU Arrays (per-tree, flattened for Δ-batching)

| Array | Type | Size per tree | Notes |
|---|---|---|---|
| `eulerDepths`       | `short` | E_max    | depth at each Euler position |
| `eulerF`            | `double` | E_max   | F(node) at each position |
| `eulerLeftChildS`   | `short` | E_max    | left-child size at intermediates |
| `eulerLeftChildF`   | `double` | E_max   | left-child F at intermediates |
| `eulerRightChildS`  | `short` | E_max    | right-child size at intermediates |
| `eulerRightChildF`  | `double` | E_max   | right-child F at intermediates |
| `sparseArgmin`      | `char`/`uint16` | LOG×E_max | left-biased argmin Euler position |
| `firstOcc`          | `int`   | n         | first Euler position per leaf |
| `leafCount`         | `int`   | 1         | k_t per tree |

The bridge-formula implementation's base Euler payloads occupy 30 bytes per
position and the compact sparse table occupies 2 bytes per cell:

```
per-tree bytes = E_max × (30 + 2·LOG) + 4n + 4
```

The wide path replaces `sparseArgmin` with a byte micro table plus a small int
macro table and widens the three `short` Euler/child-size arrays to `int`.
For a 25,000-leaf binary tree, the exact Euler length is 74,998; with 1,000
trees, wide tree data is roughly 3–4 GiB instead of the roughly 12.8 GiB that a
power-of-two-padded full 32-bit sparse table would require.

---

## 8. Memory Complexity

### GPU VRAM (for one tile + Δ-tree batch)

```
VRAM = 2 · B² · 8          (numTile + denTile, double)
     + Δ · per_tree_bytes   (tree data on GPU)
```

With B = sqrt(n·k):  `B² = n·k`  →  tile VRAM = **O(n·k)**.  
Tree-batch VRAM = O(Δ·n·log n).  
Total GPU VRAM = **O(n·k + Δ·n·log n)** — no O(n²) term.

### CPU RAM

The matrix is symmetric, so every pair is stored exactly once: the upper
triangle (a ≤ b, n(n+1)/2 cells) in 512-MiB `double[]` segments addressed by
64-bit logical indices (`SegmentedDoubleArray`), for every n. Precision remains
`double`; no pair is sampled or omitted. For 50,000 taxa the triangle is
1,250,025,000 doubles, or 9.31 GiB.

Host footprint over the matrix lifetime (bytes, n taxa):

| Phase | Resident | Size |
|---|---|---|
| accumulation (CPU scatter or CUDA tile merge) | numerator + denominator triangles | 8·n(n+1)/2 ≈ 8 n² |
| after `normalize()` | similarity triangle (denominator released) | ≈ 4 n² |
| completion-anchor rows (`SortedRows`, only when some tree is incomplete) | one nearest-neighbour row per taxon | 2 n² (n ≤ 65,536) or 4 n² |
| UPGMA | the similarity triangle itself, consumed in place | +0 |

The nearest-neighbour rows are built before UPGMA runs, because they are the
only matrix-derived input tree completion reads (similarity values are never
read during insertion). When nothing reads the matrix after the guide tree is
built — every path except S3 consensus and the consensus verifiers — UPGMA
mutates the triangle in place instead of copying it, and the
`SimilarityMatrix` object becomes unreadable (any later `getSim` throws).
Otherwise UPGMA clusters a private copy (+4 n²).

The CUDA host merge never pins or copies a whole accumulator: for each tile it
reads and writes only the row slices the tile touches through
`Get/SetDoubleArrayRegion` (O(B) scratch per call). HotSpot's
`Get<Prim>ArrayElements` always copies the entire array into native memory, so
the former whole-matrix acquisition cost a further 16 n² bytes per host batch.

The matrix is O(n²) — unavoidable while UPGMA needs every pairwise value.
Preprocessing is O(k·n·log n) in total work. Its peak host storage is also
O(k·n·log n) while the selected flattened arrays fit; beyond Java's single-array
limit, bounded original-order wide streaming reduces the preprocessing peak to
O(Δ_host·n·log n) without changing the accumulated values.

Previously (through 46,340 taxa) the matrix was held as four dense n×n
`double[]` (numerator, denominator, similarity, distance; 32 n² bytes, never
released before Phase 5), JNI copied numerator and denominator per host batch
(+16 n²), UPGMA copied the triangle (+4 n²), and sorted rows were `int` (+4 n²):
roughly 48 n² at the GPU peak and 36 n² in Phase 5, against ≈ 8 n² and
≈ 2–4 n² now.

---

## 9. Normalization and Output

After all tiles are processed, normalization runs in place on the packed
numerator triangle and releases the packed denominator:

```java
num[i] = (den[i] > 0) ? num[i] / den[i] : 0.0;   // every triangle cell
num[index(a,a)] = 1.0;                           // diagonal
dist(a,b) = 1.0 - sim(a,b)                       // evaluated on demand
```

Symmetry is structural (one cell per unordered pair), and CUDA computes only
the upper triangle, so each pair is written once. Per-tree contributions are
integer-valued doubles, so the accumulated sums are exact whatever the tile,
batch or thread order; the normalized values are bit-identical to the former
dense n×n arrays (`PackedSimilarityParityTest` checks this against an in-test
dense reference).

---

## 10. Integration with Tree Completion

`TreeCompleter` (and the lineage-scoped `ScopedTreeCompleter`) read only
`SortedRows`: for each taxon, every other taxon ordered by ascending distance
with ASTRAL-MP's descending-taxon-ID tie break. The anchor of a missing taxon is
the first present taxon in its row. `SortedRowsBuilder` produces the identical
order with a primitive merge sort (the key (distance, −ID) is unique, so the
order does not depend on the sorting algorithm) and stores rows as unsigned
16-bit IDs whenever n ≤ 65,536.

In `Main.java` (S2/S4/S7 DISCO enrichment):
1. Build `SimilarityMatrix` from the retained single-copy DISCO components.
2. If any component needs completion, build `SortedRows` from the matrix.
3. Build the UPGMA guide, consuming the matrix in place unless S3 consensus or
   a consensus verifier still needs similarity values.
4. Complete components (immediately, or batch-wise in Phase 5 for the streamed
   S2/S4/S7 path) using the rows; release the rows afterwards.

This matches ASTRAL-MP's default behavior (similarity-guided completion).
