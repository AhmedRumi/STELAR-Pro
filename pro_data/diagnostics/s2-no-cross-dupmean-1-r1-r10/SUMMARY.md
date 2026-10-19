# S2 without cross-tree transitions: 200-taxon duplication-mean-1 dataset

Ten replicates (R1–R10), 1,000 gene trees each. This temporary **S2 local-only**
variant retains S1 candidates, DISCO decomposition, completion, UPGMA, and
local DP transitions; it skips only Phase 5b cross-tree transitions. The
original rooted/tagged multicopy trees still provide the triplet weights.
No change was made to the main STELAR-Pro source tree.

Both variants were run sequentially with the **same temporary build**, strict
CUDA, 28 CPU threads, and the default GPU smaller-side-traversal weight
backend. Wall time is full process time (`/usr/bin/time`). RF is unrooted RF
against each replicate's true single-copy species tree, divided by 394; lower
is better. “Splits” means unique DP transitions before reachability pruning,
not candidate clusters. The full-S2 rerun reproduced all ten previously saved
full-S2 topologies exactly (pairwise RF 0).

| Replicate | RF local-only | RF full S2 | Time local-only (s) | Time full (s) | Splits local-only | Splits full |
|:--|--:|--:|--:|--:|--:|--:|
| R1 | 16/394 | 16/394 | 19.63 | 73.37 | 72,649 | 471,101 |
| R2 | 6/394 | 6/394 | 9.06 | 18.44 | 32,382 | 179,684 |
| R3 | 22/394 | 20/394 | 19.85 | 84.30 | 78,057 | 581,620 |
| R4 | 26/394 | 22/394 | 23.12 | 68.62 | 80,388 | 409,515 |
| R5 | 12/394 | 12/394 | 11.21 | 37.81 | 49,503 | 308,395 |
| R6 | 6/394 | 6/394 | 8.64 | 27.20 | 45,273 | 300,214 |
| R7 | 22/394 | 20/394 | 18.67 | 65.92 | 72,733 | 408,403 |
| R8 | 12/394 | 12/394 | 12.38 | 41.59 | 50,991 | 310,808 |
| R9 | 18/394 | 18/394 | 17.90 | 72.02 | 72,284 | 495,910 |
| R10 | 6/394 | 8/394 | 7.86 | 21.43 | 33,888 | 190,193 |
| **Mean** | **3.706%** | **3.553%** | **14.83** | **51.07** | **58,815** | **365,584** |

Removing cross-tree transitions reduced mean unique DP splits by **306,770
(83.9%)** and mean wall time by **36.24 s (71.0%; 3.44× speedup)**. RF was
unchanged in six replicates, worse in three (R3, R4, R7), and better in one
(R10). Mean RF increased by **0.152 percentage points** (0.6 raw RF splits).
The extracted candidate clusters are unchanged in every matched pair; only
the cross-tree DP transitions are omitted.

Count audit: an earlier comparison mistakenly used the DISCO-completion
progress count before the UPGMA proposal was emitted. UPGMA adds one R3 cluster
and two R4 clusters, bringing streamed X to 56,154 and 65,068 respectively.
The final streamed X count matches the earlier materialized run in **all ten**
replicates. For R3 and R4, materialized-versus-streamed completed-tree and
sorted cluster dumps are byte-for-byte identical. The S2 progress log now
labels DISCO-only counts and reports the UPGMA contribution separately; a
runtime check asserts that the reported final count matches X.
The runtime comparison is a measured difference between separately timed runs
on the same GPU and build, not a hardware-independent speedup guarantee.

The local-only and full-S2 unique split totals differ by exactly the number
of cross-tree splits reported in each full log. Reachability pruning removed
14 local-only and 22 full-S2 splits in R3; all other splits were scored.
The configuration header's generic S2 description still mentions cross-tree,
but `Search mode: local` and the absence of Phase 5b in the local-only logs
verify that cross-tree transitions were skipped.

Per-replicate output trees and logs are in `R1/` through `R10/`. The run
settings above and `results.csv` provide the numerical record.
