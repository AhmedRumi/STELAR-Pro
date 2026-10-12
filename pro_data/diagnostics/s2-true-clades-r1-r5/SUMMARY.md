# S2 true-clade oracle diagnostic: duplication mean 1, R1–R5

**Result:** adding the rooted true species tree as a candidate-only source did not change any of the first five S2 outputs. Mean normalized unrooted RF stayed **3.86%** (76 raw RF units across five 200-taxon trees, each with maximum RF 394).

| Replicate | S2 RF | S2 + true tree RF | New clusters | New full-DP splits | Verification |
|:--|--:|--:|--:|--:|:--|
| R1 | 16/394 (4.06%) | 16/394 (4.06%) | 0 | 0 | Full inference; output byte-identical |
| R2 | 6/394 (1.52%) | 6/394 (1.52%) | 0 | 0 | Full inference; output byte-identical |
| R3 | 20/394 (5.08%) | 20/394 (5.08%) | 0 | 0 | Full inference; output byte-identical |
| R4 | 22/394 (5.58%) | 22/394 (5.58%) | 3 | 35 | Full inference; output byte-identical |
| R5 | 12/394 (3.05%) | 12/394 (3.05%) | 0 | 0 | Exact search-space/order equivalence; full weight run omitted |
| **Mean** | **3.86%** | **3.86%** | — | — | **0 improvement** |

Method: an isolated experimental `Main.class` appended each replicate's `s_tree.trees` after S2's original, completed-DISCO, and UPGMA candidate trees. Its internal nodes were treated as speciation candidate nodes. The original 1,000 rooted/tagged multicopy gene trees remained the *only* weight evidence. S2's full cross-tree transition stage and ordinary STELAR-Pro triplet DP were retained. Baseline RF uses the saved S2 outputs in `pro_data/benchmark-stelar-pro-s2-vs-astral-pro3-gpu-dupmean-1/R*/S2.tre`; the oracle outputs for R1–R4 and paired candidate-space logs are in this directory. The GPU was unavailable on this host, so the full oracle checks used the unchanged default CPU intersection scorer. These runtimes are **not** paired speed comparisons with the earlier GPU benchmark.

The candidate-space comparisons are exact: unique clusters changed by 0, 0, 0, **+3**, 0 in R1–R5; final full-DP splits changed by 0, 0, 0, **+35**, 0. R1 and R3 each gained one *local* split that was already present after cross-tree expansion; their full runs checked the possible tie-order effect. R5 gained neither a cluster nor a local split. Because the true tree was appended last and every R5 proposal was already present, its transition order, gene partitions, and weights are unchanged; a second full weight run could only repeat the saved output. The original-gene partition counts matched exactly between conditions in every replicate.

R4's true tree became feasible after its three missing clusters were added, but DP still chose the original output. On the original gene evidence, that output scores **570,193,085**, above the true tree's **545,940,960**. R3 similarly scores **855,673,912** for the selected S2 tree versus **850,722,501** for truth. Thus missing true clades are not sufficient to explain the S2-versus-ASTRAL-Pro RF gap on these five replicates; the triplet objective can prefer a different topology even when truth is available. This is an oracle diagnostic, **not** a legitimate inference option.

The temporary source edits were removed before the runs; the experiment-only class lived outside the repository. The repository's `src/` and `build/` outputs remain unchanged. A separate attempted CPU tree-walk speed-up was excluded after a small multicopy parity check disagreed with the default scorer (reported scores 27 versus 91); it did not contribute to the numbers above.
