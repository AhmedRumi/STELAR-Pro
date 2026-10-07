# Preliminary adaptive S5 versus S3/S4: 1-lb dataset, R6–R10

This is a superseded diagnostic run: its later-round subset always included
the placement-anchor taxon. A corrected run that samples from all taxa is
being saved separately under `benchmark-stelar-pro-s5-adaptive-unbiased-gpu-dupmean-1/`.
Do not treat the numbers below as the final adaptive S5 comparison.

Dataset: `simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss`.
RF is standard unrooted RF against each true 200-taxon species tree,
normalized by 394. Times are recorded wall seconds. S3/S4 and old S5
columns use prior benchmark trees; updated S5 was run with GPU, 28 CPU
threads, `--s5-rounds 4`, `--s5-seed 233`, and a 240-second per-run cutoff.

| Replicate | S3 RF | S4 RF | Old S5 RF | Adaptive S5 RF | S3 time | S4 time | Old S5 time | Adaptive S5 time |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| R6 | 6/394 (1.52%) | 6/394 (1.52%) | 6/394 (1.52%) | 6/394 (1.52%) | 28.04 s | 30.34 s | 66.58 s | 109.86 s |
| R7 | 20/394 (5.08%) | 20/394 (5.08%) | 24/394 (6.09%) | 24/394 (6.09%) | 67.46 s | 78.87 s | 102.07 s | 187.28 s |
| R8 | 12/394 (3.05%) | 12/394 (3.05%) | 12/394 (3.05%) | 12/394 (3.05%) | 43.47 s | 43.86 s | 78.45 s | 136.39 s |
| R9 | 18/394 (4.57%) | 18/394 (4.57%) | 18/394 (4.57%) | 18/394 (4.57%) | 73.23 s | 73.62 s | 105.72 s | 193.01 s |
| R10 | 8/394 (2.03%) | 8/394 (2.03%) | 8/394 (2.03%) | 8/394 (2.03%) | 22.54 s | 22.44 s | 62.68 s | 102.51 s |
| **Mean** | **3.25%** | **3.25%** | **3.45%** | **3.45%** | **46.95 s** | **49.83 s** | **83.10 s** | **145.81 s** |

Each updated S5 run executed four initial two-step greedy proposals, a full
STELAR triplet DP, four additional proposals guided by that DP tree, and a
second full DP retaining all original S1 candidates and all eight S5
proposals. In every replicate the rooted tree and triplet score were unchanged
after the first additional pass, so the search stopped. Old versus updated S5
unrooted RF is zero in all five replicates. All runs completed below 240 s.

The bundled ASTRAL-Pro3 R6 cross-check (`-r 4 -s 4 --seed 233`) also ran one
additional four-search batch and stopped after no score improvement. Its
quartet-objective RF to truth was 4/394; updated S5's triplet-objective RF was
6/394. These objectives and their stopping tests differ: S5 follows the
ASTRAL-Pro2 supplement's rooted-tree-stability rule with at most five added
passes, whereas bundled ASTRAL-Pro3 uses score improvement and a larger cap.

Updated S5 trees, logs, and wall times are under
`benchmark-stelar-pro-s5-adaptive-gpu-dupmean-1/R6` through `R10`.
