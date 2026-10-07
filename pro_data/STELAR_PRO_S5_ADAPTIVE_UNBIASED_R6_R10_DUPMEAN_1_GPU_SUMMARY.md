# Updated S5 overall search: 1-lb dataset, R6–R10

Dataset: `simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss` (200 taxa,
1,000 gene trees per replicate). RF is standard unrooted RF against the true
single-copy species tree, divided by 394; lower is better. Runtime is wall time.
All updated S5 runs used GPU (`--gpu-strict`), 28 CPU threads, four initial
greedy rounds, four new greedy proposals per adaptive pass, seed 233, and a
hard 240-second cutoff. All five completed below the cutoff.

| Replicate | S3 RF | S4 RF | Original S5 RF | Updated S5 RF | S3 time | S4 time | Original S5 time | Updated S5 time |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| R6 | 6/394 (1.52%) | 6/394 (1.52%) | 6/394 (1.52%) | 6/394 (1.52%) | 28.04 s | 30.34 s | 66.58 s | 109.95 s |
| R7 | 20/394 (5.08%) | 20/394 (5.08%) | 24/394 (6.09%) | 24/394 (6.09%) | 67.46 s | 78.87 s | 102.07 s | 184.34 s |
| R8 | 12/394 (3.05%) | 12/394 (3.05%) | 12/394 (3.05%) | 12/394 (3.05%) | 43.47 s | 43.86 s | 78.45 s | 136.42 s |
| R9 | 18/394 (4.57%) | 18/394 (4.57%) | 18/394 (4.57%) | 18/394 (4.57%) | 73.23 s | 73.62 s | 105.72 s | 192.30 s |
| R10 | 8/394 (2.03%) | 8/394 (2.03%) | 8/394 (2.03%) | 8/394 (2.03%) | 22.54 s | 22.44 s | 62.68 s | 102.78 s |
| **Mean** | **3.25%** | **3.25%** | **3.45%** | **3.45%** | **46.95 s** | **49.83 s** | **83.10 s** | **145.16 s** |

The initial S5 search produced four greedy trees and ran the full STELAR
triplet DP. Each updated run then made one additional pass: it sampled four
new greedy proposals guided by the preceding STELAR-DP tree, accumulated all
eight proposals alongside the original 1,000 S1 gene trees, and reran the full
STELAR DP. The rooted tree was unchanged, so it stopped; the configured limit
is five additional passes. All five updated trees have unrooted RF 0 against
their original S5 counterparts, so the extra pass increased runtime without
improving RF on these replicates. R7 remains worse than S3/S4.

Cross-check: bundled ASTRAL-Pro3 on R6 with `-r 4 -s 4 --seed 233` also ran one
additional four-search batch and stopped after no score improvement. Its RF
was 4/394, versus S5's 6/394. The stopping rules differ: this S5 follows the
ASTRAL-Pro2 supplement's tree-stability/five-pass schedule, but the bundled
ASTRAL-Pro3 source continues only while its quartet score improves, up to 20
additional batches. S5 always uses the STELAR rooted-triplet objective for DP.

Corrected updated S5 trees and wall times are under
`benchmark-stelar-pro-s5-adaptive-unbiased-gpu-dupmean-1/R6` through `R10`.
The earlier `STELAR_PRO_S5_ADAPTIVE_R6_R10_DUPMEAN_1_GPU_SUMMARY.md` is a
superseded diagnostic run with a fixed subset anchor.
