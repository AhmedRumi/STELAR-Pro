# STELAR-Pro S3/S4/S5: all ten 1-lb replicates

Dataset: `simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss` (200 taxa,
1000 gene trees per replicate). RF is standard unrooted RF versus each true
single-copy species tree; maximum RF is 394. Times are recorded wall seconds.
S5 used GPU, 28 CPU threads, four rounds, and seed 233. R6–R10 S5 and the
new S4 retries had a hard 240-second limit per run; all completed below it.

| Replicate | S3 RF | S4 RF | S5 RF | S3 time (s) | S4 time (s) | S5 time (s) |
|:--|--:|--:|--:|--:|--:|--:|
| R1 | 16/394 (4.06%) | 16/394 (4.06%) | 16/394 (4.06%) | 74.58 | 87.72 | 106.27 |
| R2 | 6/394 (1.52%) | 6/394 (1.52%) | 6/394 (1.52%) | 21.40 | 23.57 | 56.02 |
| R3 | 20/394 (5.08%) | 20/394 (5.08%) | 20/394 (5.08%) | 86.06 | 103.05 | 121.69 |
| R4 | 22/394 (5.58%) | 22/394 (5.58%) | 22/394 (5.58%) | 70.72 | 84.27 | 106.17 |
| R5 | 12/394 (3.05%) | 12/394 (3.05%) | 12/394 (3.05%) | 39.02 | 46.85 | 78.26 |
| R6 | 6/394 (1.52%) | 6/394 (1.52%) | 6/394 (1.52%) | 28.04 | 30.34 | 66.58 |
| R7 | 20/394 (5.08%) | 20/394 (5.08%) | 24/394 (6.09%) | 67.46 | 78.87 | 102.07 |
| R8 | 12/394 (3.05%) | 12/394 (3.05%) | 12/394 (3.05%) | 43.47 | 43.86 | 78.45 |
| R9 | 18/394 (4.57%) | 18/394 (4.57%) | 18/394 (4.57%) | 73.23 | 73.62 | 105.72 |
| R10 | 8/394 (2.03%) | 8/394 (2.03%) | 8/394 (2.03%) | 22.54 | 22.44 | 62.68 |
| **Mean R1–R10** | **3.55%** | **3.55%** | **3.65%** | **52.65** | **59.46** | **88.39** |

S3 and S4 infer identical unrooted topologies in all ten replicates. S5
matches that topology in eight; R4 has the same RF to truth but a different
topology (pairwise RF 4), while R7 is worse by four RF counts (pairwise RF 6).

S3 outputs/times and S4 R1–R5/R7 outputs/times are from the earlier GPU
benchmarks. S4 R6/R8/R9/R10 were rerun separately under
`benchmark-stelar-pro-s4-retry-gpu-dupmean-1/` because the earlier S4 outputs
were missing. S5 R1–R10 outputs, logs, and times are under
`benchmark-stelar-pro-s5-gpu-dupmean-1/`. The new bounded runs used
`--gpu-strict --threads 28`; S5 additionally used `--s5-rounds 4 --s5-seed 233`.
