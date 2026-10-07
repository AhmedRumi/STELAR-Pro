# STELAR-Pro S3/S4/S5: first five 1-lb replicates

Dataset: `simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss`, R1–R5.
RF is standard unrooted RF against each replicate's true single-copy species
tree, divided by 394 (= 2 × (200 − 3)). Times are recorded wall seconds.
S3 and S4 trees/times are from their existing GPU benchmark runs; S5 was run
with `--search-space S5 --s5-rounds 4 --s5-seed 233 --gpu-strict --threads 28`.

| Replicate | S3 RF | S4 RF | S5 RF | S3 time (s) | S4 time (s) | S5 time (s) |
|:--|--:|--:|--:|--:|--:|--:|
| R1 | 16/394 (4.06%) | 16/394 (4.06%) | 16/394 (4.06%) | 74.58 | 87.72 | 106.27 |
| R2 | 6/394 (1.52%) | 6/394 (1.52%) | 6/394 (1.52%) | 21.40 | 23.57 | 56.02 |
| R3 | 20/394 (5.08%) | 20/394 (5.08%) | 20/394 (5.08%) | 86.06 | 103.05 | 121.69 |
| R4 | 22/394 (5.58%) | 22/394 (5.58%) | 22/394 (5.58%) | 70.72 | 84.27 | 106.17 |
| R5 | 12/394 (3.05%) | 12/394 (3.05%) | 12/394 (3.05%) | 39.02 | 46.85 | 78.26 |
| **Mean** | **3.86%** | **3.86%** | **3.86%** | **58.36** | **69.09** | **93.68** |

S3 and S4 have the same unrooted topology in all five replicates. S5 has the
same topology as S3/S4 in R1, R2, R3, and R5; R4 has the same RF to truth but
differs from S3/S4 by four RF split counts. All five S5 runs completed with
four triplet-DP backbones and GPU weight scoring; no timeout or exception was
logged. S5 outputs, logs, and wall times are in
`benchmark-stelar-pro-s5-gpu-dupmean-1/R1` through `R5`.
