# Matched S3 versus S4 GPU benchmark: 200 taxa, duplication mean 1

Ten replicates, 1,000 gene trees each. Both methods used this build, strict CUDA, 28 CPU threads, 8 GiB Java heap, and a 600-s per-run cap. Odd replicates ran S3 first; even replicates ran S4 first. Wall time includes preprocessing and inference. RF is standard unrooted RF against the true single-copy species tree (maximum 394; lower is better).

Build fingerprints: `Main.class` `d7eaa52d4ec70f21`, `libstelar_pro_sim.so` `9ddfe5428301eca6`.

| Replicate | First | S3 time (s) | S4 time (s) | S4 − S3 (s) | S3 RF | S4 RF | S3↔S4 RF |
|:--|:--|--:|--:|--:|--:|--:|--:|
| R1 | S3 | 82.94 | 81.56 | -1.39 | 16 | 16 | 0 |
| R2 | S4 | 23.20 | 22.96 | -0.24 | 6 | 6 | 0 |
| R3 | S3 | 97.04 | 96.19 | -0.86 | 20 | 20 | 0 |
| R4 | S4 | 79.58 | 78.89 | -0.69 | 22 | 22 | 0 |
| R5 | S3 | 44.24 | 44.38 | 0.14 | 12 | 12 | 0 |
| R6 | S4 | 33.24 | 33.24 | 0.00 | 6 | 6 | 0 |
| R7 | S3 | 75.30 | 75.04 | -0.25 | 20 | 20 | 0 |
| R8 | S4 | 48.66 | 48.85 | 0.20 | 12 | 12 | 0 |
| R9 | S3 | 82.22 | 82.71 | 0.49 | 18 | 18 | 0 |
| R10 | S4 | 26.46 | 26.16 | -0.30 | 8 | 8 | 0 |

| Measure | S3 | S4 | S4 − S3 |
|:--|--:|--:|--:|
| Mean wall time ± SD | 59.29 ± 27.03 s | 59.00 ± 26.80 s | -0.29 s |
| Mean RF rate | 3.55% | 3.55% | — |
| Mean DISCO/completion/matrix/UPGMA phase | 0.88 s | 0.66 s | -0.22 s |
| Mean Triplet weight-scoring phase | 50.88 s | 50.98 s | +0.10 s |

Paired S4 − S3 wall-time difference: -0.29 ± 0.56 s (mean ± SD; median -0.25 s; odd/S3-first -0.37 s; even/S4-first -0.21 s).
S3 and S4 output trees have unrooted RF 0 in 10/10 replicates.

S4 uses one GPU rooted-triplet similarity matrix for completion, UPGMA, and consensus; S3 uses its quartet matrix. Both use the same original tagged multicopy genes for final STELAR-Pro triplet weights.
