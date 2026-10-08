# S6 versus S3–S5 and ASTRAL-Pro3: 200 taxa, duplication mean 1

Standard unrooted RF is measured against each true single-copy species tree (maximum RF 394). Lower is better. Wall times include full preprocessing and inference. S6 used strict CUDA, 28 CPU threads, 4 greedy trees per DP pass, seed 233, and a 600-second per-replicate cap. Comparator times come from separate earlier runs on this machine, not a simultaneous timing trial.

| Replicate | S6 RF /394 | S6 time (s) | S5 RF /394 | S5 time (s) | A-Pro3 RF /394 | A-Pro3 time (s) |
|:--|--:|--:|--:|--:|--:|--:|
| R1 | 16 | 156.67 | 16 | 192.12 | 12 | 71.62 |
| R2 | 6 | 39.25 | 6 | 89.94 | 2 | 66.84 |
| R3 | 20 | 189.91 | 20 | 242.65 | 12 | 71.63 |
| R4 | 22 | 153.91 | 22 | 206.93 | 18 | 82.16 |
| R5 | 12 | 82.54 | 12 | 143.71 | 8 | 69.37 |
| R6 | 6 | 60.39 | 6 | 109.95 | 4 | 67.41 |
| R7 | 20 | 147.75 | 24 | 184.34 | 10 | 73.80 |
| R8 | 12 | 93.68 | 12 | 136.42 | 8 | 74.49 |
| R9 | 18 | 149.42 | 18 | 192.30 | 14 | 73.31 |
| R10 | 8 | 40.82 | 8 | 102.78 | 8 | 79.62 |

| Method | Mean RF rate ± SD | Mean wall time ± SD |
|:--|--:|--:|
| S6 | 3.55 ± 1.53% | 111.43 ± 54.49 s |
| S3 | 3.55 ± 1.53% | 52.65 ± 24.33 s |
| S4 | 3.55 ± 1.53% | 61.24 ± 28.62 s |
| S5 | 3.65 ± 1.67% | 160.11 ± 50.78 s |
| ASTRAL-Pro3 | 2.44 ± 1.19% | 73.02 ± 4.91 s |

Pairwise output-tree check: S6 versus S3 has RF 0 in 10/10 replicates; versus S4, 10/10. S6 differs from S5 in R4: RF 4, R7: RF 6.

The per-replicate table emphasizes S6 against its direct S5 predecessor and ASTRAL-Pro3; S3/S4 are included in the aggregate. S4 uses the revised single-triplet-matrix run for R6–R10 and the earlier S4 run for R1–R5.

S6 proposals use only the rooted-triplet support of the original tagged gene trees; S1 candidates are retained, and all final DP weights use the same original genes.
