# S3, S4, updated S5 versus ASTRAL-Pro3: all ten 1-lb replicates

Dataset: `simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss`, with 200
taxa and 1,000 gene trees per replicate. Accuracy is standard unrooted RF
between the inferred and true single-copy species trees; RF rate = raw RF / 394
(lower is better). Runtime is full-process wall time. SD is sample standard
deviation across the ten replicates.

| Method | Mean RF rate ± SD | Mean runtime ± SD |
|:--|--:|--:|
| STELAR-Pro S3 | 3.55 ± 1.53% | 52.65 ± 24.33 s |
| STELAR-Pro S4 | 3.55 ± 1.53% | 59.46 ± 29.44 s |
| STELAR-Pro S5, overall optimization | 3.65 ± 1.67% | 160.11 ± 50.78 s |
| ASTRAL-Pro3 | **2.44 ± 1.19%** | 73.02 ± 4.91 s |

## Per-replicate RF (raw count / 394)

| Replicate | S3 | S4 | Updated S5 | ASTRAL-Pro3 |
|:--|--:|--:|--:|--:|
| R1 | 16 | 16 | 16 | 12 |
| R2 | 6 | 6 | 6 | 2 |
| R3 | 20 | 20 | 20 | 12 |
| R4 | 22 | 22 | 22 | 18 |
| R5 | 12 | 12 | 12 | 8 |
| R6 | 6 | 6 | 6 | 4 |
| R7 | 20 | 20 | 24 | 10 |
| R8 | 12 | 12 | 12 | 8 |
| R9 | 18 | 18 | 18 | 14 |
| R10 | 8 | 8 | 8 | 8 |

## Per-replicate wall time (seconds)

| Replicate | S3 | S4 | Updated S5 | ASTRAL-Pro3 |
|:--|--:|--:|--:|--:|
| R1 | 74.58 | 87.72 | 192.12 | 71.62 |
| R2 | 21.40 | 23.57 | 89.94 | 66.84 |
| R3 | 86.06 | 103.05 | 242.65 | 71.63 |
| R4 | 70.72 | 84.27 | 206.93 | 82.16 |
| R5 | 39.02 | 46.85 | 143.71 | 69.37 |
| R6 | 28.04 | 30.34 | 109.95 | 67.41 |
| R7 | 67.46 | 78.87 | 184.34 | 73.80 |
| R8 | 43.47 | 43.86 | 136.42 | 74.49 |
| R9 | 73.23 | 73.62 | 192.30 | 73.31 |
| R10 | 22.54 | 22.44 | 102.78 | 79.62 |

ASTRAL-Pro3 has lower RF than all three STELAR-Pro paths in R1–R9 and ties
them in R10. S3 and S4 have identical unrooted topologies in all ten
replicates. The updated S5 path does not change the original S5 topology in
any replicate: each run added four new greedy proposals to the four initial
ones, retained all 1,000 S1 gene trees for the full triplet DP, and stopped
after its first added pass because the rooted tree was unchanged. S5 differs
in accuracy from S3/S4 only in R7, where its RF is 24 rather than 20.

The S3/S4 timings and ASTRAL-Pro3 timings are from earlier runs on this
machine; ASTRAL-Pro3 used 28 CPU threads and seed 42, while S3/S4 used CUDA.
The updated S5 runs used strict CUDA, 28 CPU threads, four proposals per pass,
and seed 233. These are recorded wall times from separate runs, not a
simultaneous controlled timing trial; ASTRAL-Pro3's default output work may
also differ from STELAR-Pro's. RF values were independently recomputed from
the saved trees. The initial R3 updated-S5 attempt hit a 240-second cutoff at
the end of GPU scoring; the table uses the successful 242.65-second retry
with a 300-second cutoff, and excludes the failed attempt's time.

Sources: S3 under `benchmark-stelar-pro-s3-vs-s2-vs-astral-pro3-gpu-dupmean-1/`,
S4 under `benchmark-stelar-pro-s4-vs-s3-gpu-dupmean-1/` plus four missing-run
retries under `benchmark-stelar-pro-s4-retry-gpu-dupmean-1/`, and updated S5
trees/times under `benchmark-stelar-pro-s5-adaptive-unbiased-gpu-dupmean-1/`.
The earlier ASTRAL-Pro3 results are in
`benchmark-stelar-pro-vs-astral-pro3-gpu/replicate-results.csv` and its
per-replicate tree/JSON files. The R3 S5 retry time is
`R3/S5.retry300.wall-seconds.txt`; the original `R3/S5.wall-seconds.txt`
records the failed four-minute attempt.
