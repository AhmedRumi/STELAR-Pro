# STELAR-Pro S2/S3 versus ASTRAL-Pro3: 1-duplication dataset

Ten replicates; 200 taxa and 1,000 gene trees per replicate. STELAR-Pro used
strict CUDA on the RTX 4060 with 28 CPU threads. Times are full process wall
times. RF is normalized unrooted RF (`raw RF / 394`); lower is better.

| Method | Mean RF rate ± SD | Mean runtime (s) ± SD |
|---|---:|---:|
| STELAR-Pro S2 | 0.035533 ± 0.015322 | 57.580 ± 24.552 |
| STELAR-Pro S3 | 0.035533 ± 0.015322 | 52.652 ± 24.325 |
| ASTRAL-Pro3 | **0.024365 ± 0.011917** | 73.024 ± 4.908 |

| Replicate | S2 RF | S3 RF | ASTRAL-Pro3 RF | S2 time (s) | S3 time (s) | ASTRAL-Pro3 time (s) | New S3 clusters |
|:---:|---:|---:|---:|---:|---:|---:|---:|
| R1 | 0.040609 | 0.040609 | **0.030457** | 79.320 | 74.580 | 71.616 | 0 |
| R2 | 0.015228 | 0.015228 | **0.005076** | 24.250 | 21.400 | 66.838 | 0 |
| R3 | 0.050761 | 0.050761 | **0.030457** | 90.620 | 86.060 | 71.635 | 0 |
| R4 | 0.055838 | 0.055838 | **0.045685** | 74.790 | 70.720 | 82.161 | 48 |
| R5 | 0.030457 | 0.030457 | **0.020305** | 44.370 | 39.020 | 69.374 | 0 |
| R6 | 0.015228 | 0.015228 | **0.010152** | 33.010 | 28.040 | 67.408 | 0 |
| R7 | 0.050761 | 0.050761 | **0.025381** | 72.340 | 67.460 | 73.797 | 4 |
| R8 | 0.030457 | 0.030457 | **0.020305** | 49.540 | 43.470 | 74.486 | 5 |
| R9 | 0.045685 | 0.045685 | **0.035533** | 79.900 | 73.230 | 73.308 | 0 |
| R10 | **0.020305** | **0.020305** | **0.020305** | 27.660 | 22.540 | 79.618 | 0 |

S3 synthesized 57 clusters beyond S2 (48 in R4, 4 in R7, and 5 in R8),
which produced 813 additional DP splits. None changed the inferred topology,
so S2 and S3 have identical RF on every replicate. ASTRAL-Pro3 has 31.4% lower
mean RF than both STELAR-Pro searches.

The S2 and ASTRAL-Pro3 timings are the previously recorded runs; S3 was run
later rather than in an interleaved timing experiment. Therefore the observed
8.6% lower S3 wall time is run-to-run/system variation, not an algorithmic
speedup: S3 strictly performs S2 plus consensus construction and refinement.

S3 outputs, logs, and wall times are under
`pro_data/benchmark-stelar-pro-s3-vs-s2-vs-astral-pro3-gpu-dupmean-1/`.
S2 results are under
`pro_data/benchmark-stelar-pro-s2-vs-astral-pro3-gpu-dupmean-1/`; ASTRAL-Pro3
values come from
`pro_data/benchmark-stelar-pro-vs-astral-pro3-gpu/replicate-results.csv`.
