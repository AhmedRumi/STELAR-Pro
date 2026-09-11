# STELAR-Pro S1/S2 versus ASTRAL-Pro3: 1-duplication dataset

Ten replicates; 200 taxa and 1,000 gene trees per replicate. Both methods used
28 CPU threads. STELAR-Pro S2 used strict CUDA on the RTX 4060. Times are full
process wall-clock measurements. RF is normalized unrooted RF (`raw RF / 394`);
lower is better.

| Method | Mean RF rate ± SD | Mean runtime (s) ± SD |
|---|---:|---:|
| STELAR-Pro S1 | 0.095431 ± 0.151904 | **27.749 ± 6.028** |
| STELAR-Pro S2 | 0.035533 ± 0.015322 | **57.580 ± 24.552** |
| ASTRAL-Pro3 | **0.024365 ± 0.011917** | 73.024 ± 4.908 |

| Replicate | S2 RF | ASTRAL-Pro3 RF | S2 time (s) | ASTRAL-Pro3 time (s) |
|:---:|---:|---:|---:|---:|
| R1 | 0.040609 | **0.030457** | 79.320 | **71.616** |
| R2 | 0.015228 | **0.005076** | **24.250** | 66.838 |
| R3 | 0.050761 | **0.030457** | 90.620 | **71.635** |
| R4 | 0.055838 | **0.045685** | **74.790** | 82.161 |
| R5 | 0.030457 | **0.020305** | **44.370** | 69.374 |
| R6 | 0.015228 | **0.010152** | **33.010** | 67.408 |
| R7 | 0.050761 | **0.025381** | **72.340** | 73.797 |
| R8 | 0.030457 | **0.020305** | **49.540** | 74.486 |
| R9 | 0.045685 | **0.035533** | 79.900 | **73.308** |
| R10 | **0.020305** | **0.020305** | **27.660** | 79.618 |

S2 reduced mean RF by 62.8% relative to S1, but took 2.07× as long.
ASTRAL-Pro3 had lower RF than S2 in 9 replicates and tied S2 in R10. Its mean RF was
0.011168 lower (31.4% below S2). S2 was faster in 7 replicates and used 21.2%
less wall time on average; equivalently, ASTRAL-Pro3 took 1.27× as long based
on the ratio of mean runtimes.

S2 outputs, logs, and recorded wall times are under
`pro_data/benchmark-stelar-pro-s2-vs-astral-pro3-gpu-dupmean-1/`. The existing
ASTRAL-Pro3 results were read from
`pro_data/benchmark-stelar-pro-vs-astral-pro3-gpu/replicate-results.csv`.
