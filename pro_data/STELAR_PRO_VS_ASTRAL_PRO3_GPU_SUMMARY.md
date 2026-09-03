# STELAR-Pro vs ASTRAL-Pro3 benchmark

Each cell is the mean ± sample SD over 10 replicates (200 taxa and 1,000 gene trees per replicate). Wall time used 28 CPU threads. STELAR-Pro used strict CUDA execution with S1/I1 defaults; RF rate is standard unrooted RF / 2(n−3), where lower is better.

| Mean duplicates/species | SimPhy duplication rate | STELAR-Pro time (s) | STELAR-Pro RF rate | ASTRAL-Pro3 time (s) | ASTRAL-Pro3 RF rate |
|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.8e-10 | 26.226 ± 3.827 | 0.043 ± 0.030 | 67.309 ± 7.778 | 0.024 ± 0.014 |
| 1 | 2e-10 | 27.749 ± 6.028 | 0.095 ± 0.152 | 73.024 ± 4.908 | 0.024 ± 0.012 |
| 3 | 4.1e-10 | 29.715 ± 5.942 | 0.069 ± 0.067 | 70.678 ± 18.975 | 0.023 ± 0.016 |

Detailed per-replicate measurements: `pro_data/benchmark-stelar-pro-vs-astral-pro3-gpu/replicate-results.csv`.
