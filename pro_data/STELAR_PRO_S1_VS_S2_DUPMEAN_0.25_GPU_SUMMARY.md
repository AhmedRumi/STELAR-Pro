# STELAR-Pro S1 versus S2: 0.25 duplication dataset

Ten replicates; 200 taxa and 1,000 gene trees per replicate. Both modes used
28 CPU threads and strict CUDA on the RTX 4060. Runtime is STELAR-Pro's internal
end-to-end time, including preprocessing and inference. RF is normalized
unrooted RF (`raw RF / 394`); lower is better.

| Metric | S1 mean ± SD | S2 mean ± SD | S2 change |
|---|---:|---:|---:|
| RF rate | 0.042640 ± 0.029811 | 0.034010 ± 0.017429 | -0.008629 (-20.2%) |
| Runtime (s) | 23.440 ± 3.882 | 52.552 ± 17.130 | 2.20× slower |
| Unique candidate clusters (X) | 35,904.0 ± 12,377.6 | 35,956.2 ± 12,383.9 | +52.2 (+0.15%) |
| Final DP candidate splits scored | 54,041.6 ± 13,297.6 | 338,040.3 ± 89,844.3 | +283,998.7 (6.25× total) |

| Replicate | S1 RF | S2 RF | S1 time (s) | S2 time (s) | Unique X: S1 → S2 | DP splits: S1 → S2 |
|:---:|---:|---:|---:|---:|---:|---:|
| R1 | 0.071066 | 0.071066 | 27.169 | 76.132 | 51,024 → 51,084 (+60) | 72,171 → 466,949 (+394,778) |
| R2 | 0.020305 | 0.015228 | 20.956 | 50.552 | 33,514 → 33,572 (+58) | 50,740 → 343,020 (+292,280) |
| R3 | 0.040609 | 0.045685 | 22.642 | 53.378 | 34,645 → 34,695 (+50) | 56,265 → 363,411 (+307,146) |
| R4 | 0.101523 | 0.035533 | 28.455 | 81.032 | 53,093 → 53,166 (+73) | 71,780 → 493,977 (+422,197) |
| R5 | 0.025381 | 0.025381 | 18.883 | 32.585 | 22,430 → 22,477 (+47) | 37,520 → 219,840 (+182,320) |
| R6 | 0.025381 | 0.025381 | 19.705 | 38.905 | 24,777 → 24,827 (+50) | 41,687 → 253,700 (+212,013) |
| R7 | 0.025381 | 0.025381 | 21.108 | 34.521 | 23,314 → 23,373 (+59) | 38,922 → 252,068 (+213,146) |
| R8 | 0.030457 | 0.030457 | 23.662 | 50.016 | 32,868 → 32,898 (+30) | 55,824 → 355,011 (+299,187) |
| R9 | 0.010152 | 0.015228 | 21.591 | 41.022 | 29,288 → 29,327 (+39) | 46,633 → 290,434 (+243,801) |
| R10 | 0.076142 | 0.050761 | 30.231 | 67.374 | 54,087 → 54,143 (+56) | 68,874 → 341,993 (+273,119) |

S2 improved RF in 3 replicates, tied S1 in 5, and was worse in 2. The large
runtime increase is associated with full cross-tree search: the number of
unique X candidates grew by only 0.15%, but the number of scored DP splits grew
by approximately 6.25×.

Outputs and logs are under
`pro_data/benchmark-stelar-pro-s1-vs-s2-gpu-dupmean-0.25/`. S1 runtime and RF
come from the existing strict-GPU S1 benchmark; its candidate counts were
recomputed with the current build.
