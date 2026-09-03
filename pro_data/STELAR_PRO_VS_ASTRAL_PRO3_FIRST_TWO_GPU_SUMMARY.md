# GPU benchmark and realized duplication summary

The three generated datasets do **not** realize the requested mean duplicate
copy levels of 0.25, 1, and 3 per species per gene tree. Their mean tree sizes
are only 200.2111, 200.3306, and 201.5176 leaves. Thus, the total input sizes
differ by less than 0.7%, and most gene trees remain exactly single-copy. This
explains why their running times are similar.

`Events/tree` is SimPhy's realized `Locus_Trees.n_dup / 1000`.
`Extra copies/species/tree` is `(leaves - 200) / (200 × 1000)` and directly
measures the realized copy excess relevant to the requested 0.25/1/3 levels.

## Dataset-wide realized duplication

| Requested duplicates/species/tree | SimPhy `-lb` | Mean leaves/tree | Extra copies/species/tree | Events/tree | Multicopy gene trees |
|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.8e-10 | 200.2111 | 0.0010555 | 0.0175 | 1.73% |
| 1 | 2e-10 | 200.3306 | 0.0016530 | 0.0390 | 3.80% |
| 3 | 4.1e-10 | 201.5176 | 0.0075880 | 0.0941 | 8.98% |

## Per-replicate gene-tree statistics

Each row summarizes 1,000 gene trees.

| Requested level | Replicate | Mean leaves/tree | Extra copies/species/tree | Events/tree |
|---:|:---:|---:|---:|---:|
| 0.25 | R1 | 200.058 | 0.0002900 | 0.0180 |
| 0.25 | R2 | 200.039 | 0.0001950 | 0.0160 |
| 0.25 | R3 | 200.042 | 0.0002100 | 0.0140 |
| 0.25 | R4 | 200.382 | 0.0019100 | 0.0200 |
| 0.25 | R5 | 200.293 | 0.0014650 | 0.0200 |
| 0.25 | R6 | 200.599 | 0.0029950 | 0.0130 |
| 0.25 | R7 | 200.166 | 0.0008300 | 0.0140 |
| 0.25 | R8 | 200.308 | 0.0015400 | 0.0210 |
| 0.25 | R9 | 200.184 | 0.0009200 | 0.0220 |
| 0.25 | R10 | 200.040 | 0.0002000 | 0.0170 |
| 1 | R1 | 200.245 | 0.0012250 | 0.0340 |
| 1 | R2 | 200.350 | 0.0017500 | 0.0390 |
| 1 | R3 | 200.132 | 0.0006600 | 0.0300 |
| 1 | R4 | 200.046 | 0.0002300 | 0.0160 |
| 1 | R5 | 200.848 | 0.0042400 | 0.0540 |
| 1 | R6 | 200.327 | 0.0016350 | 0.0350 |
| 1 | R7 | 200.283 | 0.0014150 | 0.0420 |
| 1 | R8 | 200.361 | 0.0018050 | 0.0480 |
| 1 | R9 | 200.164 | 0.0008200 | 0.0450 |
| 1 | R10 | 200.550 | 0.0027500 | 0.0470 |
| 3 | R1 | 200.590 | 0.0029500 | 0.0860 |
| 3 | R2 | 206.436 | 0.0321800 | 0.1440 |
| 3 | R3 | 200.796 | 0.0039800 | 0.0980 |
| 3 | R4 | 200.793 | 0.0039650 | 0.0850 |
| 3 | R5 | 200.519 | 0.0025950 | 0.0810 |
| 3 | R6 | 200.773 | 0.0038650 | 0.0740 |
| 3 | R7 | 201.941 | 0.0097050 | 0.1060 |
| 3 | R8 | 201.589 | 0.0079450 | 0.0890 |
| 3 | R9 | 201.398 | 0.0069900 | 0.0950 |
| 3 | R10 | 200.341 | 0.0017050 | 0.0830 |

## Completed GPU benchmark: first two datasets

Values are mean ± sample SD over 10 replicates. Times are full wall-clock
seconds. Both programs used 28 CPU threads; STELAR-Pro additionally used strict
CUDA on the RTX 4060 with S1/I1. RF is standard unrooted species-tree RF divided
by `2(n-3)`; lower is better.

| Requested level | STELAR-Pro time (s) | ASTRAL-Pro3 time (s) | Speedup | STELAR-Pro RF rate | ASTRAL-Pro3 RF rate |
|---:|---:|---:|---:|---:|---:|
| 0.25 | 26.226 ± 3.827 | 67.309 ± 7.778 | 2.57× | 0.042640 ± 0.029811 | 0.024365 ± 0.014318 |
| 1 | 27.749 ± 6.028 | 73.024 ± 4.908 | 2.63× | 0.095431 ± 0.151904 | 0.024365 ± 0.011917 |

## Per-replicate RF rates

Each cell reports the normalized RF rate followed by the raw RF distance out
of the maximum `394` for a 200-taxon binary tree.

| Requested level | Replicate | STELAR-Pro RF rate | ASTRAL-Pro3 RF rate |
|---:|:---:|---:|---:|
| 0.25 | R1 | 0.071066 (28/394) | 0.060914 (24/394) |
| 0.25 | R2 | 0.020305 (8/394) | 0.010152 (4/394) |
| 0.25 | R3 | 0.040609 (16/394) | 0.030457 (12/394) |
| 0.25 | R4 | 0.101523 (40/394) | 0.025381 (10/394) |
| 0.25 | R5 | 0.025381 (10/394) | 0.015228 (6/394) |
| 0.25 | R6 | 0.025381 (10/394) | 0.015228 (6/394) |
| 0.25 | R7 | 0.025381 (10/394) | 0.020305 (8/394) |
| 0.25 | R8 | 0.030457 (12/394) | 0.025381 (10/394) |
| 0.25 | R9 | 0.010152 (4/394) | 0.015228 (6/394) |
| 0.25 | R10 | 0.076142 (30/394) | 0.025381 (10/394) |
| 1 | R1 | 0.040609 (16/394) | 0.030457 (12/394) |
| 1 | R2 | 0.015228 (6/394) | 0.005076 (2/394) |
| 1 | R3 | 0.106599 (42/394) | 0.030457 (12/394) |
| 1 | R4 | 0.512690 (202/394) | 0.045685 (18/394) |
| 1 | R5 | 0.030457 (12/394) | 0.020305 (8/394) |
| 1 | R6 | 0.015228 (6/394) | 0.010152 (4/394) |
| 1 | R7 | 0.131980 (52/394) | 0.025381 (10/394) |
| 1 | R8 | 0.035533 (14/394) | 0.020305 (8/394) |
| 1 | R9 | 0.050761 (20/394) | 0.035533 (14/394) |
| 1 | R10 | 0.015228 (6/394) | 0.020305 (8/394) |

The STELAR-Pro RF variance at level 1 is dominated by R4 (`202/394 =
0.512690`); its process exited normally and its inferred tree passed the
single-copy taxa validation.
