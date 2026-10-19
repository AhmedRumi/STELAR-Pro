# S2 without cross-tree transitions: R1 experiment

This controlled variant retains S1 candidates, all DISCO trees, completion,
UPGMA, and every local DP transition. Only cross-tree DP transitions are
disabled. Triplet weights still use the original rooted/tagged multicopy gene
trees. The variant was built from a temporary source copy; the main codebase
was not changed.

RF is standard unrooted RF against the true single-copy species tree. Lower is
better. The 200-taxon comparisons below use the duplication-mean-1 R1 unless
otherwise stated.

| Dataset | Method | RF | RF rate |
|---|---|---:|---:|
| 200 taxa, 1-duplication R1 | S2 local-only (this experiment) | 16/394 | 4.06% |
| 200 taxa, 1-duplication R1 | STELAR-Pro S1, S2, S3, S4, S5, S6 (each) | 16/394 | 4.06% |
| 200 taxa, 1-duplication R1 | ASTRAL-Pro3 | 12/394 | 3.05% |
| 200 taxa, 0.25-duplication R1 | S2 local-only; S1; full S2 (each) | 28/394 | 7.11% |
| 200 taxa, 0.25-duplication R1 | ASTRAL-Pro3 | 24/394 | 6.09% |
| 2001 taxa, 1-duplication/loss-0.1 R1 | S2 local-only (this experiment) | 294/3996 | 7.36% |
| 2001 taxa, 1-duplication/loss-0.1 R1 | STELAR-Pro S1 | 508/3996 | 12.71% |
| 2001 taxa, 1-duplication/loss-0.1 R1 | ASTRAL-Pro3 | 318/3996 | 7.96% |

The large simulation is called “2000 taxa,” but the saved true and inferred
trees each have 2,001 taxa. No separate ASTRAL-Pro2 output or executable was
available; the ASTRAL result above is explicitly Pro3, not a Pro2 measurement.

| R1 dataset | Local-only unique DP splits | Local-only scored splits | Full S2 scored splits |
|---|---:|---:|---:|
| 200 taxa, 1-duplication | 72,649 | 72,649 | 471,101 |
| 2001 taxa, 1-duplication/loss-0.1 | 1,197,111 | 1,005,729 | 8,464,201 |

For the 2001-taxon run, all 60,003 retained DISCO trees were used; completion
added 276,281 unique candidate clusters. Scoring 1,005,729 reachable splits
took 2,515.1 seconds; total run time was 2,627.4 seconds (43m 47s). Peak CPU
RAM was 8,116 MiB and peak GPU VRAM was 704 MiB.

The large run used the existing exact GPU bitset-intersection backend in the
temporary build because the default smaller-side scorer projected multiple
hours. The bitset backend reproduced the default topology and exact triplet
score on the 200-taxon 1-duplication R1 and on the multicopy S25 control.
The default backend was not run to completion on the 2001-taxon input.
