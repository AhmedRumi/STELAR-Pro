# SimPhy 200-taxon multicopy dataset summary

Generated on 2026-09-01 with `simphy/run_simulator.sh`.

## Common simulation settings

- 200 species, 1,000 gene trees per replicate, and 10 replicates per dataset
- duplication rates (`-lb`): `0.8e-10`, `2e-10`, and `4.1e-10`
- loss rate (`-ld`): `0`
- speciation rate (`-sb`): `0.000001`
- effective population size: uniform from 500,000 to 1,500,000
- SimPhy seed: 42

The requested labels 0.25, 1, and 3 denote nominal mean extra copies per
gene tree. The tables also report the realized values; stochastic ancestral
duplications can affect many descendant species, so realized copy counts vary
substantially among replicates.

`Duplication events` is SimPhy's `Locus_Trees.n_dup`. `Extra copies` is
`n_leaves - 200`; because loss is zero, every gene tree contains all 200
species. `Multicopy trees` counts gene trees with at least one extra leaf.

## Aggregate statistics

| Nominal level | SimPhy `-lb` | Gene trees | Duplication events | Events/tree | Multicopy trees | Multicopy % | Extra copies | Extra copies/tree | Extra copies/species/tree | Mean leaves/tree | Max leaves | Loss events |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.8e-10 | 10,000 | 175 | 0.0175 | 173 | 1.73% | 2,111 | 0.2111 | 0.0010555 | 200.2111 | 396 | 0 |
| 1 | 2e-10 | 10,000 | 390 | 0.0390 | 380 | 3.80% | 3,306 | 0.3306 | 0.0016530 | 200.3306 | 374 | 0 |
| 3 | 4.1e-10 | 10,000 | 941 | 0.0941 | 898 | 8.98% | 15,176 | 1.5176 | 0.0075880 | 201.5176 | 399 | 0 |

Across the 10 replicates, the mean ± sample standard deviation of total
duplication events was `17.5 ± 3.21`, `39.0 ± 10.88`, and `94.1 ± 19.78` for
the three rates. Realized extra copies per gene tree were `0.2111 ± 0.1854`,
`0.3306 ± 0.2297`, and `1.5176 ± 1.8022`, respectively.

## Per-replicate statistics

| Nominal level | Replicate | Duplication events | Multicopy trees | Extra copies | Extra copies/tree | Mean leaves/tree | Max leaves | Loss events |
|---:|:---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | R1 | 18 | 18 | 58 | 0.058 | 200.058 | 211 | 0 |
| 0.25 | R2 | 16 | 16 | 39 | 0.039 | 200.039 | 206 | 0 |
| 0.25 | R3 | 14 | 14 | 42 | 0.042 | 200.042 | 214 | 0 |
| 0.25 | R4 | 20 | 20 | 382 | 0.382 | 200.382 | 396 | 0 |
| 0.25 | R5 | 20 | 19 | 293 | 0.293 | 200.293 | 367 | 0 |
| 0.25 | R6 | 13 | 13 | 599 | 0.599 | 200.599 | 360 | 0 |
| 0.25 | R7 | 14 | 14 | 166 | 0.166 | 200.166 | 343 | 0 |
| 0.25 | R8 | 21 | 21 | 308 | 0.308 | 200.308 | 320 | 0 |
| 0.25 | R9 | 22 | 21 | 184 | 0.184 | 200.184 | 289 | 0 |
| 0.25 | R10 | 17 | 17 | 40 | 0.040 | 200.040 | 209 | 0 |
| 1 | R1 | 34 | 31 | 245 | 0.245 | 200.245 | 312 | 0 |
| 1 | R2 | 39 | 38 | 350 | 0.350 | 200.350 | 272 | 0 |
| 1 | R3 | 30 | 29 | 132 | 0.132 | 200.132 | 222 | 0 |
| 1 | R4 | 16 | 16 | 46 | 0.046 | 200.046 | 220 | 0 |
| 1 | R5 | 54 | 53 | 848 | 0.848 | 200.848 | 342 | 0 |
| 1 | R6 | 35 | 35 | 327 | 0.327 | 200.327 | 328 | 0 |
| 1 | R7 | 42 | 41 | 283 | 0.283 | 200.283 | 358 | 0 |
| 1 | R8 | 48 | 48 | 361 | 0.361 | 200.361 | 366 | 0 |
| 1 | R9 | 45 | 42 | 164 | 0.164 | 200.164 | 224 | 0 |
| 1 | R10 | 47 | 47 | 550 | 0.550 | 200.550 | 374 | 0 |
| 3 | R1 | 86 | 78 | 590 | 0.590 | 200.590 | 359 | 0 |
| 3 | R2 | 144 | 134 | 6,436 | 6.436 | 206.436 | 311 | 0 |
| 3 | R3 | 98 | 97 | 796 | 0.796 | 200.796 | 399 | 0 |
| 3 | R4 | 85 | 82 | 793 | 0.793 | 200.793 | 303 | 0 |
| 3 | R5 | 81 | 76 | 519 | 0.519 | 200.519 | 337 | 0 |
| 3 | R6 | 74 | 71 | 773 | 0.773 | 200.773 | 392 | 0 |
| 3 | R7 | 106 | 102 | 1,941 | 1.941 | 201.941 | 310 | 0 |
| 3 | R8 | 89 | 85 | 1,589 | 1.589 | 201.589 | 344 | 0 |
| 3 | R9 | 95 | 93 | 1,398 | 1.398 | 201.398 | 392 | 0 |
| 3 | R10 | 83 | 80 | 341 | 0.341 | 200.341 | 236 | 0 |

## Output and validation

- `simphy-200taxa-1000gt-10rep-dupmean-0.25-lb-0.8e-10-noloss/`
- `simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss/`
- `simphy-200taxa-1000gt-10rep-dupmean-3-lb-4.1e-10-noloss/`

Each directory contains `R1` through `R10`; every replicate has a 1,000-line
`all_gt.tre` and its `s_tree.trees`. All 30,000 gene trees were checked against
the SimPhy SQLite records: leaf counts match, every tree contains exactly 200
distinct species, and all copy suffixes were collapsed so paralogs appear as
repeated species labels required by STELAR-Pro.
