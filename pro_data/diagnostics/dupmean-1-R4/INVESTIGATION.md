# R4 high-RF investigation

## Conclusion

R4's high RF rate is primarily a constrained-S1 candidate problem. The current
STELAR-Pro build uses only the unique speciation-driven subtree-bipartitions
observed in the rooted/tagged gene trees (`CB = UGB`). Several important,
nested true-tree bipartitions are absent, so their descendant combinations
cannot all occur together in any tree admitted by the S1 DP.

This run did not use S2/S3 enrichment. It also did not use the legacy
cross-tree `search-mode full`: the current STELAR-Pro scope guard accepts only
`S1 + local` and rejects `full` before preprocessing.

## R4 pipeline audit

- 200 taxa and 1,000 binary gene trees; no unresolved polytomies remained.
- 16 duplication nodes and 199,030 eligible speciation nodes were found.
- No speciation node was skipped for overlapping child-species sets.
- 199,030 observed speciation-driven bipartitions collapsed to 80,263 unique
  candidate bipartitions.
- The corresponding cluster table contained 64,990 unique clusters.
- The DP verifier passed every structural assertion: 64,791 parent states,
  80,263 unique transitions, 608 root transitions, and no transition side
  outside the cluster table.
- All 199 rooted subtree-bipartitions in the inferred STELAR-Pro tree belong to
  the candidate set, as required.

The standard unrooted RF tree has 197 internal splits. Of those true splits,
189 have at least one orientation represented in the cluster set and eight do
not. The smaller sides of the eight absent true RF splits contain 27, 36, 38,
41, 46, 55, 75, and 84 taxa. Thus, the missing set includes balanced/deep
splits such as `84|116` and `75|125`, not only terminal details.

For the rooted STELAR objective, the more relevant comparison is exact
subtree-bipartitions `(left-child taxa | right-child taxa)`. The true species
tree has 199 such bipartitions: 181 are in `CB`, while 18 are absent. Nine have
a parent cluster that is itself absent; for the other nine, the parent cluster
exists but its true child split was never observed at a speciation node.

| Missing true bipartition | Parent size | Child sizes | Parent cluster in X? |
|---:|---:|---:|:---:|
| 1 | 11 | 5 + 6 | yes |
| 2 | 13 | 5 + 8 | yes |
| 3 | 14 | 6 + 8 | yes |
| 4 | 18 | 7 + 11 | yes |
| 5 | 20 | 8 + 12 | yes |
| 6 | 23 | 10 + 13 | yes |
| 7 | 25 | 9 + 16 | yes |
| 8 | 27 | 11 + 16 | no |
| 9 | 36 | 13 + 23 | no |
| 10 | 38 | 13 + 25 | no |
| 11 | 39 | 3 + 36 | yes |
| 12 | 46 | 7 + 39 | no |
| 13 | 55 | 19 + 36 | no |
| 14 | 75 | 20 + 55 | no |
| 15 | 84 | 38 + 46 | no |
| 16 | 159 | 75 + 84 | no |
| 17 | 186 | 27 + 159 | no |
| 18 | 193 | 7 + 186 | yes |

## Comparison with neighboring replicates

`S1 RF lower bound` is an independent dynamic program that uses the exact
STELAR-Pro S1 transition graph but maximizes shared true RF splits instead of
the triplet objective.

| Replicate | Duplication nodes | Unique candidates | True rooted bipartitions in CB | True RF splits represented in X | S1 RF lower bound | Actual STELAR-Pro RF |
|:---:|---:|---:|---:|---:|---:|---:|
| R3 | 30 | 77,862 | 191/199 | 197/197 | 22/394 = 0.055838 | 42/394 = 0.106599 |
| R4 | 16 | 80,263 | 181/199 | 189/197 | 124/394 = 0.314721 | 202/394 = 0.512690 |
| R5 | 54 | 49,306 | 199/199 | 197/197 | 0/394 = 0 | 12/394 = 0.030457 |

The mean gene-tree RF rates among single-copy trees were 0.724425 (R3),
0.734204 (R4), and 0.572071 (R5). R4 therefore has the highest discordance and
the largest unique candidate set, despite having the fewest duplication nodes.
Its failure is not caused by having too few candidates overall; the observed
signal is highly dispersed and misses specific nested true bipartitions.

## Objective and implementation checks

| Species tree scored by STELAR-Pro | Rooted-triplet score |
|---|---:|
| True SimPhy species tree | 545,940,960 |
| ASTRAL-Pro3 output | 544,192,584 |
| STELAR-Pro S1 output | 534,120,754 |
| Closest-RF tree permitted by the S1 transitions | 499,449,285 |

The true and ASTRAL-Pro3 trees have better STELAR-Pro objective scores, but
they are not feasible in the current candidate set: only 181/199 and 178/199
of their rooted subtree-bipartitions, respectively, belong to `CB`. The
closest-RF feasible tree scores substantially below the inferred tree, so the
S1 DP correctly prefers the inferred tree under its optimization objective.

The GPU inference score for the selected tree (`534,120,754`) was reproduced
exactly by an independent CPU score-only run. This, together with the passing
cluster/DP structural checks, gives no evidence of a hashing, weight-kernel,
duplicate-collapse, or CUDA error in this replicate.

The source still contains an unreachable complement-side `emitType2` helper
and unrooted ASTRAL-X-era comments. Complement rotations are not part of the
rooted constrained STELAR formulation (`CB = UGB`) used here, so this dead code
did not cause R4's result, although its comments should eventually be cleaned
up to avoid confusing it with STELAR-Pro's rooted S1 path.
