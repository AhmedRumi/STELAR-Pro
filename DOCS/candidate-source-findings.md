# Candidate-source findings

Measurements of what each candidate source actually contributes to the DP search
space `X`. Recorded on branch `s2-guide-tree-candidates`.

## The UPGMA guide tree contributes nothing (2026-11-06)

Measured on all eight 200-taxon duplication conditions in
`gdl-datasets/estimated-cleaned` (1000 genes, replicate 01, S2 defaults):

| condition | new clusters | new splits | redundant |
|---|---:|---:|---:|
| `dup1_ILS25` | 0 / 200 | 0 / 200 | 100.0% |
| `dup3_ILS25` | 0 / 200 | 0 / 200 | 100.0% |
| `dup3_loss0.1` | 0 / 200 | 0 / 200 | 100.0% |
| `dup0.25_loss0.1` | 0 / 200 | 1 / 200 | 99.5% |
| `dup1_loss0.1` | 0 / 200 | 1 / 200 | 99.5% |
| `dup3_loss1` | 0 / 200 | 1 / 200 | 99.5% |
| `dup1_loss1` | 0 / 200 | 4 / 200 | 98.0% |
| `dup0.25_loss1` | 0 / 200 | 6 / 200 | 97.0% |

Against 42k-150k clusters and 63k-219k splits in `X`, that is **0% of clusters and
under 0.01% of splits**. The 37-taxon single-copy set (`all_gt_bs_rooted_37.tre`)
gives 0/36 and 0/36.

**Why.** With 1000 genes, DISCO yields 5k-31k components covering every clade UPGMA
finds. UPGMA assumes ultrametricity (a molecular clock), while the guide matrix is
`dist = 1 - sim_Q`, which is neither ultrametric nor additive — so it recovers only
consensus-level structure. What is missing from `X` is deep clades with short
branches, which no component (mean ~49 of 201 taxa) can see.

**Consequences.**

- The bar for any additional guide-tree method is zero, so additions are strict upside.
- Prefer methods that change the *objective and the transform* — NJ/BIONJ and FastME
  on `-log(sim_Q)`, spectral/Fiedler splits — over other linkages (WPGMA, Ward,
  complete), which share the bias that produces these zeros.
- Redundancy is established **at 1000 genes only**. The 250- and 500-gene subsets of
  the same replicates are the cheap check for whether it is gene-count dependent.

**How to reproduce.** `DPTable.addStreamedTransition` returns whether a split was new,
and `S2CandidateStreamer` attributes clusters and splits per source. Each run logs:

```
S2 guide tree contribution: 0/200 clusters new, 0/200 splits new
  (100.0% of its splits were already in the DP)
```

The line is emitted in Phase 5, before split scoring, so a run can be killed as soon
as it appears — 20-65 s rather than the full 110-4183 s.

## Related

- Reachability: 40-55% of clusters in `X` are unreachable (no usable split pair) under
  local-only transitions. See `search_space_summary.md`.
- `search_space_summary.md:272` states S2 does not add S1 candidates by default, but
  `Main.java` adds the originals and the run log reports them. Unresolved.

## NJ / BIONJ as additional guide trees (2026-11-06)

Added on this branch behind `--extra-guide-trees nj,bionj,nj-linear,bionj-linear`
(`nj`/`bionj` use the additive `-log(sim_Q)` transform, `-linear` variants use
`1-sim` as a control). Same eight conditions, replicate 01, 1000 genes.

### Marginal contribution to X

Guides are emitted in list order, so these are marginal counts with UPGMA first.

| condition | UPGMA | nj | bionj | nj-linear | bionj-linear |
|---|---:|---:|---:|---:|---:|
| `dup0.25_loss1` | 0 | 15 | 1 | 6 | 2 |
| `dup3_loss1` | 0 | 13 | 1 | 7 | 0 |
| `dup1_loss1` | 0 | 13 | 1 | 7 | 0 |
| `dup1_loss0.1` | 0 | 8 | 1 | 3 | 1 |
| `dup3_loss0.1` | 0 | 7 | 0 | 6 | 0 |
| `dup0.25_loss0.1` | 0 | 7 | 0 | 0 | 1 |
| `dup1_ILS25` | 0 | 6 | 0 | 0 | 0 |
| `dup3_ILS25` | 0 | 6 | 3 | 1 | 0 |
| **mean** | **0.0** | **9.4** | 0.9 | 3.8 | 0.5 |

NJ adds new clusters in every condition; UPGMA in none. Also: 0 of 20100 taxon
pairs ever have `sim = 0` at 1000 genes, so the `-log` floor (`1e-6`) never engages
and the transform comparison is clean.

### Clade-set overlap between methods (`--dump-guide-trees`)

Jaccard over the 199 rooted clades of size 2..n-1, mean across the eight conditions:

| pair | mean Jaccard |
|---|---:|
| `nj` ~ `bionj` | 0.969 |
| `nj-linear` ~ `bionj-linear` | 0.973 |
| `nj` ~ `nj-linear` | 0.741 |
| `UPGMA` ~ `nj` | 0.746 |

**The transform, not the algorithm, is the axis of variation.** NJ and BIONJ are
effectively one tree — identical (Jaccard 1.000, union 199) on `dup3_loss0.1`.
Adding `bionj` beside `nj` buys 3.1 clades on average; adding `nj-linear` buys 30.6.

Do not read BIONJ's low marginal count as poor quality: NJ is emitted first and
claims the shared clades. The two methods agree, they are not ranked.

Under high ILS all methods converge (`nj`~`nj-linear` rises to 0.94-0.95, union of
all five drops to 220-223 from 247-253) — and those are exactly the conditions where
NJ adds fewest clusters. When the methods agree there is nothing to add.

**UPGMA is the most distinctive guide tree and the least useful one**: 4-16 clades
held by no other guide, `UPGMA`~`nj` of 0.746 (lower than any NJ-family pair), and
0 new clusters in X everywhere. Distinctiveness among guides does not predict
novelty against X; here it is close to anti-correlated.

### Recommendation

`--extra-guide-trees nj,nj-linear` — two trees spanning ~230 of the ~242 clades all
five reach. Drop both BIONJ variants.

### Accuracy: measured, and the gain does not arrive (2026-11-06)

Two paired A/B sweeps, 20 replicates each (`dup1_ILS25`, `dup3_ILS25`, `dup1_loss1`,
`dup0.25_loss0.1`, reps 01-05, 1000 genes), identical inputs byte-for-byte, S2
defaults, `-t 4`:

| flag | better | worse | tied | identical trees | mean delta RF |
|---|---:|---:|---:|---:|---:|
| `nj,nj-linear` | 2 | 3 | 15 | 14/20 | +0.000758 |
| `nj-corrected,nj` | 3 | 3 | 14 | 14/20 | +0.000505 |

Baseline mean normalized RF 0.0593. Neither configuration moves it. Every single
change is +-1 bipartition (1/198 = 0.005051), except one 3-bipartition regression.

**The DP is not ignoring the clusters.** In the `nj-corrected,nj` sweep the
optimization objective improved in exactly the 6 replicates whose tree changed, and
was bit-identical in the other 14 - the changed set and the improved set are the same
set. The score never decreased, which is the expected sanity check: X only grew.

| replicate | triplet score delta | RF delta | |
|---|---:|---:|---|
| `dup0.25_loss0.1` 03 | +49758 | +0.005051 | worse |
| `dup1_loss1` 01 | +41889 | -0.005051 | better |
| `dup0.25_loss0.1` 04 | +18672 | -0.005051 | better |
| `dup0.25_loss0.1` 01 | +4156 | +0.015152 | worse |
| `dup1_loss1` 04 | +320 | +0.005051 | worse |
| `dup1_loss1` 03 | +42 | -0.005051 | better |

**So the finding is objective/accuracy misalignment, not a search-space shortfall.**
The clusters are reachable, scored, and selected; they yield strictly higher-scoring
trees that are no more accurate. The largest score gain (+49758) made the tree worse;
the smallest (+42) made it better. Pearson r over the 6 changed replicates is -0.203
at n = 6 - i.e. no relationship in either direction.

**All 10 ILS replicates have a score delta of exactly 0** despite 1-33 new clusters
each (`dup3_ILS25` 02: 33 new clusters, score unchanged to the digit). Under ILS the
incumbent optimum is already unbeatable within the enlarged space.

Contribution to X in that sweep, mean over 20 replicates: `nj-corrected` **11.7**
new clusters (max 23), `nj` marginal after it 4.3, UPGMA 0.2. `nj-corrected` is
confirmed the strongest single method, ~3x plain `nj`. It changed the answer anyway
in only 6 of 20 runs, and helped in 3.

Sweep2 and sweep3 agree on 18 of 20 output trees, so the two flag choices are nearly
interchangeable downstream despite the 3x difference in clusters donated.

**Consequence.** Candidate-set enrichment is closed as an accuracy lever: the
strongest available configuration, measured on the strongest instrumentation, moves
mean RF by +0.0005 the wrong way. Pursue instead (a) the 40-55% unreachable clusters,
which is a missing-*splits* problem of a different mechanism and three orders of
magnitude larger, or (b) the objective itself, since the DP is now demonstrably
optimizing something that accuracy does not track here.

### Caveats on the A/B

- `dup3_loss1` and `dup0.25_loss1` - the two conditions with the largest NJ
  contribution - have **zero replicates** in both sweeps. The four conditions were
  picked as the *fastest*, for a runtime estimate, before this became an A/B.
- n = 20, and 14 are ties, so the test resolves nothing below about 1 bipartition.
- The clean per-method novelty measurement (free of emission-order effects) is
  `--dump-clusters` intersected with `--dump-guide-trees` offline; still not done.
- Build note: `nj-corrected` required a rebuild. Flag-off on the new build reproduces
  the baseline tree byte-for-byte at triplet score 635352076, so the sweep3 deltas are
  attributable to the flag, not the binary.
