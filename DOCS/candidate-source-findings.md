# Candidate-source findings

Measurements of what each candidate source actually contributes to the DP search
space `X`. Recorded on branch `s2-guide-tree-candidates`.

## Bottom line (2026-11-07, after the 20-replicate four-way run)

1. **Enrichment is a conditional lever, not a dead one.** All 36 guide trees give mean
   normalized RF **0.057071 against the baseline's 0.059343** - the first configuration
   to improve. It is not significant (sign test p = 0.344, permutation p = 0.309).
2. **It helps only where the true tree outscores the inferred one.** 6 of the 9
   bipartitions recovered come from `dup1_loss1` 02, the single replicate the oracle
   test flagged that way. `--score-species-tree` is therefore a cheap *per-replicate*
   predictor of whether more candidates can help - run it first, enrich where truth wins.
3. **Loss creates the gap; duplication does not.** `dup1_loss1` accounts for the whole
   net gain (-9 bipartitions). Both ILS conditions are inert across all four arms while
   donating the *most* guide clusters (121-349 against 22-204), so DISCO already
   witnesses everything under duplication+ILS.
4. **Cost is negligible:** 35 extra guide trees add **+0.83%** to splits scored.
5. **Two or three guide trees are worse than none.** The `+nj` arms regress (+0.000758,
   +0.000505) - enough to perturb the DP onto a different equal-scoring tree, not enough
   to cover a missing clade. Partial enrichment is the worst of the three options.

Sections below are in the order they were measured; where a later run overtook an
earlier conclusion it is marked inline.

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

## The oracle test: the objective's optimum is not the true tree (2026-11-07)

This supersedes the framing above. Everything in this document measured candidate
*sources*; this measures whether the search space was ever the binding constraint.

`--score-species-tree` scores one supplied tree under the same objective the DP
optimizes. Scoring the true species tree of each of the 20 A/B replicates against the
same 1000 gene trees:

**STELAR-Pro's inferred tree outscores the true species tree in 19 of 20 replicates**,
by a mean of +0.319% (min -0.095%, max +1.675%).

| condition | replicates where inferred > true | mean margin |
|---|---:|---:|
| `dup0.25_loss0.1` | 5/5 | +0.349% |
| `dup1_ILS25` | 5/5 | +0.299% |
| `dup1_loss1` | 4/5 | +0.465% |
| `dup3_ILS25` | 5/5 | +0.162% |

Validity check: scoring an *inferred* tree through the same path reproduces the DP's
reported score exactly (438720037, 438748508, 635352076 - three for three), and the
true tree parses with 201 of 201 taxa, nothing dropped. The two numbers are the same
quantity.

### Why this makes the objective the binding constraint

> **Refined 2026-11-07.** This section's original title was "why this closes the
> candidate-source direction". The 20-replicate run below shows that overstated it: the
> argument holds in the 19 replicates where the inferred tree outscores truth, but the
> 20th was then recovered by enrichment for -6 bipartitions. Read what follows as
> scoped to those 19, and the oracle comparison as a per-replicate test.


The gap by which the inferred tree already exceeds the truth is one to five orders of
magnitude larger than anything enrichment moves:

| replicate | truth gap | enrichment delta | ratio |
|---|---:|---:|---:|
| `dup1_loss1` 01 | 405,956 | +41,889 | 10x |
| `dup0.25_loss0.1` 03 | 4,442,795 | +49,758 | 89x |
| `dup0.25_loss0.1` 04 | 2,480,451 | +18,672 | 133x |
| `dup0.25_loss0.1` 01 | 4,205,786 | +4,156 | 1,012x |
| `dup1_loss1` 04 | 573,696 | +320 | 1,793x |
| `dup1_loss1` 03 | 9,025,709 | +42 | 214,898x |
| the other 14 | 90,919 - 5,212,424 | **0** | - |

The DP is not failing to reach good trees. It is reaching trees *better than the truth*
under this objective, and enrichment nudges it a further 0.001% in the same direction.
In those 19 replicates no candidate source can fix RF, because the optimum being
searched for is not the true tree, and Tier 1 (transforms), Tier 2 (ensembles) and
Tier 3 (reachability) alike would only make the DP better at finding the same wrong
optimum. The exception is the 20th, where the truth outscores the inferred tree and
better search does help - see the four-way run below.

This is the diagnostic that should have been run before any of the implementation
below - it was proposed at the outset, and `search_space_summary.md` already warned
that oracle clade injection "is a diagnostic of search-space coverage, not a guarantee
of improved RF accuracy".

### What it does not say

It is not proof that the objective is wrong in the limit. The score is estimated from
1000 estimated gene trees, and a finite-sample optimum is expected to sit above the
truth; 0.319% is the size of that effect here. The actionable reading is that the
deviation dwarfs the search-space effects, so the live directions are the objective
(weighting, normalization, regularization) and the data, not the candidate set.
`dup1_loss1` 02 is the single replicate where the truth does outscore the inferred
tree, and the only one where better search could help.

## Tier 2: the resampling ensemble (2026-11-07)

Implemented in `ResampledGuideEnsemble` behind `--guide-ensemble method:k[:seed[:scale]]`.
k NJ trees from standard-error-perturbed resamples: `p' = clamp(p + scale*z*se)`,
`se = sqrt(p(1-p)/den)`, one draw per unordered pair per member, deduplicated by exact
clade-set equality.

**The calibrated standard error is far too small.** `den` has a median of 1.022e7
quartets per pair at 1000 genes, so `se` is ~1.3e-4 and 20 draws collapse to 12
distinct topologies (union 1.04x one tree). The binomial form treats those 1e7
quartets as independent observations when they derive from 1000 gene trees: the
effective sample size is the gene count, understating the error by
`sqrt(1e7/1e3)` ~ 100x.

| scale | distinct | union clades | vs one tree | Jaccard | guide clusters | guide splits |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 12/20 | 207 | 1.04x | 0.978 | 12 | 19 |
| 10 | 20/20 | 283 | 1.42x | 0.861 | 25 | 53 |
| **100** | 20/20 | **594** | **2.98x** | 0.565 | 92 | 179 |
| 300 | 20/20 | 1703 | 8.56x | 0.207 | 551 | 956 |

At scale 100 the union holds 178-190 of 199 true rooted clades against 127-149 for the
best single member, at 3-4 clades donated per true clade recovered - so the diversity
is signal, not noise, and the ~242-clade ceiling that five stacked methods could not
pass is gone. Scale 300 on `dup1_loss1` 01 nevertheless produced a **worse** tree than
doing nothing (25 wrong bipartitions against the baseline's 24): past some point the
perturbation manufactures wrong clades and the DP selects them. The statistically
principled scale is also the empirically best one.

On `dup1_loss1` 01 alone the configurations rank, by score: baseline 24 wrong, ensemble
scale 1 23, scale 100 **22**, scale 300 25, `nj-corrected,nj` 23, all-36-cumulative 23 -
the best tree sits at the third-lowest score and the highest score gives a worse tree.
Score and accuracy are not monotonically related even where enrichment helps.

> **Partly reversed 2026-11-07.** I concluded here that "none of it survives the oracle
> test" and that the scale-100 gain was plateau jitter. The 20-replicate run below says
> otherwise: the ensemble is part of the configuration that does beat baseline, and the
> gain concentrates exactly where the oracle test predicts it should. The n=1 reading
> was under-powered in both directions - it could not distinguish jitter from a real
> but condition-specific effect.

### Cost

Extra guide splits are disproportionately expensive *per split*. Against a 1.63 ms
average (37,268 ms / 22,898 splits, measured solo), extra guide splits cost 33-161 ms
each - 20 to 100x. Scale 300 adds 1,055 splits for +35 s on a 37 s baseline. Timings
from any run sharing the machine are useless here: the same configuration measured 155 s
against a concurrent job and would have implied 236x.

**This does not make the feature expensive overall**, because the split count barely
grows: across the 20-replicate run, 35 extra guide trees took the total from 1,086,128
splits to 1,095,128, **+0.83%**. A 100x cost on 0.83% of the work is a rounding error
against machine-load variance. See the four-way run's cost table.

### All 15 implemented specs, and one correction

Running every implemented spec cumulatively (UPGMA + 15 + 20 ensemble members = 36
guide trees) found that **`fastme-corrected` is the largest single Tier 1 contributor**
- 19 clusters on `dup1_ILS25` 03, 14 on `dup1_loss1` 01. The earlier conclusion that
FastME "adds nothing" was drawn from `fastme` and `fastme-linear` only and is wrong for
the corrected transform. The BME logs show why: under the other transforms the NJ seed
is already at a BME local optimum (2-7 NNI moves, 0.00-0.16% shorter), while
NEG_LOG_CORRECTED gives 70 moves and 13.48%. Tier 1 and Tier 2 are largely
complementary - the ensemble's marginal contribution falls only 168 -> 152 clusters when
Tier 1 is emitted first.

## The 20-replicate four-way comparison (2026-11-07)

All four configurations run on the same 20 replicates (4 conditions x 5 reps,
200 taxa, 1000 genes, `estimated-cleaned`), identical inputs byte for byte, so
differences are the configuration alone.

| arm | guides | mean normalized RF | vs baseline |
|---|---|---|---|
| baseline | UPGMA only | 0.059343 | - |
| `+nj,nj-linear` | 3 | 0.060101 | better 2, worse 3, tied 15 |
| `+nj-corrected,nj` | 3 | 0.059848 | better 3, worse 3, tied 14 |
| **all 36** | 15 distance + 20-member ensemble | **0.057071** | **better 7, worse 3, tied 10** |

All-36 is the first arm to improve mean RF: -0.002273, or 9 wrong bipartitions
saved across 20 replicates. It is **not statistically significant** - exact sign
test p = 0.344, exact sign-flip permutation on the mean p = 0.309 - and **6 of
those 9 bipartitions come from one replicate**. Dropping `dup1_loss1` 02 leaves
-3 over 19 (mean -0.000797, sign test p = 0.508).

### Which replicate, and why it matters

`dup1_loss1` 02 is exactly the replicate the oracle test singled out: the one
case in 20 where the **true species tree outscored the inferred tree** (-0.095%).
That is the signature of a genuine search gap - the objective preferred a tree
the search could not reach. Enrichment closed it:

| | score margin over truth | bipartitions missed, of 198 |
|---|---|---|
| baseline | **-0.0950%** | 23 |
| all 36 | **+0.1133%** | 17 |

Six bipartitions recovered, and the margin flips positive. Inferred now outscores
truth in **20/20** replicates, up from 19/20.

In the other 19 the inferred tree already outscored truth, so there was no search
gap to close, and RF moves by +-1 to 2 bipartitions with no pattern. **Candidate
enrichment fixes accuracy exactly where, and only where, the search is genuinely
incomplete.** The oracle test measured that at 1/20 replicates; this run recovered
precisely that one. The two results agree rather than conflict: the objective is
the binding constraint 19 times out of 20, and candidate sources are the binding
constraint the remaining once.

### Loss, not duplication, creates the gap

| condition | baseline | all 36 | delta | net bip | better/worse/tied |
|---|---|---|---|---|---|
| `dup1_loss1` | 0.079798 | 0.070707 | **-0.009091** | **-9** | 4/0/1 |
| `dup0.25_loss0.1` | 0.078788 | 0.077778 | -0.001010 | -1 | 3/2/0 |
| `dup1_ILS25` | 0.038384 | 0.038384 | 0.000000 | 0 | 0/0/5 |
| `dup3_ILS25` | 0.040404 | 0.041414 | +0.001010 | +1 | 0/1/4 |

Every net bipartition the all-36 arm gains overall comes from `dup1_loss1`: -9 there
against -1, 0 and +1 elsewhere, summing to the -9 total. The two `+nj` arms show the
mirror image - `dup0.25_loss0.1` +3 with nothing to offset it.

A counting note for anyone re-deriving these: `rf_nodep.rf` returns `mx = 2(n-3) = 396`
at n = 201, and both trees carry 198 internal bipartitions, so bipartitions missed is
`sd/2` and one missed bipartition is `2/396 = 1/198 = 0.005051` of normalized RF.

`dup1_ILS25` has a score delta of **exactly 0 in all five replicates** while
donating 121-301 guide clusters - the largest donations of any condition. The DP
sees hundreds of new clusters and none of them changes the optimum. Under
duplication with ILS but no loss, DISCO already witnesses everything the distance
methods find; under loss its decomposition is lossy and the guides fill real holes.
This is the same pattern as the original UPGMA zeros, now with a mechanism.

The score never decreased anywhere (11/20 changed, 0 decreases), as it must: the
DP maximizes over a superset. `r(score gain, wrong-bipartition change) = -0.864`
at n = 11, but that is carried by the one outlier - without it, `-0.484` at n = 10,
Spearman `-0.421` on all 11. Directionally right, under-powered.

### Cost: the earlier 2x projection was wrong

| | baseline | all 36 |
|---|---|---|
| splits scored, 20 replicates | 1,086,128 | 1,095,128 (**+0.83%**) |
| weight-table seconds, summed | 16,921 | 15,841 (-6.4%) |

**35 extra guide trees add under 1% to the split count.** Total scoring time came
out *lower* for all-36, which is machine load, not a speedup - the baseline arm was
slower on 9/20 replicates despite scoring fewer splits. The correct conclusion is
that the cost is within noise, and that cross-arm wall-clock comparisons on this
machine are worthless. A mid-run claim that runtimes had "roughly doubled" came
from exactly that mistake and was wrong; the sweep finished in 2 h 36 m. The
per-split figure of 20-100x for guide splits is not contradicted - there are simply
too few of them for it to matter.
