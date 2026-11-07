# Parked investigations

Things measured, half-measured, or proposed and then set aside, with enough detail to
restart them. Nothing here is running. Ordered by what I would pick up first.

## 1. The objective misranks trees against real competitors (strongest unfinished lead)

**Status: 67 of 160 data points collected, not analysed. Scripts ready.**

Every replicate ships precomputed rival trees, so scoring them under STELAR-Pro's own
objective costs ~3 s each and needs no inference run:

| file | method | known RF (n=20) |
|---|---|---:|
| `1000_gt_apro-cleaned.tre` | ASTRAL-Pro3 | 0.0566 |
| `1000-gt-wqfm-tree.tre` | wQFM-GDL | **0.0460** |
| `1000-cons.greedy.tree` | greedy consensus | 0.1024 |
| (ours) | STELAR-Pro S2 | 0.0593 |

The question is not "is the true tree the optimum" but "does the objective rank trees in
the same order as accuracy". wQFM-GDL beats us on RF, so if its tree scores *below* ours
the objective demonstrably prefers the less accurate tree. Greedy is the control: if
score tracks accuracy it must score lowest.

The one replicate completed before this was parked (`dup0.25_loss0.1` 01):

| tree | score | RF |
|---|---:|---:|
| ours | 813,766,614 | 0.116 |
| true tree | 809,560,828 | 0 |
| ASTRAL-Pro3 | 726,620,028 | ~0.057 |

Score order **ours > truth > apro**; RF order **truth > apro >= ours**. Our tree
outscores ASTRAL-Pro's by 12% while being no more accurate. One replicate only.

**To restart:** `scratchpad/rival/one.sh` (apro, rooted already) and
`rival/rooted_one.sh` + `rival/reroot.py` (wQFM and greedy ship **unrooted**, and the
scorer correctly refuses to root arbitrarily; `reroot.py` places the root at the
bipartition matching the true tree's own root split and reports the match Jaccard).
Partial results in `rival/results.tsv`. ~13 min for all 160 at `-P 3`.

## 2. Oracle test is complete; its comparison partner is not

**Status: done for all 160 replicates, `scratchpad/oracle8/results.tsv`.**

True-tree scores collected for all 8 duplication conditions x 20 replicates (447 s,
zero failures). But the comparison needs an *inferred* score per replicate, and
STELAR-Pro outputs exist for only 20 of the 160. On those 20 the inferred tree
outscores the truth in **19/20**, mean +0.319%. Extending that to all 160 means
running STELAR-Pro ~140 more times (~40 h), which is why item 1 exists as the cheap
substitute.

## 3. Loss-heavy sweep with fresh baselines

The conditions with the largest NJ contribution - `dup3_loss1` and `dup0.25_loss1` -
have **zero replicates** in any sweep so far, because the original four were chosen as
the *fastest* for a runtime estimate. Giving them a baseline plus a treatment arm is
~20 runs at 1,600-2,000 s each, ~6 h at two lanes. Inputs for reps 01-20 of both
conditions are already built in `scratchpad/oracle8/inputs/`.

## 4. Perturbation scale is calibrated on one replicate

`--guide-ensemble` scale 100 comes from the effective-sample-size argument
(`sqrt(den/genes)` ~ 100) and was confirmed best on exactly one replicate, where scale
300 came out *worse than baseline*. A 3-point sweep (30/100/300) across several
replicates would pin the optimum instead of assuming n=1 generalizes.

## 5. Smaller items

- **Multi-rooting cluster donation.** FastME's apparently unique clades turned out to be
  an artifact of which edge the optimized topology is rooted on. Donating clusters from
  several rootings of one topology is nearly free and was never tried.
- **Gene-level bootstrap** instead of SE perturbation: resample the 1000 gene trees with
  replacement and rebuild the matrix per member. Statistically correct, needs no scale
  fudge, costs one matrix rebuild (~0.7-0.9 s) per member.
- **`--dump-clusters` intersected with `--dump-guide-trees`** for per-method novelty free
  of emission-order effects. Still the only clean way to attribute clusters per source.
- **Gene-count dependence.** Every redundancy result is at 1000 genes. The 250- and
  500-gene subsets of the same replicates are the cheap check.
- **`DOCS/search_space_summary.md:272`** claims S2 does not add S1 candidates by default,
  but `Main.java` adds the originals and the run log reports them. Unresolved.
- **The exactly-200 (= n-1) gap** between total clusters in X and splittable clusters,
  noticed and never explained.
- **`scripts/__pycache__/*.pyc` is tracked in git** and shows as modified after every
  run. Wants untracking plus a gitignore entry.
