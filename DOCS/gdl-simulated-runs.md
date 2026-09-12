# Actual GDL simulated runs

`run-bulk-gdl-simulated.sh` is a new, inference-only runner. It does not change
the existing `run-bulk-simulated.sh` SimPhy generation workflow.

Required input layout:

```text
$PHYLOGENY_DATA_DIR/gdl-simulation/data/
  taxa10_gt10_dup1_loss0_pop50000000/
    params.txt                         optional generator provenance
    simphy_raw/simphy_raw.command       optional generator command
    simphy_raw/simphy_raw.params        optional generator parameters
    R1/
      species-tree/s_tree.trees         required true species tree
      true-genetrees/all_gt.trees       required combined gene trees
    R2/...
```

`l_trees.trees` and `true-genetrees/separate/` are not used or required.
The runner never downloads or regenerates datasets.

## Sweep selection

Edit `TAXA_LIST`, `GT_LIST`, `DUP_LIST`, `LOSS_LIST`, `POP_LIST`, and
`NUM_REPLICATES` at the top of the bulk script, or override them:

```bash
./run-bulk-gdl-simulated.sh \
  --taxa-list "10,100,1000" --gt-list "10,1000" \
  --dup-list "1,2" --loss-list "0,0.5" --pop-list "50000000" \
  --num-replicates 5 --method stelar-pro \
  --opts-list "--search-space S1;--search-space S2" --no-notify
```

The Cartesian product is matched literally to directory names. For example,
`loss0.5` and `loss5e-1` are different directory spellings. `dup` is the
directory's duplication parameter (the sample records `DUP_COPIES=1`), not an
assumed duplication rate. Replicates are R1 through the requested RN.

Missing/empty combined inputs or references produce warnings and are skipped.
No matching inputs, any inference failure, or any mirror failure returns a
nonzero exit. `--strict-missing` aborts the entire plan before running methods
if any selected replicate is missing. `--dry-run` validates the file/option
plan without writes or inference; it does not inspect every Newick tree.
Interactive execution asks for confirmation; `--yes` bypasses it.

For the small dataset on the local machine:

```bash
./run-bulk-gdl-simulated.sh \
  --gdl-data-dir /home/aaniksahaa/phylogeny/data/gdl-simulation/data \
  --method stelar-pro --num-replicates 1 --no-notify --dry-run
# Remove --dry-run to execute.
```

For the 10k-species dataset from the example:

```bash
./run-bulk-gdl-simulated.sh \
  --taxa-list 10000 --gt-list 1000 --dup-list 1 \
  --loss-list 0 --pop-list 50000000 --num-replicates 5 \
  --method stelar-pro --opts-list "--search-space S1" --yes
```

## Methods and threads

```bash
./run-bulk-gdl-simulated.sh --method astral-pro3 \
  --num-replicates 1 --opts "--thread 16" --no-notify

./run-bulk-gdl-simulated.sh --methods "stelar-pro,astral-pro3" \
  --stelar-pro-opts-list "--search-space S1;--search-space S2" \
  --astral-pro3-opts-list "--thread 16" --num-replicates 5 --no-notify
```

STELAR-Pro keeps its existing compute/thread defaults. ASTRAL-Pro3 keeps its
native default (one thread); use `--thread N` in its options for multiple CPU
threads. GPU monitoring is always disabled for this CPU-only comparison method,
so unrelated GPU activity is not reported as ASTRAL GPU usage.
`--astral-pro3-bin FILE` overrides the bundled executable. Shared `--opts` or
`--opts-list` applies to every selected method; use method-specific option lists
when their option vocabularies differ. STELAR's `--search-space` is rejected for
ASTRAL before inference. Options overriding input, output, mapping, taxon
restriction, or switching to scoring/tag-only modes are rejected.

For a single replicate without the sweep:

```bash
./test-gdl-simulated.sh --dataset taxa10_gt10_dup1_loss0_pop50000000 \
  --replicate R1 --method astral-pro3 --opts "--thread 8" --no-notify
```

## Species identity and RF

SimPhy GDL leaves such as `10_0_0` represent gene copies of species `10`.
Both methods receive the same deterministic species-labelled normalization:
exact reference species labels remain unchanged; otherwise only numeric
`species_locus_copy` labels whose prefix exists in the reference are accepted.
Distinct gene leaves are retained even when their species labels become equal.
Quoted leaf labels, branch lengths, comments, and internal labels are preserved.
Multiline gene trees are folded to one line each for STELAR's reader.
Unknown/ambiguous labels fail instead of being guessed. The original combined
file and reference are never edited.

Normalization streams trees and verifies their count against the directory's
`gt` parameter; partial datasets fail. Temporary normalized input is removed
after the run and never enters the outputs mirror. Set `TMPDIR` to a disk with
enough room for the normalized file when using large datasets.

The actual reference species set controls mapping and RF, not nominal `taxa`.
The sample `taxa10` dataset has 11 species including species `0`; they are all
retained. Unrooted RF requires identical inferred/reference species sets and
unique species leaves. If loss removes species from all gene trees, or an
inferred tree otherwise fails RF validation, `rf-rate=NA`, `rf-status=error`,
and `rf.log` explain the problem. No silent taxon pruning occurs. Inference
success and RF availability are recorded separately.

## Results, reproducibility and safety

Results are placed at:

```text
gdl-simulation/data/<dataset>/<R>/<method>-outputs/<setting>/
  out-<method>.tre
  out-<method>_stats.csv
  out-<method>.command
  stat-<method>.csv
  input-preparation.json
  input-preparation.log
  replay-gdl.sh
  rf.log
  .<method>_run.log
  .<method>.success
  .<method>.lock
```

`input-preparation.json` records original/reference/normalized SHA-256 hashes,
tree/leaf/species counts and the mapping policy. Command records retain exact
executed commands, temporary execution paths, final result location, source
paths, original invocation and git revision. `bash replay-gdl.sh` recreates
normalization and reruns safely, with notifications disabled. Stats sidecars
refer to persistent source/result paths after staging is finalized.

Completed settings are skipped only when output, success marker and stats exist
and original/reference fingerprints still match. Changed inputs automatically
rerun. `--fresh` forces a rerun. A nonblocking per-setting lock rejects concurrent
duplicate jobs. Inference is staged so failed reruns never replace current
successful output. Old successful runs move to `<setting>__previous_<stamp>`;
failed attempts move to `<setting>__failed_<stamp>` without success/lock markers.
Those leaves are mirrored too, but excluded from combined stats.

Both methods automatically copy results into:

```text
$PHYLOGENY_DATA_DIR/outputs/gdl-simulation/
  <method>-outputs/<dataset>/<R>/<setting>/...
```

The dataset's `params.txt` and `simphy_raw.command`/`simphy_raw.params` are copied
when present (the latter two are flattened into the mirrored dataset root).
Raw gene trees, true species/locus trees, separate inputs, normalized input,
databases and archives are excluded. A mirror failure is fatal for the GDL job
and explicitly reports the failure; completed source results remain available.
Override with `--simulated-outputs-dir`, or explicitly disable with
`--no-outputs-mirror`. Notifications otherwise follow the existing multiline
experiment summary and include dataset, replicate, dup/loss/pop, method options,
RF status, score, time, memory, CSV header/row, and stats location.

Existing collection, backfill and upload tools also support GDL data:

```bash
./collect-stats-simulated.sh \
  --gdl-data-dir "$PHYLOGENY_DATA_DIR/gdl-simulation/data" --out perf-gdl.csv
./sync-simulated-outputs.sh \
  --gdl-data-dir "$PHYLOGENY_DATA_DIR/gdl-simulation/data" --dry-run
./upload-bulk-simulated-outputs.sh \
  --gdl-data-dir "$PHYLOGENY_DATA_DIR/gdl-simulation/data" --sync --dry-run
```

Combined GDL CSV retains `dup`, `loss`, `pop`, `dataset`, `exit-code`, and
`rf-status`. Upload discovery accepts both original SimPhy and GDL dataset
names. It still requires simulation command provenance unless
`--allow-missing-command` is explicitly selected. The remote root remains
`ph/d/gdl-simulation/outputs`. These upload commands are previews only.
