# Reproducibility mirror for simulated runs

Simulated inputs and the original results remain under
`$PHYLOGENY_DATA_DIR/simphy/data`. Every STELAR-Pro and ASTRAL-Pro3 simulated
run also refreshes a small, method-first safety copy under:

```
$PHYLOGENY_DATA_DIR/outputs/gdl-simulation/
  stelar-pro-outputs/<dataset>/R1/<setting>/...
  astral-pro3-outputs/<dataset>/R1/<setting>/...
```

Each mirrored dataset directory also contains the SimPhy `.command` and
`.params` records. Each result leaf contains its inferred tree, statistics,
logs, hidden success/lock markers, and an `<output>.command` record containing
the exact shell command, absolute input/output paths, git revision, host, exit
code, and runtime. Incomplete datasets additionally record their pruning
fraction, seed, and minimum retained taxa.

The mirror deliberately excludes `all_gt.tre`, `s_tree.trees`,
`l_trees.trees`, `g_trees*.trees`, SimPhy databases, ZIPs, and `stat-sim.csv`.
A result leaf containing any forbidden file, symlink, or special file is
refused. Refreshes are copied into a temporary directory and renamed into
place, so a mirror leaf never retains stale files from an older run. Mirror
errors are warnings and do not change the inference exit code.

For the standard data root, the output location is always
`$PHYLOGENY_DATA_DIR/outputs/gdl-simulation`. It can be overridden with
`--simulated-outputs-dir` (also accepted as
`--gdl-simulation-outputs-dir`). A destination inside, equal to, or containing
the input data tree is rejected. Use `--no-outputs-mirror` for a deliberately
unmirrored run.

Normal bulk usage mirrors automatically:

```bash
./run-bulk-simulated.sh --method stelar-pro --num-replicates 5 \
  --opts-list "--search-space S1" --no-notify

./run-bulk-simulated.sh --method astral-pro3 --num-replicates 5 \
  --opts "--thread 16 --seed 42" --no-notify
```

The bulk runner prints its dataset/replicate/setting plan. `--dry-run` only
prints the plan, and `--yes` skips confirmation. Existing results can be
backfilled or refreshed without rerunning inference:

```bash
./sync-simulated-outputs.sh --dry-run
./sync-simulated-outputs.sh
./sync-simulated-outputs.sh --methods "stelar-pro,astral-pro3"
```

The mirror can be published folder-by-folder with the same remote layout:

```bash
./upload-bulk-simulated-outputs.sh --dry-run
./upload-bulk-simulated-outputs.sh --sync --yes
```

The uploader repeats the forbidden-file and reproducibility checks immediately
before upload. By default it writes below `ph/d/simulated/outputs/` in the
configured Hugging Face dataset repository.

The cleanup command `clear-bulk-simulated.sh` removes only `simphy/data`, so
the safety copy under `outputs/gdl-simulation` survives a data cleanup.
