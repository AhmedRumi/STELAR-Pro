# Experiment harness: S4 versus S6 `--subsetlocal`

Runs STELAR-Pro's `S4` and `S6 --subsetlocal` arms over a GDL dataset tree,
recording **elapsed wall time**, **peak host memory**, **peak GPU memory** and
**normalized unrooted RF** for every replicate.

- `run_s4_s6.py` — the driver: discovers replicates, prepares input, runs each
  arm, scores it, writes `results.tsv` plus a per-run JSON record.
- `measure.py` — runs one command and reports its time and peak memory, host
  and GPU. Usable on its own for any single run.

## Quick start

Always dry-run first. It prints the full plan and every command it would run,
and writes nothing:

```bash
python3 experiments/run_s4_s6.py --data-dir /path/to/data --dry-run
```

Then the real thing:

```bash
python3 experiments/run_s4_s6.py \
    --data-dir /path/to/data \
    --replicates R1-R5 \
    --threads 8 --xmx 32g
```

Outputs land in each replicate's own `stelar-pro/` directory. To keep the
dataset read-only, mirror the `<condition>/<replicate>/` structure elsewhere:

```bash
python3 experiments/run_s4_s6.py --data-dir /path/to/data \
    --out-dir ~/results/s4-vs-s6 --threads 8 --xmx 32g
```

## Expected dataset layout

```
<data-dir>/taxa2000_gt1000_dup1_loss0.1_pop45800000/
    R1/
        estimated-genetrees/all_gt.trees     <- inferred from, by default
        true-genetrees/all_gt.trees          <- --genetrees true
        species-tree/s_tree.trees            <- RF reference
        stelar-pro/                          <- default output location
    R2/ ...
```

Conditions and replicates are discovered, so several conditions under one
`--data-dir` all run; narrow with `--conditions` and `--replicates`.

## What each run produces

Per replicate directory:

| File | Contents |
|---|---|
| `input-<set>-species-labelled.trees` | gene trees with copy labels normalized to species |
| `input-<set>-species-labelled.json` | preparation metadata, with SHA-256 input IDs |
| `<arm>.tre` | the inferred species tree |
| `<arm>.measure.json` | time and memory record for that run |
| `<arm>.memory-trace.tsv` | every memory sample: `seconds`, `rss_mib`, `gpu_mib`, `gpu_device_mib` |
| `<arm>.stelar-pro.log` | STELAR-Pro's own log |
| `<arm>.stdout.log`, `<arm>.stderr.log` | captured streams |

And one `results.tsv` at the root of the output location, appended to across
invocations, with `rf_rate`, `raw_rf`/`max_rf`, `elapsed_seconds`,
`peak_rss_mib`, `peak_gpu_mib`, `peak_gpu_device_mib`, and the settings each row
was produced under. The two GPU columns are left **empty**, not zero, when no
GPU was visible — "this run used no device memory" and "nobody looked" must not
read the same.

## How the measurements are taken

**RF** comes from `scripts/calculate_rf_rate.py` — normalized unrooted RF
against `species-tree/s_tree.trees`. It is dependency-free, which matters here
because DendroPy is not installed. `rf_rate = raw_rf / max_rf`; on two binary
trees over the same taxa, `raw_rf / 2` is the number of true bipartitions
missed, so the paired table's `splits` column is `(raw_rf_b - raw_rf_a) / 2`.

**Elapsed time** is `time.perf_counter()` around the whole `run.sh` invocation,
so it includes JVM startup. This is the honest number for a cost comparison.
STELAR-Pro's own `Total time: N ms` log line is *not* used: it was measured in
earlier sweeps to be non-monotonic across arms under varying background load,
which makes it useless for cross-arm cost.

**Peak host memory** is taken two independent ways, because neither alone is
trustworthy, and the larger is reported as `peak_rss_bytes`:

- `rusage_peak_bytes` — `getrusage(RUSAGE_CHILDREN).ru_maxrss` after the child
  is reaped. Exact for one process, but it is the largest *single* child, not a
  sum. The unit is platform-dependent and is normalized (bytes on darwin,
  kibibytes on Linux).
- `sampled_peak_bytes` — the maximum over samples of the summed RSS of the whole
  process subtree, so `run.sh` plus its JVM are both counted. Can miss a spike
  between samples; tighten with `--sample-interval`.

**Peak GPU memory** is a *separate* measurement, not a component of the host
one. Device allocations never appear in RSS, so a CUDA run that fills 20 GiB of
VRAM is indistinguishable in `peak_rss_bytes` from one that never touches the
card. It is polled in the same sampling loop via `nvidia-smi`:

- `gpu_peak_bytes` — peak of the memory attributed to **this job's own process
  subtree**, from `nvidia-smi --query-compute-apps=pid,used_gpu_memory`. The
  JVM holds the CUDA context and is a descendant of `run.sh`, so it is matched
  by pid; another user's job on the same card is excluded.
- `gpu_device_peak_bytes` — peak of the total in use across every visible
  device, all users included. Quote `gpu_peak_bytes` as the job's cost, and
  compare it against this one to tell whether the card was shared while
  measuring.

`gpu_visible` records whether `nvidia-smi` ever answered, and
`gpu_peak_attributed` whether any pid of ours was ever found on the card. A
`--backend gpu` run where a GPU was visible but **nothing** was attributable to
it is the signature of a silent fallback from CUDA to a CPU backend, and the
driver prints a warning pointing at the STELAR-Pro log. The same pattern also
appears under a container whose pid namespace differs from the host's, where
`nvidia-smi` reports host pids that cannot be matched; there,
`peak_gpu_device_mib` is the only usable figure.

Both GPU figures depend on sampling, so `--sample-interval 0` disables them
along with `sampled_peak_bytes`. GPU sampling costs two extra `nvidia-smi`
calls per tick; at the default 1 s that is negligible beside a multi-hour
inference, but raise the interval if a run is very long and the shape of the
curve does not matter.

**The trace** (`<arm>.memory-trace.tsv`) keeps every sample rather than only the
maximum, so the shape of memory use over a run can be inspected — where the
peak falls, whether the JVM plateaued at its heap ceiling, whether GPU memory
was held for the whole run or only during one phase. One row per tick, four
columns, about 1 MB for a three-hour run at 1 s.

> **Why not `run-stelar-pro-with-monitor.sh`?** That wrapper calls
> `/usr/bin/time -v`, which is GNU-only. On macOS `/usr/bin/time` is the BSD
> build: it rejects `-v` outright (`illegal option -- v`), and spells the field
> "maximum resident set size" only under `-l`, in **bytes** where GNU reports
> **kibibytes**. So on darwin it yields no memory figure at all, and would
> mis-scale by 1024 if it ever matched. `measure.py` avoids both problems and is
> portable.

Each `measure.py` call is a fresh process, so its `RUSAGE_CHILDREN` covers
exactly one job. Calling `getrusage` repeatedly inside the long-lived driver
would instead return a high-water mark shared by every job it had ever run —
which is why the measurement is isolated in a subprocess at all.

## Input preparation

SimPhy labels gene copies `species_locus_copy` (`45_0_0`), while the species
tree uses the bare species name (`45`). `scripts/prepare_gdl_input.py` rewrites
the labels, leaving copies as separate leaves — which is how STELAR-Pro and
ASTRAL-Pro express multiple copies. Source data is never modified.

This runs **once per replicate** and both arms read the same prepared file, so
the arms cannot differ because of their input. The file is reused on later
invocations unless the raw gene trees are newer or `--fresh` is passed.

If your pipeline already produced `stelar-pro/all_gt_stripped.trees` and
`gene_species_map.txt`, those serve the same purpose; this harness regenerates
them in-repo so a run is reproducible from the dataset alone.

## Options worth knowing

| Option | Why |
|---|---|
| `--dry-run` | print the plan and every command, write nothing |
| `--genetrees estimated\|true` | which gene-tree set to infer from (default `estimated`) |
| `--arms` | any of `S4`, `S6-subsetlocal`, `S6`; a bare `S6` arm is defined for a later three-way |
| `--threads`, `--xmx` | worker threads and JVM heap; raise both for large taxon counts |
| `--backend cpu\|gpu\|auto` | default `cpu` |
| `--timeout SEC` | kill a run that overruns; the job is launched in its own process group so the signal reaches the JVM |
| `--fresh` | rerun arms whose output tree already exists (default is to reuse and only re-score) |
| `--no-build` | skip `build.sh`; **only** when `build/` already matches the checked-out branch |
| `--sample-interval` | memory sampling period, default 1.0 s; `0` disables sampling and every GPU figure |

`--no-build` is off by default on purpose. A stale `build/` silently runs the
wrong code: the tree in `build/` at the time this harness was written predated
the `stelar-pro-new` branch by about 15 hours, which would have measured the
previous commit's S6 while reporting it as the new one. Run `./build.sh` once
after switching branches, then `--no-build` is safe and saves a rebuild per run.

## Verification

The harness was validated end to end on a synthetic dataset in the real
directory layout, with gene trees drawn from a known species tree by random NNI
perturbation. Both arms recovered the true tree exactly on one replicate
(`RF 0.000000`, 0/22) and missed one bipartition on the other (`RF 0.090909`,
2/22) — a mis-wired RF comparison could not produce an exact zero. Time, peak
memory, the TSV, the reuse path, and the replicate and arm selectors were all
exercised. Accuracy figures from that synthetic set mean nothing; only the
plumbing was under test.

GPU accounting was verified against a **stub `nvidia-smi`**, since the machine
this was written on has no NVIDIA card and an untested parser is worth nothing.
With the stub reporting 2048 MiB for a pid inside the job's subtree, 4096 MiB
for another user's pid, an `[N/A]` row, and two devices holding 6144 + 128 MiB,
the harness recorded exactly 2048 MiB as the job's own and 6272 MiB device-wide.
With the job's pid absent, it recorded `gpu_peak_attributed: false` and printed
the fallback warning.
