# Experiment harness: S4, S6 `--subsetlocal`, ASTRAL-Pro3

Runs STELAR-Pro's `S4` and `S6 --subsetlocal` arms, and ASTRAL-Pro3 as an
external baseline, over a GDL dataset tree,
recording **elapsed wall time**, **CPU time**, **peak host memory**, **peak GPU
memory** and **normalized unrooted RF** for every replicate.

- `run_s4_s6.py` — the driver: discovers replicates, prepares input, runs each
  arm, scores it, writes `results.tsv` plus a per-run JSON record.
- `measure.py` — runs one command and reports its wall time, CPU time and peak
  memory, host and GPU. Usable on its own for any single run.
- `parse_phases.py` — recovers STELAR-Pro's own per-phase instrumentation from
  a run's log. Also runnable directly on any `*.stelar-pro.log`.

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

## Arms

`--arms` is ordered, and the first entry is the baseline every other arm is
compared against in the paired tables. The default is
`S4,S6-subsetlocal`; ASTRAL-Pro3 is opt-in so an existing command's behaviour
does not change:

```bash
python3 experiments/run_s4_s6.py --data-dir /path/to/data \
    --arms S4,S6-subsetlocal,astral-pro3 --threads 32 --xmx 128g
```

| Arm | Program | Defining options |
|---|---|---|
| `S4` | `./run.sh` | `--search-space S4` |
| `S6-subsetlocal` | `./run.sh` | `--search-space S6 --subsetlocal` |
| `S6` | `./run.sh` | `--search-space S6` |
| `astral-pro3` | `ASTER-Linux/bin/astral-pro3` | — |

**ASTRAL-Pro3 reads the same prepared input file** as the STELAR-Pro arms, so
the comparison cannot differ because of its input: copies already sit as
separate leaves labelled with their species name, which is exactly what
ASTRAL-Pro expects, so no `-a` mapping file is needed.

**It is given every available core by default.** `-t` is passed explicitly
because ASTRAL-Pro3 defaults to **1 thread** on its own, which would otherwise
make the arm look far more expensive than it is. It does genuinely parallelize —
measured on 6000 gene trees, 14 taxa:

| `-t` | wall s | CPU s | utilization | speedup |
|---:|---:|---:|---:|---:|
| 1 | 6.02 | 5.90 | 0.98x | 1.00x |
| 2 | 3.64 | 6.78 | 1.86x | 1.65x |
| 4 | 2.29 | 7.67 | 3.35x | 2.62x |
| 8 | 2.20 | 9.96 | 4.53x | 2.74x |

Two things to read off that. Scaling is real but sublinear, flattening past 4
threads on an 8-core machine with 4 performance and 4 efficiency cores. And
**total CPU time rises with thread count** — 5.90 s to 9.96 s — so more threads
buys wall time at the cost of CPU, which is why both are recorded.

`--threads` sets this for every arm, and its default is the number of cores
*actually available*: it shares `measure.py`'s denominator, so
`sched_getaffinity` and a cgroup `cpu.max` quota are both respected rather than
`os.cpu_count()` being taken at face value. A node that pinned the job to 8 of
its 128 CPUs would otherwise be oversubscribed 16-fold while the efficiency
columns divided by the smaller figure. The plan banner names the source
(`affinity`, `cgroup-quota`, or `cpu_count` off Linux) whenever the default is
in use.

Two options do not apply to it and are recorded as such rather than
misreported: `--xmx` (no JVM, so the `xmx` column is empty for its rows) and
`--backend` (no GPU support at all, so the `backend` column records `cpu` and
no GPU figure is collected). The GPU fallback warning is likewise suppressed
for it. It also has no `PhaseLogger`, so `internal_total_ms`, `phase_count`,
the `log_peak_*` columns and the phase table are empty for its rows — the
wall/CPU/memory columns, which are measured from outside, are fully populated.

Its own objective goes in a **separate** `quartet_score` column, never in
STELAR-Pro's `triplet_score`: they are different objectives over different
units and a shared column would invite comparing numbers that are not
comparable.

The binary is built per machine by `ensure_backends.sh` and is not tracked by
git, so a fresh checkout needs `./ensure_backends.sh --cpu-only` once. If it is
missing, that arm fails with a message saying exactly that rather than a bare
`FileNotFoundError`.

> The script is still called `run_s4_s6.py`, which is now historical — it runs
> whatever `--arms` names. Renaming it would break existing commands.

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
| `<arm>.measure.json` | time, CPU and memory record for that run |
| `<arm>.memory-trace.tsv` | every memory sample: `seconds`, `rss_mib`, `gpu_mib`, `gpu_device_mib` |
| `<arm>.stelar-pro.log` | STELAR-Pro's own log, parsed for the phase table |
| `<arm>.stdout.log`, `<arm>.stderr.log` | captured streams |

And one `results.tsv` at the root of the output location, appended to across
invocations, with `rf_rate`, `raw_rf`/`max_rf`, `elapsed_seconds`,
`cpu_seconds`, `cpu_utilization`, `cpu_seconds_per_core`, `cpu_efficiency`,
`cpu_count`, `peak_rss_mib`, `peak_gpu_mib`, `peak_gpu_device_mib`, the
log-derived columns below, and the settings each row was produced under.

A companion **`results-phases.tsv`** holds one row per phase, in long format:
`condition`, `replicate`, `arm`, `phase_index`, `phase`, `device`, `phase_ms`,
`peak_ram_mib`, `peak_ram_so_far_mib`, `peak_vram_mib`, `peak_vram_so_far_mib`.
It is a separate table rather than extra columns because the number of phases
varies and the phase *sets* differ between arms — S4 runs `Phase 1b S4 DISCO +
completion + UPGMA guide` where S6 runs `Phase 1b S6 rooted-triplet greedy
proposals` — so one column per phase name would be mostly empty. Join on the
first three columns. The two GPU columns are left **empty**, not zero, when no
GPU was visible — "this run used no device memory" and "nobody looked" must not
read the same.

## STELAR-Pro's own instrumentation

`PhaseLogger` already measures more finely in time than any external wrapper
can: it polls `/proc/self/status` `VmRSS` every **50 ms**, and during GPU phases
`cudaMemGetInfo()` every **50 ms**, reporting a per-phase and a run-wide peak
alongside each phase's wall time. It prints all of that for a human to read and
keeps nothing, so `parse_phases.py` reads it back out of the log.

Into `results.tsv`: `internal_total_ms`, `triplet_score`, `log_peak_ram_mib`,
`log_peak_vram_mib`, `phase_count`, `slowest_phase`, `slowest_phase_ms`,
`unfinished_phase`. Into `results-phases.tsv`: every phase.

`internal_total_ms` is **not** the same quantity as `elapsed_seconds`, and the
difference is the point of keeping both. It is a `System.nanoTime()` delta from
`Main.java:119` to `Main.java:914`, so it excludes JVM startup and class
loading, the `run.sh` wrapper, the `build.sh`/`ensure_backends.sh` step,
argument parsing, compute-mode probing and the banner — measured at 503 ms
against 781 ms of wall clock on a small test run. Being logged from a `finally`
block, it appears even when the run throws, so a `Total time` line is not
evidence of success; `unfinished_phase` is, naming the phase the log ends inside.

Two properties of these numbers belong to the source, not the parsing:

- **`peak_vram_mib` is device-wide, not per-process.** `PhaseLogger` derives it
  as `total - minFree` from `cudaMemGetInfo`, which reports the whole card, so
  another job sharing the GPU inflates it. `measure.py`'s `peak_gpu_mib` is the
  pid-attributed figure, and this one has the better time resolution; they
  answer different questions and both are kept.
- **`peak_ram_mib` changes meaning off Linux.** When `/proc/self/status` cannot
  be read, `readVmRssMiB()` falls back to `totalMemory() - freeMemory()` — Java
  heap only, not RSS, and much smaller. Fine on Linux; figures from a macOS run
  are not comparable.

## How the harness's own measurements are taken

**RF** comes from `scripts/calculate_rf_rate.py` — normalized unrooted RF
against `species-tree/s_tree.trees`. It is dependency-free, which matters here
because DendroPy is not installed. `rf_rate = raw_rf / max_rf`; on two binary
trees over the same taxa, `raw_rf / 2` is the number of true bipartitions
missed, so the paired table's `splits` column is `(raw_rf_b - raw_rf_a) / 2`.

**Elapsed time** is `time.perf_counter()` around the whole `run.sh` invocation,
so it includes JVM startup. `perf_counter` is monotonic, so an NTP step cannot
make a long run report a negative duration. This is the honest number for a
cost comparison. STELAR-Pro's own `Total time: N ms` log line is *not* used: it
was measured in earlier sweeps to be non-monotonic across arms under varying
background load, which makes it useless for cross-arm cost.

**CPU time** is `ru_utime + ru_stime` — the time the job was actually scheduled
on a CPU, summed over every thread, as the kernel accounts it. Wall time and
CPU time answer different questions: wall time is what you wait, CPU time is
what you spent, and `cpu_utilization` (CPU / wall) is the average number of
cores kept busy. That ratio is the useful one for comparing search spaces, since
an arm that is slower in wall time but far below `--threads` in utilization is
losing to a serial section, not to raw work.

**CPU time per CPU** is `cpu_seconds / cpu_count`, and `cpu_efficiency` is the
same thing as a fraction (`cpu_utilization / cpu_count`, 1.0 = every available
CPU busy for the whole run). Nothing has to be summed to get there:
`ru_utime + ru_stime` is *already* aggregated over every core and thread, so it
is one division. `cpu_seconds_per_core` is the wall time the job would take on
a perfect parallel machine, which makes it directly comparable to
`elapsed_seconds` — the gap between the two is the headroom lost to serial
sections.

Read it as **normalized, not as what any one core did**: 3 CPUs fully busy and
32 CPUs each 9% busy give the same figure.

The denominator matters enough to be recorded. `cpu_count` comes from
`sched_getaffinity`, not `os.cpu_count()`, so a 128-core node that pinned the
job to 8 CPUs does not make it look 6% efficient; a cgroup v2 `cpu.max` quota is
tighter still and wins when set, which is what makes the figure meaningful
inside a container. `cpu_count_source` says which applied (`affinity`,
`cgroup-quota`, or `cpu_count` off Linux). `physical_cores` is recorded
separately, from `/proc/cpuinfo`, because on a 16-core/32-thread box the choice
of denominator moves the per-CPU figure by 2x and which one is wanted depends
on the question — `cpu_count` counts *logical* CPUs.

Two things CPU time does **not** measure. It excludes time spent blocked, so a
job starved of I/O or descheduled behind other load looks cheap — which is
exactly why both numbers are kept. And it excludes GPU kernel execution: under
`--backend gpu`, work moved onto the device leaves wall time unchanged while CPU
time falls, so CPU time understates the cost of a GPU run rather than measuring
it. Pair it with `peak_gpu_mib` there, not with wall time alone.

Both figures come from **`os.wait4` on the job's own pid**, not from
`getrusage(RUSAGE_CHILDREN)`. That distinction is not pedantry: the
`RUSAGE_CHILDREN` bucket accumulates every child the process reaps, which
includes each `ps` and `nvidia-smi` the sampler spawns — of the order of 30,000
processes over a multi-hour run at the default 1 s interval. A *max* like
`ru_maxrss` shrugs that off, since those helpers are tiny beside a JVM; a *sum*
like CPU time does not, and would have billed minutes of sampling overhead to
the job. `wait4` counts the child and the descendants *it* reaped, and nothing
else, so the JVM's time rolls up through `run.sh` while the sampler's does not.

One caveat on a timed-out run: the job is killed by process group, so a shell
that dies without reaping its own JVM loses that JVM's accounting. CPU time for
a row with `status` `timeout` is a lower bound.

**Peak host memory** is taken two independent ways, because neither alone is
trustworthy, and the larger is reported as `peak_rss_bytes`:

- `rusage_peak_bytes` — `getrusage(RUSAGE_CHILDREN).ru_maxrss` after the child
  is reaped. Exact for one process, but it is the largest *single* child, not a
  sum. The unit is platform-dependent and is normalized (bytes on darwin,
  kibibytes on Linux).
- `sampled_peak_bytes` — the maximum over samples of the summed RSS of the whole
  process subtree, so `run.sh`, its JVM and the `astral-pro3` rooting
  subprocess are all counted. Can miss a spike between samples; tighten with
  `--sample-interval`.
- `hwm_sum_bytes` / `hwm_max_bytes` — on Linux, the sum and the max of each
  subtree process's own `VmHWM`, the kernel's per-process peak-RSS high-water
  mark. **Gap-free**: unlike a sampled sum these cannot miss a spike between
  ticks. The sum is an upper bound on the simultaneous total, since the
  per-process peaks need not coincide; the max is exact for the largest
  process. Values are accumulated per pid as they are seen, so a process that
  exits still contributes — which is the point, since the preprocessing
  subprocesses are the short-lived ones.
- `peak_swap_bytes` — on Linux, the peak of the subtree's summed `VmSwap`. This
  closes a hole in RSS: only *resident* pages are RSS, so a job pushed into
  swap reports **less** memory while running slower. The driver warns when it
  is non-zero, because `peak_rss_mib` is then understating real demand.

### Sampling rate

The host sampler reads `/proc` directly where it exists, so no subprocess is
spawned per tick and the default interval is **500 ms**. On a platform without
`/proc` it must spawn `ps`, which costs about **11 ms of CPU per spawn**, so the
default there stays at **1.0 s**. The record says which applied in
`memory_source`.

500 ms rather than something finer on purpose. The sampled series does not have
to catch the peak — `hwm_sum_bytes` and `hwm_max_bytes` come from the kernel's
own `VmHWM` and are gap-free at any interval — so the only job left for the
samples is showing the *shape* of memory use over the run, which half-second
resolution does perfectly well over a run of minutes. Going to 50 ms would
multiply the trace file and the sampler's CPU by ten for no gain in the figures
that get reported.

Measured cost of one `/proc` tick, which is one small read per process:

| processes on the box | per tick | at 2 Hz |
|---:|---:|---:|
| 100 | 1.2 ms | 0.25% of one core |
| 400 | 5.1 ms | 1.0% of one core |
| 1000 | 13.8 ms | 2.8% of one core |

A conservative over-estimate, since real procfs is kernel-generated rather than
read from disk. Pass `--sample-interval` to go finer if a short run needs it, or
coarser on an unusually crowded box.

**`nvidia-smi` keeps its own, slower cadence** (`--gpu-sample-interval`,
default 1.0 s). It is the one expensive probe left — tens to hundreds of
milliseconds per call, and it contends with the CUDA driver — so polling it at
the host rate would be self-defeating. Host memory is sampled 20x more often
than GPU memory by design.

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
along with `sampled_peak_bytes`. Wall time, CPU time and `rusage_peak_bytes`
come from `wait4` and are unaffected. GPU sampling costs two extra `nvidia-smi`
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
| `--arms` | any of `S4`, `S6-subsetlocal`, `S6`, `astral-pro3`; the **first listed is the baseline** every other arm is paired against |
| `--astral-exe` | ASTRAL-Pro3 binary (default `ASTER-Linux/bin/astral-pro3`) |
| `--astral-opts` | extra options for that arm, e.g. `"-u 0"` for a topology-only output tree |
| `--threads`, `--xmx` | worker threads and JVM heap; raise both for large taxon counts |
| `--backend cpu\|gpu\|auto` | default `cpu` |
| `--timeout SEC` | kill a run that overruns; the job is launched in its own process group so the signal reaches the JVM |
| `--fresh` | rerun arms whose output tree already exists (default is to reuse and only re-score) |
| `--no-build` | skip `build.sh`; **only** when `build/` already matches the checked-out branch |
| `--sample-interval` | host memory sampling period; default 500 ms with `/proc`, 1.0 s without; `0` disables sampling and every GPU figure |
| `--gpu-sample-interval` | GPU memory sampling period, default 1.0 s |

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

CPU accounting was verified three ways. A single-threaded burner reported
`0.44 s` CPU against `0.51 s` wall, the gap being process startup. Four burners
launched in parallel by a shell and reaped by `wait` — the same shape as
`run.sh` launching a JVM — reported `1.59 s` CPU against `0.45 s` wall,
`3.49x`, which is what confirms that descendant accounting rolls up at all.
And the measurement that justifies using `wait4`: a 3-second `sleep` watched at
a 10 ms interval, so the sampler spawned 118 helper processes, reported
**`0.01 s`** of CPU — while 300 bare `ps` spawns on their own cost **`3.43 s`**
in `RUSAGE_CHILDREN`. The timeout, nonzero-exit and killed-by-signal paths were
re-exercised after the reaping rewrite, including confirming no process survives
the process-group kill.

The per-CPU figures were checked against a known workload: four parallel burners
on an 8-CPU machine reported `3.56x` utilization and `45% busy`, against the
`50%` that four of eight CPUs predicts, the shortfall being process startup; one
burner reported `11%` against a predicted `12.5%`. Both Linux-only denominators
are unreachable on the machine this was written on, so they take their path as
an argument and were unit-tested against fixtures instead of shipped unrun:
`/proc/cpuinfo` with 4 logical CPUs over 2 SMT cores → `2`, two single-core
sockets → `2`, an ARM-style file with no `core id` → `None`; `cpu.max` of
`400000 100000` → `4.0`, `max 100000` → `None`, unparseable → `None`, missing
→ `None`.

The ASTRAL-Pro3 arm was verified end to end on the same synthetic set: three
arms in one run produced two paired tables, both against the `S4` baseline, and
its rows carried wall time, CPU time, per-CPU figures and peak RSS (5 MiB
against the JVM arms' ~70 MiB) with `xmx` empty, `backend` `cpu`,
`triplet_score` empty and `quartet_score` `54463`. Its output tree scored
identically through the RF path in both ASTRAL's default annotated form and
under `-u 0`, so no output option is forced. The missing-binary path was
exercised with `--astral-exe /nonexistent`, which fails that arm with the build
instruction and leaves the others running.

The `/proc` readers cannot run on the machine this was written on, so — as with
the cgroup and cpuinfo parsers — they take their root path as an argument and
were unit-tested against a fixture tree rather than shipped unrun. A fake
`/proc` with `run.sh` → `java` → `astral (pro) 3` plus one unrelated process
confirmed: the subtree walk returns exactly the three related pids, the summed
RSS matches the page arithmetic, `VmHWM` sums to the three processes' peaks
while the unrelated process's much larger `VmHWM` is excluded, `VmSwap` sums
correctly, a `comm` field containing spaces *and* parentheses parses (the
fields are taken from after the last `)`), and a missing pid or root yields
`None`/`{}` rather than an exception. The GPU-rate decoupling was verified by
counting invocations of a stub `nvidia-smi`: 25 host samples over 4 s produced
4 GPU polls and 8 calls, against the 50 calls an undecoupled sampler would have
made.

The phase parser was checked against real logs from both arms, which recovered
8 phases each with times and RAM peaks, and against fixtures for the cases real
logs here cannot produce: a GPU phase carrying VRAM peaks, ANSI colour escapes
(normally absent when stderr is redirected, but `FORCE_COLOR` or a tty would
leave them), a log that stops inside a phase with no Run Summary — where the
run-wide peaks are recovered from the phase blocks instead and
`unfinished_phase` names the phase — and empty and missing files, which yield
an empty summary rather than an exception. One bug was found this way and
fixed: the first version left the previous phase marked open after its
completion tick, so every finished phase was reported as unfinished.

GPU accounting was verified against a **stub `nvidia-smi`**, since the machine
this was written on has no NVIDIA card and an untested parser is worth nothing.
With the stub reporting 2048 MiB for a pid inside the job's subtree, 4096 MiB
for another user's pid, an `[N/A]` row, and two devices holding 6144 + 128 MiB,
the harness recorded exactly 2048 MiB as the job's own and 6272 MiB device-wide.
With the job's pid absent, it recorded `gpu_peak_attributed: false` and printed
the fallback warning.
