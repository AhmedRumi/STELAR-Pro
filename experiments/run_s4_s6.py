#!/usr/bin/env python3
"""Run STELAR-Pro's S4 and S6 (--subsetlocal) arms over a GDL dataset tree.

For every selected replicate this driver:

  1. normalizes copy labels to species labels with ``scripts/prepare_gdl_input.py``
     (``45_0_0`` -> ``45``; copies stay separate leaves, source data untouched),
     once per replicate and shared by both arms so the arms see identical input;
  2. runs each arm through ``experiments/measure.py``, recording elapsed wall
     time, CPU time, peak host memory and -- on a GPU run -- peak device memory;
  3. scores the result against ``species-tree/s_tree.trees`` with
     ``scripts/calculate_rf_rate.py`` (normalized unrooted RF);
  4. appends one row to ``results.tsv`` and writes a per-job JSON record.

Expected dataset layout, as produced by the SimPhy GDL simulations::

    <data-dir>/taxa2000_gt1000_dup1_loss0.1_pop45800000/
        R1/
            estimated-genetrees/all_gt.trees
            true-genetrees/all_gt.trees
            species-tree/s_tree.trees
            stelar-pro/                 <- default output location
        R2/ ...

Outputs go into each replicate's own ``stelar-pro/`` directory by default. Pass
``--out-dir DIR`` to mirror the ``<condition>/<replicate>/`` structure under a
separate root instead, leaving the dataset read-only.

The arms are deliberately a flag, not hard-coded, so the same harness serves
later comparisons; ``--arms`` takes any of the keys in ``ARMS`` below.

Examples
--------
Dry run first -- it prints the full plan and every command without writing::

    python3 experiments/run_s4_s6.py --data-dir /path/to/data --dry-run

Both arms, five replicates, outputs beside the data::

    python3 experiments/run_s4_s6.py --data-dir /path/to/data \
        --replicates R1-R5 --threads 8 --xmx 32g

Same, but keeping the dataset pristine::

    python3 experiments/run_s4_s6.py --data-dir /path/to/data \
        --out-dir ~/results/s4-vs-s6 --threads 8 --xmx 32g
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
MEASURE = Path(__file__).resolve().parent / "measure.py"

sys.path.insert(0, str(SCRIPTS))
from calculate_rf_rate import compare_newicks, read_tree  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from parse_phases import parse_log  # noqa: E402

@dataclass(frozen=True)
class Arm:
    """How to run one arm.

    ``kind`` picks the command builder, because the arms are no longer all the
    same program: ``stelar`` goes through ``./run.sh``, ``astral`` invokes the
    ASTRAL-Pro3 binary directly. ``options`` is only ever what distinguishes
    this arm from its siblings -- threads, heap and backend are added per run.
    """
    kind: str
    options: tuple[str, ...] = ()


ARMS: dict[str, Arm] = {
    "S4": Arm("stelar", ("--search-space", "S4")),
    "S6-subsetlocal": Arm("stelar", ("--search-space", "S6", "--subsetlocal")),
    # Available for a later three-way comparison; not run unless asked for.
    "S6": Arm("stelar", ("--search-space", "S6")),
    # ASTRAL-Pro3 reads the same prepared input, so the comparison cannot differ
    # because of its input. It has no GPU support, so --backend does not reach
    # it, and no JVM, so --xmx does not either.
    "astral-pro3": Arm("astral"),
}

DEFAULT_ARMS = ["S4", "S6-subsetlocal"]
ASTRAL_DEFAULT = ROOT / "ASTER-Linux" / "bin" / "astral-pro3"

GENETREE_SOURCES = {
    "estimated": Path("estimated-genetrees") / "all_gt.trees",
    "true": Path("true-genetrees") / "all_gt.trees",
}

TRUE_SPECIES_TREE = Path("species-tree") / "s_tree.trees"

RESULT_FIELDS = [
    "condition", "replicate", "arm", "status",
    "rf_rate", "raw_rf", "max_rf", "taxa",
    "output_splits", "true_splits", "shared_splits",
    "elapsed_seconds", "cpu_seconds", "cpu_utilization",
    "cpu_seconds_per_core", "cpu_efficiency", "cpu_count",
    "peak_rss_mib", "peak_gpu_mib", "peak_gpu_device_mib",
    "genetrees", "threads", "xmx", "backend", "exit_code",
    "internal_total_ms", "triplet_score", "quartet_score",
    "log_peak_ram_mib", "log_peak_vram_mib",
    "phase_count", "slowest_phase", "slowest_phase_ms", "unfinished_phase",
    "output_tree", "measurement_json", "finished_at",
]

# One row per phase, in long format. A variable number of phases cannot go into
# results.tsv without either exploding its width or inventing a column per phase
# name, and the phase sets differ between arms (S4 has "Phase 1b S4 DISCO +
# completion + UPGMA guide" where S6 has "Phase 1b S6 rooted-triplet greedy
# proposals"), so they get their own table to join on the first four columns.
PHASE_FIELDS = [
    "condition", "replicate", "arm",
    "phase_index", "phase", "device", "phase_ms",
    "peak_ram_mib", "peak_ram_so_far_mib",
    "peak_vram_mib", "peak_vram_so_far_mib",
]


# ── dataset discovery ────────────────────────────────────────────────────────

def parse_replicates(spec: str | None, available: list[str]) -> list[str]:
    """Resolve ``R1-R5``, ``R1,R3``, ``1-3`` or ``all`` against what exists."""
    if not spec or spec.lower() == "all":
        return available
    wanted: list[str] = []
    for piece in re.split(r"[,\s]+", spec.strip()):
        if not piece:
            continue
        match = re.fullmatch(r"[Rr]?(\d+)\s*-\s*[Rr]?(\d+)", piece)
        if match:
            low, high = int(match.group(1)), int(match.group(2))
            if low > high:
                low, high = high, low
            wanted.extend(f"R{index}" for index in range(low, high + 1))
            continue
        match = re.fullmatch(r"[Rr]?(\d+)", piece)
        if match:
            wanted.append(f"R{int(match.group(1))}")
            continue
        raise ValueError(f"could not read replicate selector {piece!r}")
    missing = [name for name in wanted if name not in available]
    if missing:
        raise ValueError("requested replicate(s) not present: " + ", ".join(missing))
    # Preserve the requested order, without duplicates.
    seen: set[str] = set()
    return [name for name in wanted if not (name in seen or seen.add(name))]


def replicate_sort_key(name: str) -> tuple[int, str]:
    match = re.fullmatch(r"[Rr](\d+)", name)
    return (int(match.group(1)), name) if match else (1 << 30, name)


def discover(data_dir: Path, condition_spec: str | None,
             replicate_spec: str | None) -> list[tuple[str, str]]:
    """Return the (condition, replicate) pairs to run, in a stable order."""
    if not data_dir.is_dir():
        raise ValueError(f"--data-dir is not a directory: {data_dir}")

    # A common slip is to point at one condition instead of the folder holding
    # them, which would otherwise read each R* as a condition and report only
    # that nothing matched.
    if any(entry.is_dir() and re.fullmatch(r"[Rr]\d+", entry.name)
           for entry in data_dir.iterdir()):
        raise ValueError(
            f"--data-dir points at a single condition: {data_dir}\n"
            "It expects the folder that CONTAINS the taxa*_gt*_dup*_loss*_pop* "
            "condition directories, so pass the parent:\n"
            f"    --data-dir {data_dir.parent}\n"
            f"and add --conditions {data_dir.name} to run only this one.")

    conditions = sorted(
        entry.name for entry in data_dir.iterdir()
        if entry.is_dir() and not entry.name.startswith(".")
    )
    if condition_spec and condition_spec.lower() != "all":
        requested = [piece for piece in re.split(r"[,\s]+", condition_spec.strip()) if piece]
        unknown = [name for name in requested if name not in conditions]
        if unknown:
            raise ValueError(
                "condition(s) not found under " + str(data_dir) + ": "
                + ", ".join(unknown) + "\navailable: " + (", ".join(conditions) or "(none)"))
        conditions = requested
    if not conditions:
        raise ValueError(f"no condition directories found under {data_dir}")

    jobs: list[tuple[str, str]] = []
    for condition in conditions:
        available = sorted(
            (entry.name for entry in (data_dir / condition).iterdir()
             if entry.is_dir() and re.fullmatch(r"[Rr]\d+", entry.name)),
            key=replicate_sort_key,
        )
        if not available:
            print(f"warning: {condition} has no R* replicate directories; skipping",
                  file=sys.stderr)
            continue
        for replicate in parse_replicates(replicate_spec, available):
            jobs.append((condition, replicate))
    if not jobs:
        raise ValueError("nothing to run after applying the condition/replicate filters")
    return jobs


# ── per-job plumbing ─────────────────────────────────────────────────────────

def output_root(args: argparse.Namespace, condition: str, replicate: str) -> Path:
    """Where this replicate's artifacts are written."""
    if args.out_dir:
        return Path(args.out_dir).expanduser().resolve() / condition / replicate
    return Path(args.data_dir).expanduser().resolve() / condition / replicate / "stelar-pro"


def prepare_input(args: argparse.Namespace, condition: str, replicate: str,
                  destination: Path) -> Path | None:
    """Normalize copy labels to species labels, reusing the file if current.

    Both arms must read byte-identical input, so this runs once per replicate.
    """
    replicate_dir = Path(args.data_dir).expanduser().resolve() / condition / replicate
    raw = replicate_dir / GENETREE_SOURCES[args.genetrees]
    reference = replicate_dir / TRUE_SPECIES_TREE
    prepared = destination / f"input-{args.genetrees}-species-labelled.trees"
    metadata = destination / f"input-{args.genetrees}-species-labelled.json"

    for path, what in ((raw, "gene trees"), (reference, "true species tree")):
        if not path.is_file() or path.stat().st_size == 0:
            print(f"  missing {what}: {path}", file=sys.stderr)
            return None

    if prepared.is_file() and prepared.stat().st_size > 0 \
            and prepared.stat().st_mtime >= raw.stat().st_mtime and not args.fresh:
        print(f"  input: reusing {prepared.name}")
        return prepared

    # No -I here: prepare_gdl_input.py imports its siblings in scripts/, and
    # isolated mode drops the script's own directory from sys.path.
    command = [
        sys.executable, str(SCRIPTS / "prepare_gdl_input.py"),
        "--input", str(raw),
        "--reference", str(reference),
        "--output", str(prepared),
        "--metadata", str(metadata),
    ]
    print(f"  input: {shlex.join(command)}")
    if args.dry_run:
        return prepared
    destination.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(command, cwd=SCRIPTS)
    if result.returncode != 0 or not prepared.is_file() or prepared.stat().st_size == 0:
        print(f"  input preparation failed for {condition}/{replicate}", file=sys.stderr)
        return None
    return prepared


def absorb_measurement(row: dict[str, object], measurement: Path) -> dict | None:
    """Copy time and memory out of a measure.py record into ``row``.

    Returns the parsed record so callers can inspect fields the row does not
    carry, or None if there was nothing readable to absorb.
    """
    if not measurement.is_file():
        return None
    try:
        saved = json.loads(measurement.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

    def mib(key: str) -> float | None:
        value = saved.get(key)
        return None if value is None else value / (1024 * 1024)

    row["elapsed_seconds"] = saved.get("elapsed_seconds")
    # CPU time comes from wait4 on the job itself, so it is present even when
    # sampling is off; it is absent only if the child could not be reaped.
    row["cpu_seconds"] = saved.get("cpu_seconds")
    row["cpu_utilization"] = saved.get("cpu_utilization")
    row["cpu_seconds_per_core"] = saved.get("cpu_seconds_per_core")
    row["cpu_efficiency"] = saved.get("cpu_efficiency")
    row["cpu_count"] = saved.get("cpu_count")
    row["peak_rss_mib"] = mib("peak_rss_bytes")
    row["exit_code"] = saved.get("exit_code")
    # GPU columns stay empty rather than 0 when no GPU was visible: "this run
    # used no device memory" and "nobody looked" must not read the same.
    if saved.get("gpu_visible"):
        row["peak_gpu_mib"] = mib("gpu_peak_bytes")
        row["peak_gpu_device_mib"] = mib("gpu_device_peak_bytes")
    return saved


def absorb_log(row: dict[str, object], condition: str, replicate: str, arm: str,
               stelar_log: Path) -> list[dict[str, object]]:
    """Copy STELAR-Pro's own phase instrumentation into ``row`` and phase rows.

    PhaseLogger samples VmRSS every 50 ms and, in GPU phases, cudaMemGetInfo
    every 50 ms -- twenty times finer than this harness can sample from
    outside. It prints that for a human and keeps nothing, so it is recovered
    from the log rather than re-measured.
    """
    summary = parse_log(stelar_log)
    row["internal_total_ms"] = summary.internal_total_ms
    row["triplet_score"] = summary.triplet_score
    row["log_peak_ram_mib"] = summary.log_peak_ram_mib
    row["log_peak_vram_mib"] = summary.log_peak_vram_mib
    row["phase_count"] = len(summary.phases) or None
    row["unfinished_phase"] = summary.unfinished_phase
    slowest = summary.slowest
    if slowest is not None:
        row["slowest_phase"] = slowest.label
        row["slowest_phase_ms"] = slowest.ms
    return [
        {"condition": condition, "replicate": replicate, "arm": arm,
         "phase_index": phase.index, "phase": phase.label, "device": phase.device,
         "phase_ms": phase.ms,
         "peak_ram_mib": phase.peak_ram_mib,
         "peak_ram_so_far_mib": phase.peak_ram_so_far_mib,
         "peak_vram_mib": phase.peak_vram_mib,
         "peak_vram_so_far_mib": phase.peak_vram_so_far_mib}
        for phase in summary.phases
    ]


ASTRAL_SCORE = re.compile(r"^\s*Score:\s*(?P<value>-?\d+(?:\.\d+)?)\s*$")


def absorb_astral_log(row: dict[str, object], stderr_log: Path) -> None:
    """Recover ASTRAL-Pro3's own objective from its stderr.

    Deliberately a *separate* column from STELAR-Pro's ``triplet_score``: this
    is a quartet score over a different objective, and putting the two in one
    column would invite comparing numbers that are not comparable.
    """
    try:
        text = stderr_log.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return
    # The last Score: line is the final one; earlier ones are search rounds.
    for line in reversed(text.splitlines()):
        found = ASTRAL_SCORE.match(line)
        if found:
            row["quartet_score"] = found.group("value")
            return


def build_inference(args: argparse.Namespace, arm: str, prepared: Path,
                    tree: Path, destination: Path) -> tuple[list[str], Path | None]:
    """The command for one arm, plus its own log file if it writes one."""
    spec = ARMS[arm]
    if spec.kind == "astral":
        astral = Path(args.astral_exe).expanduser()
        if not astral.is_file() or not os.access(astral, os.X_OK):
            raise FileNotFoundError(
                f"ASTRAL-Pro3 not found or not executable: {astral}\n"
                "It is built per machine and not tracked by git, so a fresh\n"
                "checkout has to build it once:\n"
                "    ./ensure_backends.sh --cpu-only\n"
                "or pass --astral-exe PATH to use one you already have.")
        # -t is needed explicitly: ASTRAL-Pro3 defaults to 1 thread, which
        # would quietly make the arm look far more expensive than it is.
        command = [str(astral), "-i", str(prepared), "-o", str(tree),
                   "-t", str(args.threads), *spec.options]
        if args.astral_opts:
            command.extend(shlex.split(args.astral_opts))
        return command, None

    stelar_log = destination / f"{arm}.stelar-pro.log"
    command = [
        "./run.sh", "--input", str(prepared), "--output", str(tree),
        "--threads", str(args.threads), "--xmx", args.xmx,
        *spec.options, f"--{args.backend}",
        "--log-file", str(stelar_log),
    ]
    if args.no_build:
        command.append("--no-build")
    if args.extra_opts:
        command.extend(shlex.split(args.extra_opts))
    return command, stelar_log


def run_arm(args: argparse.Namespace, condition: str, replicate: str, arm: str,
            prepared: Path, destination: Path) -> dict[str, object]:
    """Run one arm under measure.py and score the tree it produced."""
    spec = ARMS[arm]
    tree = destination / f"{arm}.tre"
    measurement = destination / f"{arm}.measure.json"
    stderr_log = destination / f"{arm}.stderr.log"
    row: dict[str, object] = {
        "condition": condition, "replicate": replicate, "arm": arm,
        "genetrees": args.genetrees, "threads": args.threads,
        # Neither heap nor backend applies to ASTRAL-Pro3; recording the
        # harness's values against it would misreport how it was run.
        "xmx": args.xmx if spec.kind == "stelar" else "",
        "backend": args.backend if spec.kind == "stelar" else "cpu",
        "output_tree": str(tree), "measurement_json": str(measurement),
        "finished_at": "", "status": "pending",
    }

    def absorb_own_log() -> None:
        """Pull in whatever instrumentation this kind of arm leaves behind."""
        if spec.kind == "astral":
            absorb_astral_log(row, stderr_log)
        else:
            row["_phases"] = absorb_log(
                row, condition, replicate, arm,
                destination / f"{arm}.stelar-pro.log")

    if tree.is_file() and tree.stat().st_size > 0 and not args.fresh:
        print(f"  {arm}: output exists, scoring only (use --fresh to rerun)")
        row["status"] = "reused"
        absorb_measurement(row, measurement)
        absorb_own_log()
        return score(args, condition, replicate, tree, row)

    try:
        inference, stelar_log = build_inference(
            args, arm, prepared, tree, destination)
    except FileNotFoundError as error:
        row["status"] = "failed"
        print(f"  {arm}: {error}", file=sys.stderr)
        return row

    command = [
        sys.executable, str(MEASURE),
        "--json-out", str(measurement),
        "--label", f"{condition}/{replicate}/{arm}",
        "--stdout", str(destination / f"{arm}.stdout.log"),
        "--stderr", str(destination / f"{arm}.stderr.log"),
        "--cwd", str(ROOT),
        "--sample-interval", str(args.sample_interval),
        "--trace", str(destination / f"{arm}.memory-trace.tsv"),
    ]
    if args.timeout:
        command.extend(["--timeout", str(args.timeout)])
    command.append("--")
    command.extend(inference)

    print(f"  {arm}: {shlex.join(inference)}")
    if args.dry_run:
        row["status"] = "dry-run"
        return row

    destination.mkdir(parents=True, exist_ok=True)
    if stelar_log is not None:
        # run.sh realpaths both --output and --log-file, which fails on a path
        # that does not exist yet, so both are created empty before it starts.
        # ASTRAL-Pro3 does neither, so it needs no pre-touching.
        tree.touch(exist_ok=True)
        stelar_log.touch(exist_ok=True)
    result = subprocess.run(command)
    row["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")

    saved = absorb_measurement(row, measurement)
    absorb_own_log()
    if row.get("unfinished_phase"):
        print(f"  {arm}: STELAR-Pro's log ends inside "
              f"{row['unfinished_phase']!r}", file=sys.stderr)
    if saved is not None:
        if saved.get("timed_out"):
            row["status"] = "timeout"
            print(f"  {arm}: TIMED OUT after {args.timeout} s", file=sys.stderr)
            return row
        # A GPU that was there but held nothing for us is the signature of a
        # silent fallback to a CPU backend, which would otherwise be reported
        # as a clean GPU timing.
        if (spec.kind == "stelar" and args.backend in ("gpu", "auto")
                and saved.get("gpu_visible")
                and not saved.get("gpu_peak_attributed")):
            print(f"  {arm}: warning: a GPU was visible but no device memory "
                  f"was attributable to this run -- check "
                  f"{destination / (arm + '.stelar-pro.log')} for a fallback "
                  f"to a CPU backend", file=sys.stderr)

    produced_tree = tree.is_file() and tree.stat().st_size > 0
    if not produced_tree:
        row["status"] = "failed"
        print(f"  {arm}: FAILED (exit {result.returncode}); see "
              f"{destination / (arm + '.stderr.log')}", file=sys.stderr)
        return row

    if result.returncode != 0:
        # The inference produced a tree, so the wrapper itself failed. Score the
        # tree rather than discarding hours of work over a lost measurement.
        row["status"] = "ok/measurement-failed"
        print(f"  {arm}: tree produced, but the measurement wrapper exited "
              f"{result.returncode}; scoring anyway, time and memory may be "
              f"missing", file=sys.stderr)
    else:
        row["status"] = "ok"
    return score(args, condition, replicate, tree, row)


def score(args: argparse.Namespace, condition: str, replicate: str,
          tree: Path, row: dict[str, object]) -> dict[str, object]:
    """Normalized unrooted RF of ``tree`` against the replicate's true tree."""
    truth = Path(args.data_dir).expanduser().resolve() \
        / condition / replicate / TRUE_SPECIES_TREE
    try:
        result = compare_newicks(read_tree(tree), read_tree(truth))
    except (OSError, ValueError) as error:
        print(f"  RF failed for {row['arm']}: {error}", file=sys.stderr)
        row["status"] = f"{row['status']}/rf-failed"
        return row
    row.update(
        rf_rate=result.rf_rate, raw_rf=result.raw_rf, max_rf=result.max_rf,
        taxa=result.taxa, output_splits=result.output_splits,
        true_splits=result.true_splits, shared_splits=result.shared_splits,
    )
    elapsed = row.get("elapsed_seconds")
    cpu = row.get("cpu_seconds")
    cores = row.get("cpu_utilization")
    per_core = row.get("cpu_seconds_per_core")
    eff = row.get("cpu_efficiency")
    peak = row.get("peak_rss_mib")
    gpu = row.get("peak_gpu_mib")
    print(f"  {row['arm']}: RF {result.rf_rate:.6f} "
          f"({result.raw_rf}/{result.max_rf})"
          + (f", {elapsed:.1f} s wall" if isinstance(elapsed, (int, float)) else "")
          + (f", {cpu:.1f} s CPU" if isinstance(cpu, (int, float)) else "")
          + (f" ({cores:.1f}x)" if isinstance(cores, (int, float)) else "")
          + (f", {per_core:.1f} s/CPU" if isinstance(per_core, (int, float)) else "")
          + (f" = {eff * 100:.0f}% busy" if isinstance(eff, (int, float)) else "")
          + (f", {peak:.0f} MiB host" if isinstance(peak, (int, float)) else "")
          + (f", {gpu:.0f} MiB GPU" if isinstance(gpu, (int, float)) else ""))
    return row


# ── reporting ────────────────────────────────────────────────────────────────

def append_phase_rows(path: Path, rows: list[dict[str, object]]) -> int:
    """Append the per-phase table, writing the header only on first creation."""
    phases = [phase for row in rows for phase in row.get("_phases") or ()]
    if not phases:
        return 0
    fresh = not path.is_file() or path.stat().st_size == 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PHASE_FIELDS,
                                delimiter="\t", extrasaction="ignore")
        if fresh:
            writer.writeheader()
        for phase in phases:
            writer.writerow({key: phase.get(key, "") for key in PHASE_FIELDS})
    return len(phases)


def append_rows(path: Path, rows: list[dict[str, object]]) -> None:
    """Append to results.tsv, writing the header only on first creation."""
    fresh = not path.is_file() or path.stat().st_size == 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS,
                                delimiter="\t", extrasaction="ignore")
        if fresh:
            writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in RESULT_FIELDS})


def short_name(condition: str, replicate: str) -> str:
    """Compact row label: drop the constant pop suffix, keep the replicate."""
    trimmed = re.sub(r"_pop\d+$", "", condition)
    return f"{trimmed}/{replicate}"


def summarise(rows: list[dict[str, object]], arms: list[str]) -> None:
    """Per-arm means, then the paired per-replicate table if both arms ran."""
    scored = [row for row in rows if isinstance(row.get("rf_rate"), (int, float))]
    if not scored:
        print("\nNo scored runs to summarize.")
        return

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    # Show the GPU column only when a run actually measured device memory, so
    # a CPU sweep is not padded with a column of zeros.
    show_gpu = any(isinstance(row.get("peak_gpu_mib"), (int, float)) for row in scored)
    header = (f"{'arm':<18}{'n':>4}{'mean RF':>12}{'wall s':>11}"
              f"{'CPU s':>11}{'s/CPU':>9}{'cores':>8}{'busy%':>8}"
              f"{'host MiB':>11}")
    if show_gpu:
        header += f"{'GPU MiB':>11}"
    width = len(header)
    print(header)
    print("-" * width)
    for arm in arms:
        group = [row for row in scored if row["arm"] == arm]
        if not group:
            continue
        mean = lambda key: (                                      # noqa: E731
            sum(float(row[key]) for row in group if isinstance(row.get(key), (int, float)))
            / max(1, sum(1 for row in group if isinstance(row.get(key), (int, float)))))
        line = (f"{arm:<18}{len(group):>4}{mean('rf_rate'):>12.6f}"
                f"{mean('elapsed_seconds'):>11.1f}{mean('cpu_seconds'):>11.1f}"
                f"{mean('cpu_seconds_per_core'):>9.1f}"
                f"{mean('cpu_utilization'):>8.1f}"
                f"{mean('cpu_efficiency') * 100:>8.0f}"
                f"{mean('peak_rss_mib'):>11.0f}")
        if show_gpu:
            line += f"{mean('peak_gpu_mib'):>11.0f}"
        print(line)
    print("-" * width)

    if len(arms) < 2:
        return
    # Every other arm is compared against the first, which is the baseline by
    # construction: --arms is ordered, so the baseline is whatever is listed
    # first. Previously anything other than exactly two arms printed nothing.
    index = {(row["condition"], row["replicate"], row["arm"]): row for row in scored}
    pairs = sorted({(row["condition"], row["replicate"]) for row in scored})
    for other in arms[1:]:
        paired_table(index, pairs, arms[0], other)


def paired_table(index: dict, pairs: list, left: str, right: str) -> None:
    """One per-replicate table of ``right`` against ``left``."""
    complete = [key for key in pairs
                if (key[0], key[1], left) in index and (key[0], key[1], right) in index]
    if not complete:
        print(f"\nNo replicate has both {left} and {right}; no paired comparison.")
        return

    name_width = max(len("replicate"),
                     *(len(short_name(c, r)) for c, r in complete)) + 2
    arm_width = max(14, len(left) + 2, len(right) + 2)
    rule = "-" * (name_width + 2 * arm_width + 20)
    print(f"\n{'replicate':<{name_width}}{left:>{arm_width}}{right:>{arm_width}}"
          f"{'delta':>12}{'splits':>8}")
    print(rule)
    deltas = []
    for condition, replicate in complete:
        a = index[(condition, replicate, left)]
        b = index[(condition, replicate, right)]
        delta = float(b["rf_rate"]) - float(a["rf_rate"])
        deltas.append(delta)
        # raw_rf counts split differences in both directions on binary trees.
        splits = (int(b["raw_rf"]) - int(a["raw_rf"])) // 2
        print(f"{short_name(condition, replicate):<{name_width}}"
              f"{float(a['rf_rate']):>{arm_width}.6f}"
              f"{float(b['rf_rate']):>{arm_width}.6f}{delta:>+12.6f}{splits:>+8}")
    print(rule)
    better = sum(1 for d in deltas if d < 0)
    worse = sum(1 for d in deltas if d > 0)
    print(f"{'mean':<{name_width}}"
          f"{sum(float(index[(c, r, left)]['rf_rate']) for c, r in complete)/len(complete):>{arm_width}.6f}"
          f"{sum(float(index[(c, r, right)]['rf_rate']) for c, r in complete)/len(complete):>{arm_width}.6f}"
          f"{sum(deltas)/len(deltas):>+12.6f}")
    print(f"\n{right} better on {better}, worse on {worse}, "
          f"tied on {len(deltas) - better - worse} of {len(deltas)} paired replicate(s).")
    print("Negative delta favours " + right + ".")


# ── entry point ──────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", required=True,
                        help="directory holding the taxa*_gt*_dup*_loss*_pop* conditions")
    parser.add_argument("--conditions", default="all",
                        help="comma-separated condition directory names (default: all)")
    parser.add_argument("--replicates", default="all",
                        help="R1-R5, R1,R3 or all (default: all present)")
    parser.add_argument("--arms", default=",".join(DEFAULT_ARMS),
                        help="comma-separated arms from " + ", ".join(ARMS)
                             + f" (default: {','.join(DEFAULT_ARMS)}); the "
                             "first listed is the baseline every other arm is "
                             "compared against in the paired tables")
    parser.add_argument("--genetrees", choices=sorted(GENETREE_SOURCES), default="estimated",
                        help="which gene-tree set to infer from (default: estimated)")
    parser.add_argument("--out-dir",
                        help="mirror <condition>/<replicate>/ under this root instead of "
                             "writing into each replicate's stelar-pro/ directory")
    parser.add_argument("--results",
                        help="results TSV path (default: <out-dir or data-dir>/results.tsv)")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="STELAR-Pro worker threads (default: all cores)")
    parser.add_argument("--xmx", default="16g", help="JVM max heap (default: 16g)")
    parser.add_argument("--backend", choices=["cpu", "gpu", "auto"], default="cpu",
                        help="compute backend for the STELAR-Pro arms; "
                             "ASTRAL-Pro3 has no GPU support and ignores it "
                             "(default: cpu)")
    parser.add_argument("--astral-exe", default=str(ASTRAL_DEFAULT),
                        help=f"ASTRAL-Pro3 binary (default: {ASTRAL_DEFAULT})")
    parser.add_argument("--astral-opts", default="",
                        help="extra options for the astral-pro3 arm, e.g. "
                             "\"-u 0\" for a topology-only output tree")
    parser.add_argument("--timeout", type=float,
                        help="per-run timeout in seconds (default: none)")
    parser.add_argument("--sample-interval", type=float, default=1.0,
                        help="memory sampling period in seconds (default: 1.0)")
    parser.add_argument("--extra-opts", default="",
                        help="extra STELAR-Pro options added to every arm")
    parser.add_argument("--no-build", action="store_true",
                        help="pass --no-build to run.sh; only safe when build/ is current")
    parser.add_argument("--fresh", action="store_true",
                        help="rerun arms whose output tree already exists")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the plan and every command, writing nothing")
    args = parser.parse_args(argv)

    arms = [piece for piece in re.split(r"[,\s]+", args.arms.strip()) if piece]
    unknown = [arm for arm in arms if arm not in ARMS]
    if unknown:
        parser.error("unknown arm(s): " + ", ".join(unknown)
                     + "; choose from " + ", ".join(ARMS))
    if not arms:
        parser.error("--arms selected nothing")

    try:
        jobs = discover(Path(args.data_dir).expanduser().resolve(),
                        args.conditions, args.replicates)
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    results_path = Path(args.results).expanduser().resolve() if args.results else (
        (Path(args.out_dir).expanduser().resolve() if args.out_dir
         else Path(args.data_dir).expanduser().resolve()) / "results.tsv")

    print("=" * 78)
    print("STELAR-Pro experiment: " + " vs ".join(arms))
    print("=" * 78)
    print(f"  repo              {ROOT}")
    print(f"  data              {Path(args.data_dir).expanduser().resolve()}")
    print(f"  gene trees        {args.genetrees} ({GENETREE_SOURCES[args.genetrees]})")
    print(f"  outputs           {args.out_dir or 'in each replicate: stelar-pro/'}")
    print(f"  results           {results_path}")
    print(f"  phase table       {results_path.stem}-phases{results_path.suffix}")
    print(f"  backend/threads   {args.backend}, {args.threads} thread(s), heap {args.xmx}")
    print(f"  timeout           {args.timeout or 'none'}")
    print(f"  jobs              {len(jobs)} replicate(s) x {len(arms)} arm(s) "
          f"= {len(jobs) * len(arms)} run(s)")
    if args.no_build:
        print("  build             --no-build: build/ must already match this branch")
    if args.dry_run:
        print("  MODE              dry run; nothing will be written")
    print()

    rows: list[dict[str, object]] = []
    for number, (condition, replicate) in enumerate(jobs, start=1):
        print(f"[{number}/{len(jobs)}] {condition}/{replicate}")
        destination = output_root(args, condition, replicate)
        prepared = prepare_input(args, condition, replicate, destination)
        if prepared is None:
            for arm in arms:
                rows.append({"condition": condition, "replicate": replicate, "arm": arm,
                             "status": "skipped/input-missing", "genetrees": args.genetrees})
            print()
            continue
        for arm in arms:
            rows.append(run_arm(args, condition, replicate, arm, prepared, destination))
        print()

    if not args.dry_run:
        append_rows(results_path, rows)
        print(f"Appended {len(rows)} row(s) to {results_path}")
        phases_path = results_path.with_name(
            results_path.stem + "-phases" + results_path.suffix)
        written = append_phase_rows(phases_path, rows)
        if written:
            print(f"Appended {written} phase row(s) to {phases_path}")
        summarise(rows, arms)
    else:
        print("Dry run complete; no files were written.")

    failed = [row for row in rows
              if str(row.get("status", "")).startswith(("failed", "timeout", "skipped"))]
    if failed:
        print(f"\n{len(failed)} run(s) did not produce a scored tree:", file=sys.stderr)
        for row in failed:
            print(f"  {row['condition']}/{row['replicate']}/{row['arm']}: "
                  f"{row['status']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
