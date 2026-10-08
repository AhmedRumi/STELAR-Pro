#!/usr/bin/env python3
"""Run S6 on all ten 1-lb replicates and compare with saved S5/A-Pro3 runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics
import subprocess
import time

from calculate_rf_rate import compare_newicks, read_tree


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "pro_data/simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss"
OUTPUT = ROOT / "pro_data/benchmark-stelar-pro-s6-gpu-dupmean-1"
S3 = ROOT / "pro_data/benchmark-stelar-pro-s3-vs-s2-vs-astral-pro3-gpu-dupmean-1"
S4_FIRST = ROOT / "pro_data/benchmark-stelar-pro-s4-vs-s3-gpu-dupmean-1"
S4_LAST = ROOT / "pro_data/benchmark-stelar-pro-s4-single-triplet-gpu-dupmean-1"
S5 = ROOT / "pro_data/benchmark-stelar-pro-s5-adaptive-unbiased-gpu-dupmean-1"
ASTRAL = ROOT / "pro_data/benchmark-stelar-pro-vs-astral-pro3-gpu" / DATA.name
SUMMARY = ROOT / "pro_data/STELAR_PRO_S6_VS_S5_ASTRAL_PRO3_DUPMEAN_1_GPU_SUMMARY.md"


def run_replicate(replicate: int, threads: int, timeout: int, rounds: int) -> dict:
    name = f"R{replicate}"
    source = DATA / name
    folder = OUTPUT / name
    folder.mkdir(parents=True, exist_ok=True)
    tree = folder / "S6.tre"
    metadata = folder / "S6.json"
    command = [
        str(ROOT / "run.sh"), "--no-build", "--no-notify", "--gpu-strict",
        "--threads", str(threads), "--xmx", "8g", "--search-space", "S6",
        "--s6-rounds", str(rounds), "--s6-seed", "233", "--quiet",
        "--input", str(source / "all_gt.tre"), "--output", str(tree),
    ]
    if metadata.exists() and tree.exists() and tree.stat().st_size:
        saved = json.loads(metadata.read_text(encoding="utf-8"))
        if saved.get("exit_code") == 0 and saved.get("command") == command:
            return saved
    started = time.perf_counter()
    with (folder / "S6.stdout.log").open("w", encoding="utf-8") as stdout, \
         (folder / "S6.stderr.log").open("w", encoding="utf-8") as stderr:
        try:
            result = subprocess.run(command, cwd=ROOT, stdout=stdout,
                                    stderr=stderr, timeout=timeout, check=False)
            code = result.returncode
        except subprocess.TimeoutExpired:
            code = 124
    saved = {
        "command": command,
        "exit_code": code,
        "wall_seconds": time.perf_counter() - started,
        "timeout_seconds": timeout,
    }
    metadata.write_text(json.dumps(saved, indent=2) + "\n", encoding="utf-8")
    return saved


def baseline(replicate: int, method: str) -> tuple[int, float]:
    name = f"R{replicate}"
    truth = read_tree(DATA / name / "s_tree.trees")
    if method == "S3":
        folder = S3 / name
        seconds = float((folder / "S3.wall-seconds.txt").read_text(
            encoding="utf-8").strip())
        tree = folder / "S3.tre"
    elif method == "S4":
        folder = (S4_FIRST if replicate <= 5 else S4_LAST) / name
        seconds = float((folder / "S4.wall-seconds.txt").read_text(
            encoding="utf-8").strip())
        tree = folder / "S4.tre"
    elif method == "S5":
        folder = S5 / name
        seconds_file = folder / ("S5.retry300.wall-seconds.txt"
                                 if replicate == 3 else "S5.wall-seconds.txt")
        seconds = float(seconds_file.read_text(encoding="utf-8").strip())
        tree = folder / "S5.tre"
    else:
        folder = ASTRAL / name
        seconds = json.loads((folder / "astral-pro3.json").read_text(
            encoding="utf-8"))["wall_seconds"]
        tree = folder / "astral-pro3.tre"
    rf = compare_newicks(read_tree(tree), truth)
    return rf.raw_rf, float(seconds)


def mean_sd(values: list[float]) -> str:
    return f"{statistics.mean(values):.2f} ± {statistics.stdev(values):.2f}"


def write_summary(rows: list[dict], rounds: int, threads: int,
                  timeout: int) -> None:
    complete = [row for row in rows if row["S6 RF"] is not None]
    text = [
        "# S6 versus S3–S5 and ASTRAL-Pro3: 200 taxa, duplication mean 1",
        "",
        "Standard unrooted RF is measured against each true single-copy species "
        "tree (maximum RF 394). Lower is better. Wall times include full preprocessing "
        "and inference. S6 used strict CUDA, " + str(threads) + " CPU threads, "
        + str(rounds) + " greedy trees per DP pass, seed 233, and a "
        + str(timeout) + "-second "
        "per-replicate cap. Comparator times come from separate earlier runs "
        "on this machine, not a simultaneous timing trial.",
        "",
        "| Replicate | S6 RF /394 | S6 time (s) | S5 RF /394 | S5 time (s) "
        "| A-Pro3 RF /394 | A-Pro3 time (s) |",
        "|:--|--:|--:|--:|--:|--:|--:|",
    ]
    for row in rows:
        s6_rf = "—" if row["S6 RF"] is None else str(row["S6 RF"])
        s6_time = "—" if row["S6 time"] is None else f"{row['S6 time']:.2f}"
        text.append(f"| {row['replicate']} | {s6_rf} | {s6_time} "
            f"| {row['S5 RF']} | {row['S5 time']:.2f} "
            f"| {row['A-Pro3 RF']} | {row['A-Pro3 time']:.2f} |")
    if len(complete) == 10:
        text.extend([
            "", "| Method | Mean RF rate ± SD | Mean wall time ± SD |",
            "|:--|--:|--:|",
        ])
        for label, rf_key, time_key in [
            ("S6", "S6 RF", "S6 time"),
            ("S3", "S3 RF", "S3 time"),
            ("S4", "S4 RF", "S4 time"),
            ("S5", "S5 RF", "S5 time"),
            ("ASTRAL-Pro3", "A-Pro3 RF", "A-Pro3 time"),
        ]:
            rf = [100.0 * row[rf_key] / 394 for row in rows]
            seconds = [row[time_key] for row in rows]
            text.append(f"| {label} | {mean_sd(rf)}% | {mean_sd(seconds)} s |")
        s6_s3, s6_s4, s6_s5 = [], [], []
        for replicate in range(1, 11):
            name = f"R{replicate}"
            s6_tree = read_tree(OUTPUT / name / "S6.tre")
            s6_s3.append(compare_newicks(s6_tree,
                read_tree(S3 / name / "S3.tre")).raw_rf)
            s4_folder = (S4_FIRST if replicate <= 5 else S4_LAST) / name
            s6_s4.append(compare_newicks(s6_tree,
                read_tree(s4_folder / "S4.tre")).raw_rf)
            s6_s5.append(compare_newicks(s6_tree,
                read_tree(S5 / name / "S5.tre")).raw_rf)
        differences = ", ".join(f"R{i}: RF {distance}"
            for i, distance in enumerate(s6_s5, 1) if distance)
        text.extend(["", f"Pairwise output-tree check: S6 versus S3 has "
            f"RF 0 in {s6_s3.count(0)}/10 replicates; versus S4, "
            f"{s6_s4.count(0)}/10. S6 differs from S5 in "
            f"{differences or 'none'}." ])
    else:
        text.extend(["", f"S6 completed {len(complete)}/10 replicates; missing "
                     "values are not included in an aggregate."])
    text.extend(["", "The per-replicate table emphasizes S6 against its direct "
                 "S5 predecessor and ASTRAL-Pro3; S3/S4 are included in the "
                 "aggregate. S4 uses the revised single-triplet-matrix run for "
                 "R6–R10 and the earlier S4 run for R1–R5.", "",
                 "S6 proposals use only the rooted-triplet support of the "
                 "original tagged gene trees; S1 candidates are retained, and "
                 "all final DP weights use the same original genes.", ""])
    SUMMARY.write_text("\n".join(text), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=28)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--rounds", type=int, default=4)
    args = parser.parse_args()
    rows = []
    for replicate in range(1, 11):
        name = f"R{replicate}"
        print(f"RUN {name}", flush=True)
        result = run_replicate(replicate, args.threads, args.timeout, args.rounds)
        truth = read_tree(DATA / name / "s_tree.trees")
        tree = OUTPUT / name / "S6.tre"
        s6_rf = None
        if result["exit_code"] == 0 and tree.exists() and tree.stat().st_size:
            s6_rf = compare_newicks(read_tree(tree), truth).raw_rf
        s3_rf, s3_time = baseline(replicate, "S3")
        s4_rf, s4_time = baseline(replicate, "S4")
        s5_rf, s5_time = baseline(replicate, "S5")
        apro_rf, apro_time = baseline(replicate, "ASTRAL")
        rows.append({
            "replicate": name, "S6 RF": s6_rf,
            "S6 time": result["wall_seconds"] if s6_rf is not None else None,
            "S3 RF": s3_rf, "S3 time": s3_time,
            "S4 RF": s4_rf, "S4 time": s4_time,
            "S5 RF": s5_rf, "S5 time": s5_time,
            "A-Pro3 RF": apro_rf, "A-Pro3 time": apro_time,
        })
        write_summary(rows, args.rounds, args.threads, args.timeout)
        print(f"DONE {name}: exit={result['exit_code']} RF={s6_rf} "
              f"time={result['wall_seconds']:.2f}s", flush=True)
    return 0 if all(row["S6 RF"] is not None for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
