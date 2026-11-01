#!/usr/bin/env python3
"""Benchmark the default one-pass, local S5/S6 paths on sim200_dup1_loss1.

Uses the exact species-labelled inputs from the earlier ten-replicate comparison.
Results are written to a new directory; source data and earlier runs are untouched.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess

from calculate_rf_rate import compare_newicks, read_tree


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT.parent / "data/sim200_500/sim200_dup1_loss1"
TIME_FORMAT = "wall_seconds=%e\npeak_rss_kib=%M\nexit_code=%x"
FIELDS = ["replicate", "method", "raw_rf", "max_rf", "rf_rate",
          "wall_seconds", "peak_rss_mib", "peak_gpu_vram_mib", "scored_splits"]


def metrics(path: Path) -> dict[str, str]:
    result = {}
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            key, sep, value = line.partition("=")
            if sep:
                result[key] = value
    return result


def run_one(args: argparse.Namespace, rep: int, method: str) -> dict:
    folder = args.output_dir / f"{rep:02d}" / method
    folder.mkdir(parents=True, exist_ok=True)
    tree = folder / "species_tree.tre"
    log = folder / "terminal.log"
    time_file = folder / "resource.txt"
    input_tree = args.data_dir / "stelar-astral-first10-20261020" / f"{rep:02d}" / "input_species_labels.tre"
    truth = args.data_dir / f"{rep:02d}" / "s_trees.tree"
    if not input_tree.is_file() or not truth.is_file():
        raise FileNotFoundError(f"replicate {rep:02d} input or truth is missing")
    command = ["java", "-Xms256m", f"-Xmx{args.xmx}",
               f"-Dstelarpro.home={ROOT}", f"-Djava.library.path={args.native_path}",
               "-cp", str(args.classes), "stelarx.Main",
               "--input", str(input_tree), "--output", str(tree),
               "--gpu-strict", "--threads", str(args.threads),
               "--search-space", method]
    (folder / "command.json").write_text(json.dumps(command, indent=2) + "\n",
                                         encoding="utf-8")
    timed = ["/usr/bin/time", "-f", TIME_FORMAT, "-o", str(time_file), *command]
    with log.open("w", encoding="utf-8") as stream:
        process = subprocess.Popen(timed, cwd=ROOT, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            code = 124
    if code != 0 or not tree.is_file() or not tree.stat().st_size:
        raise RuntimeError(f"{method} replicate {rep:02d} failed (exit {code}); see {log}")
    rf = compare_newicks(read_tree(tree), read_tree(truth))
    timing = metrics(time_file)
    if timing.get("exit_code") != "0":
        raise RuntimeError(f"GNU time did not record a successful run: {time_file}")
    text = log.read_text(encoding="utf-8", errors="replace")
    if ("Phase 5b Cross-tree transitions" in text or "DP cross-tree allocations" in text
            or "Phase 7b" in text):
        raise RuntimeError(f"{method} did not use the requested one-pass local defaults: {log}")
    vram = re.search(r"Max GPU VRAM\s+(\d+) MiB", text)
    scored = re.findall(r"Weight table: (\d+) splits scored", text)
    row = {
        "replicate": f"{rep:02d}", "method": method,
        "raw_rf": rf.raw_rf, "max_rf": rf.max_rf, "rf_rate": rf.rf_rate,
        "wall_seconds": float(timing["wall_seconds"]),
        "peak_rss_mib": int(timing["peak_rss_kib"]) / 1024,
        "peak_gpu_vram_mib": int(vram[1]) if vram else "",
        "scored_splits": int(scored[-1]) if scored else "",
    }
    (folder / "row.json").write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    return row


def write_summary(args: argparse.Namespace) -> None:
    rows = []
    for rep in args.replicates:
        for method in args.methods:
            path = args.output_dir / f"{rep:02d}" / method / "row.json"
            if path.is_file():
                rows.append(json.loads(path.read_text(encoding="utf-8")))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "results.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    report = ["# sim200_dup1_loss1: default S5/S6 (local, one pass)", "",
              "Ten replicate inputs from `stelar-astral-first10-20261020`; "
              f"strict CUDA, {args.threads} CPU threads, {args.xmx} Java heap. "
              "Cross-tree splits are disabled and max search passes is one. "
              "Times include preprocessing and inference. RSS is GNU time peak; "
              "VRAM is STELAR's reported peak allocation. "
              "RF is standard unrooted RF / 396 (201 taxa).",
              *([args.run_note] if args.run_note else []), "",
              "| Method | Runs | Mean RF rate | Mean time (s) | Mean peak RSS (MiB) | Mean peak GPU VRAM (MiB) |",
              "|:--|--:|--:|--:|--:|--:|"]
    for method in args.methods:
        own = [row for row in rows if row["method"] == method]
        if not own:
            continue
        mean = statistics.fmean
        vram = [float(row["peak_gpu_vram_mib"]) for row in own
                if row["peak_gpu_vram_mib"] != ""]
        vram_mean = f"{mean(vram):.1f}" if vram else "—"
        report.append(f"| {method} | {len(own)} | {100 * mean(row['rf_rate'] for row in own):.2f}% "
                      f"| {mean(row['wall_seconds'] for row in own):.2f} "
                      f"| {mean(row['peak_rss_mib'] for row in own):.1f} "
                      f"| {vram_mean} |")
    old_defaults = args.data_dir / "stelar-s5-s6-defaults-20261021"
    old_results = old_defaults / "results.csv"
    if old_results.is_file():
        with old_results.open(newline="", encoding="utf-8") as stream:
            previous = list(csv.DictReader(stream))
        report += ["", "Earlier one-pass local runs (separate compiled build):", "",
                   "| Method | Runs | Mean RF rate | Mean time (s) | Mean peak RSS (MiB) | Mean peak GPU VRAM (MiB) |",
                   "|:--|--:|--:|--:|--:|--:|"]
        for method in args.methods:
            earlier = [row for row in previous if row["method"] == method]
            if earlier:
                report.append(f"| {method} | {len(earlier)} | "
                              f"{100 * statistics.fmean(float(row['rf_rate']) for row in earlier):.2f}% | "
                              f"{statistics.fmean(float(row['wall_seconds']) for row in earlier):.2f} | "
                              f"{statistics.fmean(float(row['peak_rss_mib']) for row in earlier):.1f} | "
                              f"{statistics.fmean(float(row['peak_gpu_vram_mib']) for row in earlier):.1f} |")
        report.append("")
        for method in args.methods:
            complete = [row for row in rows if row["method"] == method]
            if not complete:
                continue
            matches = sum(
                (args.output_dir / row["replicate"] / method / "species_tree.tre").read_bytes()
                == (old_defaults / row["replicate"] / method / "species_tree.tre").read_bytes()
                for row in complete
                if (old_defaults / row["replicate"] / method / "species_tree.tre").is_file()
            )
            report.append(f"{method} tree outputs identical to earlier defaults: "
                          f"{matches}/{len(complete)}.")
        report += ["", "The S5 commands used the same inputs and options, but different "
                   "compiled class directories. Earlier and current S5 topologies "
                   "match exactly; their wall times differ mainly during preprocessing "
                   "and greedy proposal construction. The archived logs do not establish "
                   "the cause of that timing difference."]
    prior = args.data_dir / "stelar-astral-first10-20261020" / "results.csv"
    if prior.is_file():
        with prior.open(newline="", encoding="utf-8") as stream:
            prior_rows = list(csv.DictReader(stream))
        report += ["", "Different configuration: earlier benchmark with full cross-tree "
                   "search and adaptive passes (ASTRAL-Pro3 used CPU). Its S5 RF rate "
                   "must not be compared as a one-pass local result:", "",
                   "| Method | Runs | Mean RF rate | Mean time (s) | Mean peak RSS (MiB) |",
                   "|:--|--:|--:|--:|--:|"]
        for method in ["S5", "S6", "ASTRAL-Pro3"]:
            own = [row for row in prior_rows if row["method"] == method]
            if own:
                report.append(f"| {method} | {len(own)} | "
                              f"{100 * statistics.fmean(float(row['rf_rate']) for row in own):.2f}% | "
                              f"{statistics.fmean(float(row['wall_seconds']) for row in own):.2f} | "
                              f"{statistics.fmean(int(row['peak_host_rss_kib']) / 1024 for row in own):.1f} |")
    report += ["", "| Replicate | Method | RF (raw/396) | Time (s) | Peak RSS (MiB) | GPU VRAM (MiB) | Scored splits |",
               "|:--:|:--|--:|--:|--:|--:|--:|"]
    for row in rows:
        report.append(f"| {row['replicate']} | {row['method']} | "
                      f"{100 * row['rf_rate']:.2f}% ({row['raw_rf']}/{row['max_rf']}) | "
                      f"{row['wall_seconds']:.2f} | {row['peak_rss_mib']:.1f} | "
                      f"{row['peak_gpu_vram_mib'] or '—'} | {row['scored_splits'] or '—'} |")
    (args.output_dir / "SUMMARY.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-dir", type=Path,
                        default=DEFAULT_DATA / "stelar-s5-s6-defaults-20261021")
    parser.add_argument("--classes", type=Path, required=True)
    parser.add_argument("--native-path", default=str(ROOT / "native"),
                        help="Java native library path; directories are colon-separated")
    parser.add_argument("--run-note", default="", help="Additional summary provenance")
    parser.add_argument("--replicates", type=int, nargs="+", default=list(range(1, 11)))
    parser.add_argument("--methods", nargs="+", choices=["S5", "S6"], default=["S5", "S6"])
    parser.add_argument("--threads", type=int, default=28)
    parser.add_argument("--xmx", default="6g")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.classes = args.classes.resolve()
    if not args.classes.is_dir() or any(not 1 <= rep <= 10 for rep in args.replicates):
        parser.error("classes directory must exist; replicates must be 1–10")
    if args.output_dir == args.data_dir or args.data_dir not in args.output_dir.parents:
        parser.error("output directory must be a new child of the dataset directory")
    for rep in args.replicates:
        for method in args.methods:
            path = args.output_dir / f"{rep:02d}" / method / "row.json"
            if path.is_file() and (path.parent / "species_tree.tre").is_file():
                saved = json.loads(path.read_text(encoding="utf-8"))
                print(f"{rep:02d} {method}: already complete, {saved['wall_seconds']:.2f}s", flush=True)
                continue
            row = run_one(args, rep, method)
            write_summary(args)
            print(f"{rep:02d} {method}: RF={row['raw_rf']}/{row['max_rf']} "
                  f"time={row['wall_seconds']:.2f}s RSS={row['peak_rss_mib']:.1f}MiB "
                  f"VRAM={row['peak_gpu_vram_mib']}MiB", flush=True)
    write_summary(args)


if __name__ == "__main__":
    main()
