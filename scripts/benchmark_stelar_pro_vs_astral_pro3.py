#!/usr/bin/env python3
"""Run and summarize the 30-replicate STELAR-Pro/ASTRAL-Pro3 benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
import time

from calculate_rf_rate import compare_newicks, read_tree


DATASETS = (
    ("0.25", "0.8e-10", "simphy-200taxa-1000gt-10rep-dupmean-0.25-lb-0.8e-10-noloss"),
    ("1", "2e-10", "simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss"),
    ("3", "4.1e-10", "simphy-200taxa-1000gt-10rep-dupmean-3-lb-4.1e-10-noloss"),
)
METHODS = ("stelar-pro", "astral-pro3")


def natural_key(path: Path) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", path.name)
    )


def command_for(
    method: str,
    root: Path,
    gene_trees: Path,
    output_tree: Path,
    threads: int,
) -> list[str]:
    if method == "stelar-pro":
        return [
            str(root / "stelar-pro"),
            "--no-build",
            "--gpu-strict",
            "--threads",
            str(threads),
            "--search-space",
            "S1",
            "--search-mode",
            "local",
            "--quiet",
            "--input",
            str(gene_trees),
            "--output",
            str(output_tree),
        ]
    return [
        str(root / "ASTER-Linux/bin/astral-pro3"),
        "--thread",
        str(threads),
        "--seed",
        "42",
        "--verbose",
        "1",
        "--output",
        str(output_tree),
        str(gene_trees),
    ]


def run_one(command: list[str], run_dir: Path, method: str) -> float:
    output_tree = run_dir / f"{method}.tre"
    metadata_path = run_dir / f"{method}.json"
    if output_tree.stat().st_size > 0 if output_tree.exists() else False:
        if metadata_path.exists():
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            if metadata.get("exit_code") == 0:
                print(f"SKIP {run_dir.parent.name}/{run_dir.name} {method}", flush=True)
                return float(metadata["wall_seconds"])

    stdout_path = run_dir / f"{method}.stdout.log"
    stderr_path = run_dir / f"{method}.stderr.log"
    print(f"RUN  {run_dir.parent.name}/{run_dir.name} {method}", flush=True)
    started = time.perf_counter()
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(command, cwd=run_dir, stdout=stdout, stderr=stderr)
    wall_seconds = time.perf_counter() - started
    metadata = {
        "command": command,
        "exit_code": completed.returncode,
        "wall_seconds": wall_seconds,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"{method} failed for {run_dir}: see {stderr_path} and {stdout_path}"
        )
    if not output_tree.is_file() or output_tree.stat().st_size == 0:
        raise RuntimeError(f"{method} produced no tree for {run_dir}")
    print(f"DONE {run_dir.parent.name}/{run_dir.name} {method}: {wall_seconds:.3f}s", flush=True)
    return wall_seconds


def mean_sd(values: list[float]) -> str:
    mean = statistics.fmean(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return f"{mean:.3f} ± {sd:.3f}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=METHODS,
        default=list(METHODS),
        help="methods to run; completed runs are resumed automatically",
    )
    args = parser.parse_args(argv)

    root = args.root.resolve()
    results_dir = (args.results_dir or root / "pro_data/benchmark-stelar-pro-vs-astral-pro3-gpu").resolve()
    summary_path = (args.summary or root / "pro_data/STELAR_PRO_VS_ASTRAL_PRO3_GPU_SUMMARY.md").resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for duplicate_mean, simphy_rate, dirname in DATASETS:
        dataset = root / "pro_data" / dirname
        replicates = sorted((p for p in dataset.glob("R*") if p.is_dir()), key=natural_key)
        if len(replicates) != 10:
            raise RuntimeError(f"expected 10 replicates in {dataset}, found {len(replicates)}")
        for replicate in replicates:
            gene_trees = replicate / "all_gt.tre"
            true_tree = replicate / "s_tree.trees"
            run_dir = results_dir / dirname / replicate.name
            run_dir.mkdir(parents=True, exist_ok=True)
            for method in args.methods:
                output_tree = run_dir / f"{method}.tre"
                command = command_for(method, root, gene_trees, output_tree, args.threads)
                wall_seconds = run_one(command, run_dir, method)
                rf = compare_newicks(read_tree(output_tree), read_tree(true_tree))
                rows.append(
                    {
                        "duplicate_mean": duplicate_mean,
                        "simphy_duplication_rate": simphy_rate,
                        "dataset": dirname,
                        "replicate": replicate.name,
                        "method": method,
                        "threads": args.threads,
                        "wall_seconds": wall_seconds,
                        "raw_rf": rf.raw_rf,
                        "max_rf": rf.max_rf,
                        "rf_rate": rf.rf_rate,
                    }
                )
                print(
                    f"RF   {dirname}/{replicate.name} {method}: "
                    f"{rf.rf_rate:.6f} ({rf.raw_rf}/{rf.max_rf})",
                    flush=True,
                )

    csv_path = results_dir / "replicate-results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary_lines = [
        "# STELAR-Pro vs ASTRAL-Pro3 benchmark",
        "",
        (
            f"Each cell is the mean ± sample SD over 10 replicates (200 taxa and "
            f"1,000 gene trees per replicate). Wall time used {args.threads} CPU threads. "
            "STELAR-Pro used strict CUDA execution with S1/I1 defaults; RF rate is "
            "standard unrooted RF / 2(n−3), where lower is better."
        ),
        "",
        "| Mean duplicates/species | SimPhy duplication rate | STELAR-Pro time (s) | STELAR-Pro RF rate | ASTRAL-Pro3 time (s) | ASTRAL-Pro3 RF rate |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for duplicate_mean, simphy_rate, _ in DATASETS:
        selected = [row for row in rows if row["duplicate_mean"] == duplicate_mean]
        by_method = {
            method: [row for row in selected if row["method"] == method]
            for method in METHODS
        }
        if any(len(method_rows) != 10 for method_rows in by_method.values()):
            raise RuntimeError("cannot create final summary until both methods have 10 replicates")
        stelar_times = [float(row["wall_seconds"]) for row in by_method["stelar-pro"]]
        stelar_rf = [float(row["rf_rate"]) for row in by_method["stelar-pro"]]
        astral_times = [float(row["wall_seconds"]) for row in by_method["astral-pro3"]]
        astral_rf = [float(row["rf_rate"]) for row in by_method["astral-pro3"]]
        summary_lines.append(
            f"| {duplicate_mean} | {simphy_rate} | {mean_sd(stelar_times)} | "
            f"{mean_sd(stelar_rf)} | {mean_sd(astral_times)} | {mean_sd(astral_rf)} |"
        )
    summary_lines.extend(
        [
            "",
            f"Detailed per-replicate measurements: `{csv_path.relative_to(root)}`.",
        ]
    )
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(f"SUMMARY {summary_path}")
    print(f"DETAILS {csv_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(2)
