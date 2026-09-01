#!/usr/bin/env python3
"""Calculate inferred-vs-true species-tree RF rates across replicates."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from pathlib import Path
import re
import statistics

from calculate_rf_rate import compare_newicks, read_tree


RESULT_FIELDS = [
    "replicate",
    "rf_rate",
    "raw_rf",
    "max_rf",
    "taxa",
    "output_splits",
    "true_splits",
    "shared_splits",
    "output_only_splits",
    "true_only_splits",
    "output_tree",
    "true_tree",
]

SUMMARY_FIELDS = [
    "replicates",
    "mean_rf_rate",
    "sample_sd_rf_rate",
    "min_rf_rate",
    "max_rf_rate",
    "mean_raw_rf",
]


def natural_key(path: Path) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", path.name)
    )


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare one inferred STELAR-Pro species tree with the true species "
            "tree in each DATASET/R* replicate directory."
        )
    )
    parser.add_argument("dataset", type=Path, help="dataset containing replicate folders")
    parser.add_argument("--replicate-glob", default="R*", help="replicate folder glob")
    parser.add_argument(
        "--output-tree-file",
        default="out-stelar-pro.tre",
        help=(
            "inferred-tree path relative to each replicate "
            "(default: out-stelar-pro.tre)"
        ),
    )
    parser.add_argument(
        "--true-tree-file",
        default="s_tree.trees",
        help="true-tree path relative to each replicate (default: s_tree.trees)",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        help="results directory (default: DATASET/rf-rate-results)",
    )
    parser.add_argument("--force", action="store_true", help="replace existing CSV outputs")
    args = parser.parse_args(argv)

    if not args.dataset.is_dir():
        parser.error(f"dataset directory does not exist: {args.dataset}")
    replicates = sorted(
        (path for path in args.dataset.glob(args.replicate_glob) if path.is_dir()),
        key=natural_key,
    )
    if not replicates:
        parser.error(
            f"no replicate directories matching {args.replicate_glob!r} "
            f"under {args.dataset}"
        )

    results_dir = args.results_dir or args.dataset / "rf-rate-results"
    rates_path = results_dir / "rf-rates.csv"
    summary_path = results_dir / "rf-rates-summary.csv"
    if not args.force:
        existing = [path for path in (rates_path, summary_path) if path.exists()]
        if existing:
            parser.error("output already exists; use --force: " + ", ".join(map(str, existing)))

    rows: list[dict[str, object]] = []
    for replicate in replicates:
        output_path = replicate / args.output_tree_file
        true_path = replicate / args.true_tree_file
        if not output_path.is_file():
            parser.error(f"inferred species tree does not exist: {output_path}")
        if not true_path.is_file():
            parser.error(f"true species tree does not exist: {true_path}")
        try:
            result = compare_newicks(read_tree(output_path), read_tree(true_path))
        except (OSError, ValueError) as error:
            parser.error(f"cannot compare {replicate.name}: {error}")

        row = {"replicate": replicate.name}
        row.update(asdict(result))
        row["output_tree"] = str(output_path)
        row["true_tree"] = str(true_path)
        rows.append(row)
        print(
            f"{replicate.name}: RF rate={result.rf_rate:.6f} "
            f"({result.raw_rf}/{result.max_rf})"
        )

    rates = [float(row["rf_rate"]) for row in rows]
    raw_distances = [int(row["raw_rf"]) for row in rows]
    summary = {
        "replicates": len(rows),
        "mean_rf_rate": statistics.fmean(rates),
        "sample_sd_rf_rate": statistics.stdev(rates) if len(rates) > 1 else 0.0,
        "min_rf_rate": min(rates),
        "max_rf_rate": max(rates),
        "mean_raw_rf": statistics.fmean(raw_distances),
    }

    results_dir.mkdir(parents=True, exist_ok=True)
    write_csv(rates_path, RESULT_FIELDS, rows)
    write_csv(summary_path, SUMMARY_FIELDS, [summary])
    print(
        f"ALL: replicates={len(rows)}, mean RF rate={summary['mean_rf_rate']:.6f}, "
        f"sample SD={summary['sample_sd_rf_rate']:.6f}"
    )
    print(f"Replicate CSV: {rates_path}")
    print(f"Summary CSV:   {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
