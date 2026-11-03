#!/usr/bin/env python3
"""Compare local, one-pass S2/S5/S7 on the same ten sim200 gene-tree inputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import statistics

from benchmark_sim200_s5_s6_defaults import ROOT, DEFAULT_DATA, run_one


METHODS = ("S2", "S5", "S7")
FIELDS = ("replicate", "method", "raw_rf", "max_rf", "rf_rate",
          "wall_seconds", "peak_rss_mib", "peak_gpu_vram_mib",
          "candidate_clusters", "dp_score", "scored_splits")


def candidates(log: Path) -> int:
    text = log.read_text(encoding="utf-8", errors="replace")
    streamed = re.findall(r"streamed candidates:.*?\b(\d+) total in X", text)
    extracted = re.findall(r"Cluster extraction: \d+ candidates -> (\d+) unique clusters", text)
    found = streamed or extracted
    if not found:
        raise RuntimeError(f"candidate count missing from {log}")
    return int(found[-1])


def dp_score(log: Path) -> int:
    text = log.read_text(encoding="utf-8", errors="replace")
    found = re.findall(r"Inference DP: optimization-objective triplet score = (\d+)", text)
    if not found:
        raise RuntimeError(f"optimization score missing from {log}")
    return int(found[-1])


def summarize(args: argparse.Namespace) -> None:
    rows = []
    for rep in args.replicates:
        for method in METHODS:
            path = args.output_dir / f"{rep:02d}" / method / "row.json"
            if path.is_file():
                row = json.loads(path.read_text(encoding="utf-8"))
                if "candidate_clusters" not in row or "dp_score" not in row:
                    log = path.parent / "terminal.log"
                    row["candidate_clusters"] = candidates(log)
                    row["dp_score"] = dp_score(log)
                    path.write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
                rows.append(row)
    with (args.output_dir / "results.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    report = ["# S7 versus S2 and S5: sim200_dup1_loss1", "",
              "Selected from the first ten replicates; identical species-labelled inputs, strict CUDA, "
              f"{args.threads} CPU threads, {args.xmx} Java heap, local DP, one pass, "
              "no cross-tree splits. RF = unrooted RF / 396. Times include "
              "preprocessing; memory is peak host RSS and application-reported GPU VRAM.",
              "", "| Method | Runs | Mean RF | Mean time (s) | Mean RSS (MiB) | "
              "Mean VRAM (MiB) | Mean candidates |",
              "|:--|--:|--:|--:|--:|--:|--:|"]
    for method in METHODS:
        own = [row for row in rows if row["method"] == method]
        if not own:
            continue
        mean = lambda key: statistics.fmean(float(row[key]) for row in own)
        report.append(f"| {method} | {len(own)} | {100 * mean('rf_rate'):.2f}% "
                      f"| {mean('wall_seconds'):.2f} | {mean('peak_rss_mib'):.1f} "
                      f"| {mean('peak_gpu_vram_mib'):.1f} "
                      f"| {mean('candidate_clusters'):.0f} |")
    complete = [
        {row["method"]: row for row in rows if row["replicate"] == f"{rep:02d}"}
        for rep in args.replicates
    ]
    complete = [own for own in complete if len(own) == len(METHODS)]
    if complete:
        extra = statistics.fmean(own["S7"]["candidate_clusters"]
                                 - own["S2"]["candidate_clusters"] for own in complete)
        objective_gains = sum(own["S7"]["dp_score"] > own["S2"]["dp_score"]
                              for own in complete)
        better_rf = sum(own["S7"]["raw_rf"] < own["S2"]["raw_rf"]
                        for own in complete)
        worse_rf = sum(own["S7"]["raw_rf"] > own["S2"]["raw_rf"]
                       for own in complete)
        report += ["", f"S7 adds {extra:.1f} unique clusters over S2 per replicate "
                   f"on average. Its optimized triplet score improves in "
                   f"{objective_gains}/{len(complete)} replicates, while RF improves "
                   f"in {better_rf} and worsens in {worse_rf} relative to S2. "
                   "More candidates therefore did not improve mean RF here."]
    report += ["", "| Rep | S2 RF | S5 RF | S7 RF | S2 s | S5 s | S7 s | "
               "S2 X | S5 X | S7 X |",
               "|:--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
    for rep in args.replicates:
        own = {r["method"]: r for r in rows if r["replicate"] == f"{rep:02d}"}
        if len(own) != len(METHODS):
            continue
        a, b, c = (own[m] for m in METHODS)
        if c["candidate_clusters"] < max(a["candidate_clusters"], b["candidate_clusters"]):
            raise AssertionError(f"S7 lost S2/S5 candidates in replicate {rep:02d}")
        if c["dp_score"] < max(a["dp_score"], b["dp_score"]):
            raise AssertionError(f"S7 optimized score fell below S2/S5 in replicate {rep:02d}")
        report.append(f"| {rep:02d} | {a['raw_rf']}/396 | {b['raw_rf']}/396 "
                      f"| {c['raw_rf']}/396 | {a['wall_seconds']:.2f} "
                      f"| {b['wall_seconds']:.2f} | {c['wall_seconds']:.2f} "
                      f"| {a['candidate_clusters']} | {b['candidate_clusters']} "
                      f"| {c['candidate_clusters']} |")
    report += ["", "S7 combines S5 greedy edge orientations with S2 streamed "
               "DISCO completion and UPGMA candidates. All weights come from "
               "the original rooted/tagged gene trees. Candidate counts are "
               "unique cluster identities in X; DP transition and scored-split "
               "counts are in the per-run logs and CSV. For every completed "
               "replicate, S7's candidate count and optimized triplet score "
               "are at least those of both baselines.", ""]
    (args.output_dir / "SUMMARY.md").write_text("\n".join(report), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-dir", type=Path,
                        default=DEFAULT_DATA / "stelar-s7-vs-s2-s5-first10-20261103")
    parser.add_argument("--classes", type=Path, default=ROOT / "build")
    parser.add_argument("--native-path", default=str(ROOT / "native"))
    parser.add_argument("--threads", type=int, default=28)
    parser.add_argument("--replicates", type=int, nargs="+", default=list(range(1, 11)))
    parser.add_argument("--xmx", default="6g")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.classes = args.classes.resolve()
    if (args.data_dir not in args.output_dir.parents or not args.classes.is_dir()
            or any(not 1 <= rep <= 10 for rep in args.replicates)):
        parser.error("output must be inside the dataset and classes must exist")
    for rep in args.replicates:
        for method in METHODS:
            folder = args.output_dir / f"{rep:02d}" / method
            path = folder / "row.json"
            if path.is_file() and (folder / "species_tree.tre").is_file():
                row = json.loads(path.read_text(encoding="utf-8"))
            else:
                row = run_one(args, rep, method)
                row["candidate_clusters"] = candidates(folder / "terminal.log")
                path.write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
            if "candidate_clusters" not in row or "dp_score" not in row:
                row["candidate_clusters"] = candidates(folder / "terminal.log")
                row["dp_score"] = dp_score(folder / "terminal.log")
                path.write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
            summarize(args)
            print(f"{rep:02d} {method}: RF {row['raw_rf']}/{row['max_rf']}, "
                  f"{row['wall_seconds']:.2f}s, X={row['candidate_clusters']}, "
                  f"VRAM={row['peak_gpu_vram_mib']}MiB", flush=True)


if __name__ == "__main__":
    main()
