#!/usr/bin/env python3
"""Sequential S4 GPU benchmark on the first ten sim200_dup1_loss1 replicates.

Uses the same species-labelled input and true species tree as the archived
S2/S5/S6/ASTRAL-Pro3 comparison. Compilation is outside the timed process.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import statistics
import subprocess

from calculate_rf_rate import compare_newicks, read_tree


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT.parent / "data/sim200_500/sim200_dup1_loss1"
OUTPUT = DATA / "stelar-s4-streamed-first10-20261103"
PRIOR = DATA / "stelar-astral-first10-20261020/results.csv"
DEFAULTS = DATA / "stelar-s5-s6-defaults-20261021/results.csv"
FIELDS = ("replicate", "method", "raw_rf", "max_rf", "rf_rate",
          "wall_seconds", "peak_rss_mib", "peak_gpu_vram_mib",
          "candidate_clusters", "scored_splits")


def run_one(rep: int, timeout: int) -> dict:
    name = f"{rep:02d}"
    folder = OUTPUT / name / "S4"
    folder.mkdir(parents=True, exist_ok=True)
    tree = folder / "species_tree.tre"
    log = folder / "terminal.log"
    resource = folder / "resource.txt"
    input_tree = DATA / "stelar-astral-first10-20261020" / name / "input_species_labels.tre"
    truth = DATA / name / "s_trees.tree"
    if not input_tree.is_file() or not truth.is_file():
        raise FileNotFoundError(f"missing input or true tree for {name}")
    command = ["java", "-Xms256m", "-Xmx6g",
               f"-Dstelarpro.home={ROOT}", f"-Djava.library.path={ROOT / 'native'}",
               "-cp", str(ROOT / "build"), "stelarx.Main",
               "--input", str(input_tree), "--output", str(tree),
               "--gpu-strict", "--threads", "28", "--search-space", "S4"]
    (folder / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    timed = ["/usr/bin/time", "-f",
             "wall_seconds=%e\npeak_rss_kib=%M\nexit_code=%x", "-o", str(resource),
             *command]
    with log.open("w") as stream:
        subprocess.run(timed, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                       timeout=timeout, check=True)
    if not tree.is_file() or not tree.stat().st_size:
        raise RuntimeError(f"S4 {name} produced no species tree; see {log}")
    text = log.read_text(errors="replace")
    if ("S4 streamed candidates:" not in text
            or "Phase 3.5 Greedy consensus" in text
            or "Phase 5b Cross-tree transitions" in text):
        raise RuntimeError(f"S4 {name} did not follow the new local streamed path: {log}")
    timing = dict(line.split("=", 1) for line in resource.read_text().splitlines()
                  if "=" in line)
    if timing.get("exit_code") != "0":
        raise RuntimeError(f"S4 {name} failed; see {log}")
    rf = compare_newicks(read_tree(tree), read_tree(truth))
    vram = re.search(r"Max GPU VRAM\s+(\d+) MiB", text)
    clusters = re.search(r"S4 streamed candidates:.*?; (\d+) total in X", text, re.S)
    splits = re.findall(r"Weight table: (\d+) splits scored", text)
    row = dict(replicate=name, method="S4-local-triplet", raw_rf=rf.raw_rf,
               max_rf=rf.max_rf, rf_rate=rf.rf_rate,
               wall_seconds=float(timing["wall_seconds"]),
               peak_rss_mib=int(timing["peak_rss_kib"]) / 1024,
               peak_gpu_vram_mib=int(vram[1]) if vram else "",
               candidate_clusters=int(clusters[1]) if clusters else "",
               scored_splits=int(splits[-1]) if splits else "")
    (folder / "row.json").write_text(json.dumps(row, indent=2) + "\n")
    return row


def report() -> None:
    rows = [json.loads(path.read_text()) for path in
            (OUTPUT / f"{rep:02d}" / "S4" / "row.json" for rep in range(1, 11))
            if path.is_file()]
    with (OUTPUT / "results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    if not rows:
        return
    prior = list(csv.DictReader(PRIOR.open(newline="")))
    defaults = list(csv.DictReader(DEFAULTS.open(newline="")))
    methods = [("S4-local-triplet", rows, "current, streamed"),
               ("S2-no-cross", [r for r in prior if r["method"] == "S2-no-cross"],
                "archived 2026-10-20"),
               ("S5-local", [r for r in defaults if r["method"] == "S5"],
                "archived 2026-10-21"),
               ("S6-local", [r for r in defaults if r["method"] == "S6"],
                "archived 2026-10-21"),
               ("ASTRAL-Pro3", [r for r in prior if r["method"] == "ASTRAL-Pro3"],
                "archived 2026-10-20")]
    lines = ["# S4 streamed triplet guide: sim200_dup1_loss1, estimated gene trees", "",
             f"Completed {len(rows)}/10 replicates. S4 ran alone, strict CUDA, 28 CPU threads, "
             "6 GiB Java heap. Timed Java process includes rooting/tagging and inference, "
             "but excludes compilation and input preparation. Peak RSS includes Java/native "
             "host memory; GPU VRAM is the program-reported allocation peak. RF is "
             "standard unrooted RF / 396 (201 taxa). Archived comparator timings were "
             "measured on earlier builds/dates, so use them as indicative, not a controlled "
             "speed comparison. S4 retains S1 candidates; its extra DISCO/UPGMA candidates "
             "use a triplet pair guide, bounded completion, local transitions, and no "
             "consensus or cross-tree splits. Final weights use only original tagged genes. "
             "Completion still has a four-point insertion decision.", "",
             "| Method | Runs | Mean RF rate | Mean time (s) | Mean peak RSS (MiB) | "
             "Mean peak GPU VRAM (MiB) | Source |",
             "|:--|--:|--:|--:|--:|--:|:--|"]
    for name, group, source in methods:
        if not group:
            continue
        def mean(field: str) -> float:
            return statistics.fmean(float(r[field]) for r in group if r[field] not in ("", None))
        vram = (f"{mean('peak_gpu_vram_mib'):.1f}"
                if any(r.get("peak_gpu_vram_mib") not in ("", None) for r in group) else "—")
        rss_field = "peak_rss_mib" if name in ("S4-local-triplet", "S5-local", "S6-local") \
            else "peak_host_rss_kib"
        rss = mean(rss_field) if rss_field == "peak_rss_mib" else mean(rss_field) / 1024
        lines.append(f"| {name} | {len(group)} | {100 * mean('rf_rate'):.2f}% | "
                     f"{mean('wall_seconds'):.2f} | {rss:.1f} | {vram} | {source} |")
    lines += ["", "| Replicate | S4 RF | S4 time (s) | S4 RSS (MiB) | "
              "S4 VRAM (MiB) | Candidate clusters | Scored splits |",
              "|:--:|--:|--:|--:|--:|--:|--:|"]
    for row in rows:
        lines.append(f"| {row['replicate']} | {100 * row['rf_rate']:.2f}% | "
                     f"{row['wall_seconds']:.2f} | {row['peak_rss_mib']:.1f} | "
                     f"{row['peak_gpu_vram_mib']} | {row['candidate_clusters']} | "
                     f"{row['scored_splits']} |")
    by_method = {name: {r["replicate"]: r for r in group}
                 for name, group, _ in methods}
    lines += ["", "| Replicate | S4 RF | S2-no-cross RF | S5-local RF | "
              "S6-local RF | ASTRAL-Pro3 RF |",
              "|:--:|--:|--:|--:|--:|--:|"]
    for row in rows:
        rep = row["replicate"]
        rates = [f"{100 * float(by_method[name][rep]['rf_rate']):.2f}%"
                 for name in ("S4-local-triplet", "S2-no-cross", "S5-local",
                              "S6-local", "ASTRAL-Pro3")]
        lines.append(f"| {rep} | " + " | ".join(rates) + " |")
    (OUTPUT / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, nargs="+", default=list(range(1, 11)))
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if any(not 1 <= rep <= 10 for rep in args.replicates):
        parser.error("replicates must be 1..10")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for rep in args.replicates:
        saved = OUTPUT / f"{rep:02d}" / "S4" / "row.json"
        if saved.is_file():
            row = json.loads(saved.read_text())
        else:
            row = run_one(rep, args.timeout)
        report()
        print(f"{rep:02d}: RF {row['raw_rf']}/{row['max_rf']}, "
              f"{row['wall_seconds']:.2f}s, RSS {row['peak_rss_mib']:.1f} MiB", flush=True)


if __name__ == "__main__":
    main()
