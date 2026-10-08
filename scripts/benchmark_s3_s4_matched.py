#!/usr/bin/env python3
"""Matched S3/S4 GPU benchmark on all ten duplication-mean-1 replicates."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics
import subprocess
import time

from calculate_rf_rate import compare_newicks, read_tree


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "pro_data/simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss"
OUTPUT = ROOT / "pro_data/benchmark-stelar-pro-s3-s4-matched-gpu-dupmean-1"
SUMMARY = ROOT / "pro_data/STELAR_PRO_S3_S4_MATCHED_GPU_DUPMEAN_1_SUMMARY.md"
MAX_RF = 2 * (200 - 3)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def phase_seconds(log: Path) -> dict[str, float]:
    current = None
    result = {}
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        start = re.search(r"▶\s+(Phase [\d.]+[^[]+?)\s+\[", line)
        if start:
            current = start.group(1).strip()
        end = re.search(r"✓\s+(\d+)\s+ms", line)
        if end and current:
            result[current] = int(end.group(1)) / 1000.0
            current = None
    return result


def run_one(replicate: int, method: str, threads: int, timeout: int,
            java_hash: str, cuda_hash: str) -> dict:
    name = f"R{replicate}"
    folder = OUTPUT / name
    folder.mkdir(parents=True, exist_ok=True)
    tree = folder / f"{method}.tre"
    metadata = folder / f"{method}.json"
    command = [
        str(ROOT / "run.sh"), "--no-build", "--no-notify", "--gpu-strict",
        "--threads", str(threads), "--xmx", "8g", "--search-space", method,
        "--quiet", "--input", str(DATA / name / "all_gt.tre"),
        "--output", str(tree),
    ]
    if metadata.exists() and tree.exists() and tree.stat().st_size:
        saved = json.loads(metadata.read_text(encoding="utf-8"))
        if (saved.get("exit_code") == 0 and saved.get("command") == command
                and saved.get("java_sha256") == java_hash
                and saved.get("cuda_sha256") == cuda_hash):
            return saved

    started = time.perf_counter()
    with (folder / f"{method}.stdout.log").open("w", encoding="utf-8") as stdout, \
         (folder / f"{method}.stderr.log").open("w", encoding="utf-8") as stderr:
        try:
            completed = subprocess.run(command, cwd=ROOT, stdout=stdout,
                                       stderr=stderr, timeout=timeout, check=False)
            code = completed.returncode
        except subprocess.TimeoutExpired:
            code = 124
    saved = {
        "command": command,
        "exit_code": code,
        "wall_seconds": time.perf_counter() - started,
        "timeout_seconds": timeout,
        "java_sha256": java_hash,
        "cuda_sha256": cuda_hash,
        "phase_seconds": phase_seconds(folder / f"{method}.stderr.log"),
    }
    metadata.write_text(json.dumps(saved, indent=2) + "\n", encoding="utf-8")
    return saved


def fmt(value: float | None) -> str:
    return "—" if value is None else f"{value:.2f}"


def write_summary(rows: list[dict], threads: int, timeout: int,
                  java_hash: str, cuda_hash: str) -> None:
    complete = [r for r in rows if r["S3 RF"] is not None and r["S4 RF"] is not None]
    lines = [
        "# Matched S3 versus S4 GPU benchmark: 200 taxa, duplication mean 1",
        "",
        "Ten replicates, 1,000 gene trees each. Both methods used this build, "
        f"strict CUDA, {threads} CPU threads, 8 GiB Java heap, and a {timeout}-s "
        "per-run cap. Odd replicates ran S3 first; even replicates ran S4 first. "
        "Wall time includes preprocessing and inference. RF is standard unrooted "
        f"RF against the true single-copy species tree (maximum {MAX_RF}; lower is better).",
        "",
        f"Build fingerprints: `Main.class` `{java_hash}`, "
        f"`libstelar_pro_sim.so` `{cuda_hash}`.",
        "",
        "| Replicate | First | S3 time (s) | S4 time (s) | S4 − S3 (s) "
        "| S3 RF | S4 RF | S3↔S4 RF |",
        "|:--|:--|--:|--:|--:|--:|--:|--:|",
    ]
    for row in rows:
        delta = None if row["S3 time"] is None or row["S4 time"] is None \
            else row["S4 time"] - row["S3 time"]
        lines.append(f"| {row['replicate']} | {row['first']} | "
            f"{fmt(row['S3 time'])} | {fmt(row['S4 time'])} | "
            f"{fmt(delta)} | {row['S3 RF'] if row['S3 RF'] is not None else '—'} "
            f"| {row['S4 RF'] if row['S4 RF'] is not None else '—'} "
            f"| {row['pair RF'] if row['pair RF'] is not None else '—'} |")

    if len(complete) == 10:
        s3 = [r["S3 time"] for r in complete]
        s4 = [r["S4 time"] for r in complete]
        delta = [b - a for a, b in zip(s3, s4)]
        guides = {method: [r[f"{method} phases"].get("guide") for r in complete]
                  for method in ("S3", "S4")}
        weights = {method: [r[f"{method} phases"].get("weight") for r in complete]
                   for method in ("S3", "S4")}
        lines.extend([
            "", "| Measure | S3 | S4 | S4 − S3 |",
            "|:--|--:|--:|--:|",
            f"| Mean wall time ± SD | {statistics.mean(s3):.2f} ± {statistics.stdev(s3):.2f} s "
            f"| {statistics.mean(s4):.2f} ± {statistics.stdev(s4):.2f} s "
            f"| {statistics.mean(delta):+.2f} s |",
            f"| Mean RF rate | {statistics.mean(r['S3 RF'] for r in complete) / MAX_RF * 100:.2f}% "
            f"| {statistics.mean(r['S4 RF'] for r in complete) / MAX_RF * 100:.2f}% | — |",
        ])
        for label, values in (("DISCO/completion/matrix/UPGMA phase", guides),
                              ("Triplet weight-scoring phase", weights)):
            a, b = values["S3"], values["S4"]
            if all(v is not None for v in a + b):
                lines.append(f"| Mean {label} | {statistics.mean(a):.2f} s "
                    f"| {statistics.mean(b):.2f} s "
                    f"| {statistics.mean(b) - statistics.mean(a):+.2f} s |")
        odd_delta = statistics.mean(delta[::2])
        even_delta = statistics.mean(delta[1::2])
        lines.extend([
            "", f"Paired S4 − S3 wall-time difference: {statistics.mean(delta):+.2f} "
            f"± {statistics.stdev(delta):.2f} s (mean ± SD; median "
            f"{statistics.median(delta):+.2f} s; odd/S3-first "
            f"{odd_delta:+.2f} s; even/S4-first {even_delta:+.2f} s).",
            f"S3 and S4 output trees have unrooted RF 0 in "
            f"{sum(r['pair RF'] == 0 for r in complete)}/10 replicates.",
        ])
    else:
        lines.extend(["", f"Completed paired runs: {len(complete)}/10. "
                      "Aggregate statistics are withheld until all ten finish."])
    lines.extend(["", "S4 uses one GPU rooted-triplet similarity matrix for "
                  "completion, UPGMA, and consensus; S3 uses its quartet matrix. "
                  "Both use the same original tagged multicopy genes for final "
                  "STELAR-Pro triplet weights.", ""])
    SUMMARY.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=28)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    java_hash = digest(ROOT / "build/stelarx/Main.class")
    cuda_hash = digest(ROOT / "native/libstelar_pro_sim.so")
    rows = []
    for replicate in range(1, 11):
        name = f"R{replicate}"
        first = "S3" if replicate % 2 else "S4"
        order = (first, "S4" if first == "S3" else "S3")
        print(f"RUN {name}: {' then '.join(order)}", flush=True)
        metadata = {}
        for method in order:
            metadata[method] = run_one(replicate, method, args.threads,
                                       args.timeout, java_hash, cuda_hash)
            print(f"  {method}: exit={metadata[method]['exit_code']} "
                  f"wall={metadata[method]['wall_seconds']:.2f}s", flush=True)
        truth = read_tree(DATA / name / "s_tree.trees")
        inferred = {}
        for method in ("S3", "S4"):
            tree = OUTPUT / name / f"{method}.tre"
            inferred[method] = read_tree(tree) if metadata[method]["exit_code"] == 0 \
                and tree.exists() and tree.stat().st_size else None
        row = {"replicate": name, "first": first}
        for method in ("S3", "S4"):
            row[f"{method} RF"] = None if inferred[method] is None \
                else compare_newicks(inferred[method], truth).raw_rf
            row[f"{method} time"] = None if inferred[method] is None \
                else metadata[method]["wall_seconds"]
            phases = metadata[method]["phase_seconds"]
            row[f"{method} phases"] = {
                "guide": next((v for k, v in phases.items() if k.startswith("Phase 1b ")), None),
                "weight": next((v for k, v in phases.items() if k.startswith("Phase 6 ")), None),
            }
        row["pair RF"] = None if None in inferred.values() \
            else compare_newicks(inferred["S3"], inferred["S4"]).raw_rf
        rows.append(row)
        write_summary(rows, args.threads, args.timeout, java_hash, cuda_hash)
    return 0 if all(r["pair RF"] is not None for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
