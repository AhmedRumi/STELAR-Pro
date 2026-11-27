#!/usr/bin/env python3
"""Resumable, sequential Java S4/S6 versus C++ s2/s3 benchmark.

Preparation, compilation, RF and parity checks are outside each GNU-time
interval. No candidate dumps are requested during timed inference. GPU memory
is process-attributed nvidia-smi sampling, hence an observed lower bound.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from calculate_rf_rate import compare_newicks, read_tree  # noqa: E402
from parity import topology, score  # noqa: E402

METHODS = (("cpp-s2", "cpp", "s2"), ("java-S4", "java", "S4"),
           ("cpp-s3", "cpp", "s3"), ("java-S6", "java", "S6"))
TIME_FORMAT = "wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\npeak_rss_kib=%M\nexit_code=%x"


def save(path: Path, value) -> None:
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n")
    temp.replace(path)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_digest() -> str:
    h = hashlib.sha256()
    for path in sorted((ROOT / "src").rglob("*.java")):
        h.update(str(path.relative_to(ROOT)).encode())
        h.update(path.read_bytes())
    return h.hexdigest()


def descendants(pid: int) -> set[int]:
    result, pending = {pid}, [pid]
    while pending:
        current = pending.pop()
        try:
            children = Path(f"/proc/{current}/task/{current}/children").read_text().split()
        except (OSError, ProcessLookupError):
            continue
        for child in map(int, children):
            if child not in result:
                result.add(child)
                pending.append(child)
    return result


def gpu_sample(pids: set[int]) -> int | None:
    try:
        p = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,used_gpu_memory",
                            "--format=csv,noheader,nounits"], capture_output=True,
                           text=True, timeout=3)
        if p.returncode:
            return None
        total = 0
        for line in p.stdout.splitlines():
            words = [word.strip() for word in line.split(",")]
            if len(words) == 2 and all(word.isdigit() for word in words):
                if int(words[0]) in pids:
                    total += int(words[1])
        return total
    except (OSError, subprocess.TimeoutExpired):
        return None


def foreign_inference() -> list[dict]:
    found = []
    for folder in Path("/proc").iterdir():
        if not folder.name.isdigit():
            continue
        try:
            name = (folder / "comm").read_text().strip()
            if name in ("java", "astral-pro3", "stelar-pro-cpp", "stelar-pro-cpp-"):
                stat = (folder / "stat").read_text().split(") ", 1)[1]
                if not stat.startswith("Z"):
                    found.append({"pid": int(folder.name), "comm": name})
        except OSError:
            pass
    return found


def datasets(tag: str) -> list[dict]:
    small = ROOT.parent / "data/sim200_500/sim200_dup1_loss1"
    large = ROOT.parent / "data/taxa2000_gt1000_dup1_loss0.1_pop45800000/R1"
    rows = []
    for i in range(1, 11):
        rep = f"{i:02d}"
        rows.append(dict(dataset="sim200_dup1_loss1", replicate=rep,
                         input=str(small / "stelar-astral-first10-20261020" / rep / "input_species_labels.tre"),
                         truth=str(small / rep / "s_trees.tree"),
                         output=str(small / tag / rep)))
    rows.append(dict(dataset="taxa2000_R1", replicate="R1",
                     input=str(large / "stelar-pro/all_gt_stripped.trees"),
                     truth=str(large / "species-tree/s_tree.trees"),
                     output=str(large / "stelar-pro" / tag)))
    for row in rows:
        row["input_sha256"] = digest(Path(row["input"]))
        row["truth_sha256"] = digest(Path(row["truth"]))
    return rows


def parse_metrics(path: Path) -> dict:
    values = {}
    if path.exists():
        for line in path.read_text().splitlines():
            key, sep, value = line.partition("=")
            if sep:
                try:
                    values[key] = float(value) if key.endswith("seconds") else int(value)
                except ValueError:
                    pass
    return values


def run_one(case: dict, method: tuple, args, manifest: dict) -> dict:
    name, engine, mode = method
    folder = Path(case["output"]) / name
    folder.mkdir(parents=True, exist_ok=True)
    tree = folder / "species.tree"
    command = ([str(args.native)] if engine == "cpp" else
               ["java", f"-Xmx{args.xmx}", f"-Dstelarpro.home={ROOT}",
                f"-Djava.library.path={ROOT / 'native'}", "-cp", str(args.classes), "stelarx.Main"])
    command += ["--input", case["input"], "--output", str(tree), "--search-space", mode,
                "--gpu-strict", "--threads", str(args.threads)]
    signature = dict(command=command, input_sha256=case["input_sha256"], truth_sha256=case["truth_sha256"],
                     engine_sha256=manifest["cpp_binary_sha256"] if engine == "cpp" else manifest["java_source_sha256"],
                     native_backend_sha256=manifest["native_backend_sha256"])
    meta = folder / "row.json"
    if meta.exists() and tree.exists():
        saved = json.loads(meta.read_text())
        if saved.get("signature") == signature and saved.get("exit_code") == 0:
            print("REUSE", case["dataset"], case["replicate"], name, flush=True)
            return saved
    others = foreign_inference()
    if others:
        raise RuntimeError(f"Other inference processes are active; refusing overlap: {others}")
    timefile, logfile = folder / "resource.txt", folder / "terminal.log"
    timed = ["/usr/bin/time", "-f", TIME_FORMAT, "-o", str(timefile), *command]
    save(folder / "command.json", timed)
    save(args.work / "current.json", dict(dataset=case["dataset"], replicate=case["replicate"],
                                         method=name, log=str(logfile), started=datetime.now().isoformat(), status="starting"))
    print("START", case["dataset"], case["replicate"], name, str(logfile), flush=True)
    started = time.monotonic()
    peak_gpu, samples = None, 0
    timeout = False
    with logfile.open("w") as stream, (folder / "gpu_samples.csv").open("w") as gpu_file:
        gpu_writer = csv.writer(gpu_file)
        gpu_writer.writerow(["elapsed_seconds", "process_tree_gpu_mib"])
        process = subprocess.Popen(timed, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        save(args.work / "current.json", dict(dataset=case["dataset"], replicate=case["replicate"],
                                             method=name, log=str(logfile), wrapper_pid=process.pid,
                                             started=datetime.now().isoformat(), status="running"))
        try:
            while process.poll() is None:
                elapsed = time.monotonic() - started
                if elapsed > args.timeout:
                    timeout = True
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    break
                sample = gpu_sample(descendants(process.pid))
                if sample is not None:
                    samples += 1
                    peak_gpu = max(peak_gpu or 0, sample)
                    gpu_writer.writerow([f"{elapsed:.4f}", sample])
                    gpu_file.flush()
                time.sleep(args.sample_interval)
        except BaseException:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait()
            raise
        code = 124 if timeout else process.wait()
    result = dict(dataset=case["dataset"], replicate=case["replicate"], method=name, engine=engine,
                  mode=mode, signature=signature, gpu_peak_sampled_mib=peak_gpu,
                  gpu_samples=samples, gpu_sample_interval_seconds=args.sample_interval,
                  tree=str(tree), log=str(logfile), **parse_metrics(timefile))
    # A GNU-time exit code must not hide a timeout handled by this wrapper.
    result["exit_code"] = code
    if code == 0:
        log = logfile.read_text()
        result["rf"] = asdict(compare_newicks(read_tree(tree), read_tree(Path(case["truth"]))))
        result["final_score"] = str(score(log))
        if engine == "cpp":
            stats = json.loads(Path(str(tree) + ".stats.json").read_text())
            result["candidate_clusters"] = stats["clusters"]
            result["weight_calls"] = list(map(int, re.findall(r"Weight calculation: (\d+) scored splits", log)))
            result["weight_seconds"] = list(map(float, re.findall(r"Weight calculation: \d+ scored splits; \w+ ([0-9.]+) s", log)))
            result["reported_seconds"] = stats["seconds"]
        else:
            clusters = re.findall(r"Cluster extraction: \d+ candidates -> (\d+) unique clusters", log)
            streamed = re.findall(r"(\d+) total in X", log)
            result["candidate_clusters"] = int(streamed[-1] if streamed else clusters[-1])
            result["weight_calls"] = list(map(int, re.findall(r"Weight table: (\d+) splits scored", log)))
            result["weight_seconds"] = [int(value) / 1000 for value in re.findall(
                r"Weight table: \d+ splits scored[^\n]* in (\d+) ms", log)]
            seconds = re.findall(r"Total time: (\d+) ms", log)
            result["reported_seconds"] = int(seconds[-1]) / 1000 if seconds else None
        result["final_scored_splits"] = result["weight_calls"][-1]
    save(meta, result)
    print("DONE", case["dataset"], case["replicate"], name,
          json.dumps({k: result.get(k) for k in ("exit_code", "wall_seconds", "peak_rss_kib", "gpu_peak_sampled_mib", "final_score")}), flush=True)
    return result


def parity(case: dict, rows: list[dict]) -> list[dict]:
    pairs = []
    for cpp, java in (("cpp-s2", "java-S4"), ("cpp-s3", "java-S6")):
        c = next((r for r in rows if r["method"] == cpp and r["replicate"] == case["replicate"] and r["dataset"] == case["dataset"]), None)
        j = next((r for r in rows if r["method"] == java and r["replicate"] == case["replicate"] and r["dataset"] == case["dataset"]), None)
        if not c or not j or c["exit_code"] or j["exit_code"]:
            continue
        pairs.append(dict(dataset=case["dataset"], replicate=case["replicate"], cpp=cpp, java=java,
                          same_rooted_topology=topology(read_tree(Path(c["tree"]))) == topology(read_tree(Path(j["tree"]))),
                          same_final_score=c["final_score"] == j["final_score"],
                          same_weight_calls=c["weight_calls"] == j["weight_calls"],
                          same_candidate_count=c["candidate_clusters"] == j["candidate_clusters"],
                          same_rf=c["rf"]["raw_rf"] == j["rf"]["raw_rf"]))
    return pairs


def report(cases: list[dict], rows: list[dict], args) -> None:
    checks = [p for case in cases for p in parity(case, rows)]
    save(args.work / "results.json", rows)
    save(args.work / "parity.json", checks)
    fields = ["dataset", "replicate", "method", "rf_rate", "raw_rf", "max_rf", "wall_seconds",
              "user_seconds", "system_seconds", "peak_rss_mib", "gpu_peak_sampled_mib",
              "candidate_clusters", "final_scored_splits", "weight_seconds_total", "final_score", "exit_code"]
    folders = {args.work, *(Path(case["output"]).parent if case["dataset"] == "sim200_dup1_loss1" else Path(case["output"]) for case in cases)}
    for folder in folders:
        selected = rows if folder == args.work else [r for r in rows if str(folder) in r["tree"]]
        folder.mkdir(parents=True, exist_ok=True)
        with (folder / "results.csv").open("w") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for row in selected:
                out = {key: row.get(key) for key in fields}
                out.update({key: row.get("rf", {}).get(key) for key in ("rf_rate", "raw_rf", "max_rf")})
                out["peak_rss_mib"] = row.get("peak_rss_kib", 0) / 1024
                out["weight_seconds_total"] = sum(row.get("weight_seconds", []))
                writer.writerow(out)
        lines = ["# Java S4/S6 versus C++ s2/s3", "", f"Completed runs: {len(selected)}. Strict GPU, {args.threads} CPU workers, default search settings; Java -Xmx{args.xmx}.",
                 "All runs are fresh for this experiment, sequential and on the same local host. Compilation, input preparation, RF/parity checks and candidate dumps are excluded from timing. Time is GNU-time wall seconds; RAM is GNU-time peak RSS. VRAM is sampled process-attributed memory, not a guaranteed exact peak.", "",
                 "| Dataset | Implementation | n | Mean RF | Mean time (s) | Mean peak RSS (MiB) | Mean sampled GPU RAM (MiB) | Mean final scored splits |",
                 "|---|---|---:|---:|---:|---:|---:|---:|"]
        for dataset in dict.fromkeys(row["dataset"] for row in selected):
            for name, _, _ in METHODS:
                group = [r for r in selected if r["dataset"] == dataset and r["method"] == name and r["exit_code"] == 0]
                if not group:
                    continue
                mean = lambda key: statistics.fmean(r[key] for r in group)
                gpu = [r["gpu_peak_sampled_mib"] for r in group if r.get("gpu_peak_sampled_mib") is not None]
                lines.append(f"| {dataset} | {name} | {len(group)} | {statistics.fmean(r['rf']['rf_rate'] for r in group):.4%} | {mean('wall_seconds'):.2f} | {mean('peak_rss_kib') / 1024:.1f} | {statistics.fmean(gpu) if gpu else float('nan'):.1f} | {mean('final_scored_splits'):.1f} |")
        relevant = checks if folder == args.work else [p for p in checks if any(r["dataset"] == p["dataset"] and r["replicate"] == p["replicate"] for r in selected)]
        good = sum(all(p[k] for k in ("same_rooted_topology", "same_final_score", "same_weight_calls", "same_candidate_count", "same_rf")) for p in relevant)
        lines += ["", f"Paired parity checks: {good}/{len(relevant)} match rooted topology, exact score, every backbone/final scored-split count, candidate count and RF.",
                  "RF retains the simulator outgroup: normalized unrooted RF = raw RF / [2(n−3)], denominator 396 for 201 labels and 3,996 for 2,001 labels.",
                  "", "Per-run commands, logs, resource reports, sampled VRAM CSV, output tree and row.json are in each replicate/method folder. Full per-replicate measurements are in results.csv."]
        (folder / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--tag", default="stelar-cpp-java-s2-s3-20261127")
    p.add_argument("--native", type=Path, default=ROOT / "cpp/bin/stelar-pro-cpp")
    p.add_argument("--threads", type=int, default=28)
    p.add_argument("--xmx", default="8g")
    p.add_argument("--timeout", type=float, default=7200)
    p.add_argument("--sample-interval", type=float, default=.25)
    args = p.parse_args()
    args.work = args.work.resolve()
    args.native = args.native.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    args.classes = args.work / "java-classes"
    java_hash = source_digest()
    args.classes.mkdir(exist_ok=True)
    with (args.work / "javac.log").open("w") as f:
        subprocess.run(["javac", "-sourcepath", str(ROOT / "src"), "-d", str(args.classes), str(ROOT / "src/stelarx/Main.java")],
                       stdout=f, stderr=subprocess.STDOUT, check=True)
    cases = datasets(args.tag)
    manifest = dict(created=datetime.now().isoformat(), java_source_sha256=java_hash,
                    cpp_binary_sha256=digest(args.native), threads=args.threads, xmx=args.xmx,
                    native_backend_sha256={str(path.name): digest(path) for path in sorted((ROOT / "native").glob("*.so"))},
                    cases=cases)
    save(args.work / "manifest.json", manifest)
    rows = []
    report(cases, rows, args)
    for case in cases:
        for method in METHODS:
            rows.append(run_one(case, method, args, manifest))
            report(cases, rows, args)
    checks = [p for case in cases for p in parity(case, rows)]
    failed = [p for p in checks if not all(p[k] for k in ("same_rooted_topology", "same_final_score", "same_weight_calls", "same_candidate_count", "same_rf"))]
    save(args.work / "current.json", dict(status="complete", runs=len(rows), failed_runs=sum(r["exit_code"] != 0 for r in rows), parity_mismatches=failed))
    print("COMPLETE", len(rows), "runs; parity mismatches:", len(failed), flush=True)
    return bool(failed or any(r["exit_code"] for r in rows))


if __name__ == "__main__":
    raise SystemExit(main())
