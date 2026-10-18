#!/usr/bin/env python3
"""Measure peak process RSS (and sampled CUDA VRAM) on the 1-lb dataset.

Results go to a new directory; earlier timing/RF benchmarks are never changed.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import signal
import statistics
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "pro_data/simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss"
RESULTS = ROOT / "pro_data/benchmark-memory-s4-s6-astral-pro3-dupmean-1"
TIME_FORMAT = "wall_seconds=%e\nmax_rss_kb=%M\nexit_status=%x"


def descendant_pids(parent: int) -> set[int]:
    """Find the timed command, including Java after run.sh execs into it."""
    found = {parent}
    pending = [parent]
    while pending:
        pid = pending.pop()
        try:
            children = Path(f"/proc/{pid}/task/{pid}/children").read_text().split()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        for child in children:
            value = int(child)
            if value not in found:
                found.add(value)
                pending.append(value)
    return found


def gpu_memory_mib(pids: set[int]) -> int | None:
    """Sample process-attributed VRAM; polling cannot guarantee the exact peak."""
    try:
        result = subprocess.run([
            "nvidia-smi", "--query-compute-apps=pid,used_gpu_memory",
            "--format=csv,noheader,nounits",
        ], capture_output=True, text=True, timeout=3, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    total = 0
    for row in result.stdout.splitlines():
        parts = [part.strip() for part in row.split(",")]
        if len(parts) != 2 or not parts[0].isdigit() or not parts[1].isdigit():
            continue
        if int(parts[0]) in pids:
            total += int(parts[1])
    return total


def command(method: str, replicate: int, output: Path, backend: str,
            threads: int, xmx: str) -> list[str]:
    source = DATA / f"R{replicate}" / "all_gt.tre"
    if method == "ASTRAL-Pro3":
        return [str(ROOT / "ASTER-Linux/bin/astral-pro3"),
                "--thread", str(threads), "--seed", "42", "--verbose", "1",
                "--output", str(output), str(source)]
    gpu_flag = "--gpu-strict" if backend == "gpu" else "--cpu"
    result = [str(ROOT / "run.sh"), "--no-build", "--no-notify", gpu_flag,
              "--threads", str(threads), "--xmx", xmx,
              "--search-space", method, "--quiet"]
    if method == "S6":
        result += ["--s6-rounds", "4", "--s6-seed", "233"]
    return result + ["--input", str(source), "--output", str(output)]


def parse_time(path: Path) -> dict[str, float | int]:
    if not path.exists():
        return {}
    parsed: dict[str, float | int] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        key, _, value = line.partition("=")
        if key == "wall_seconds":
            parsed[key] = float(value)
        elif key in ("max_rss_kb", "exit_status"):
            parsed[key] = int(value)
    return parsed


def run_one(method: str, replicate: int, args: argparse.Namespace) -> dict:
    folder = RESULTS / args.backend / f"R{replicate}"
    folder.mkdir(parents=True, exist_ok=True)
    stem = method.lower().replace("-", "_")
    tree = folder / f"{stem}.tre"
    metadata = folder / f"{stem}.json"
    time_file = folder / f"{stem}.time"
    cmd = command(method, replicate, tree, args.backend, args.threads, args.xmx)
    if metadata.exists() and tree.exists() and tree.stat().st_size:
        saved = json.loads(metadata.read_text(encoding="utf-8"))
        if (saved.get("command") == cmd and saved.get("exit_code") == 0
                and saved.get("max_rss_kb", 0) > 0):
            return saved

    # GNU time records the child process's high-water RSS. For STELAR-Pro,
    # run.sh execs Java, so the timed child is the actual JVM, not a shell.
    timed = ["/usr/bin/time", "-f", TIME_FORMAT, "-o", str(time_file), *cmd]
    started = time.monotonic()
    peak_gpu = None
    samples = 0
    with (folder / f"{stem}.stdout.log").open("w", encoding="utf-8") as out, \
         (folder / f"{stem}.stderr.log").open("w", encoding="utf-8") as err:
        process = subprocess.Popen(timed, cwd=ROOT, stdout=out, stderr=err,
                                   start_new_session=True)
        timed_out = False
        while process.poll() is None:
            if time.monotonic() - started > args.timeout:
                timed_out = True
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                break
            if method != "ASTRAL-Pro3" and args.backend == "gpu":
                sample = gpu_memory_mib(descendant_pids(process.pid))
                if sample is not None:
                    samples += 1
                    peak_gpu = max(peak_gpu or 0, sample)
            time.sleep(0.5)
        code = 124 if timed_out else process.wait()
    metrics = parse_time(time_file)
    saved = {
        "command": cmd,
        "exit_code": code,
        "wall_seconds": metrics.get("wall_seconds", round(time.monotonic() - started, 3)),
        "max_rss_kb": metrics.get("max_rss_kb"),
        "gpu_peak_sampled_mib": peak_gpu,
        "gpu_samples": samples,
        "gpu_sampling_seconds": 0.5 if method != "ASTRAL-Pro3" and args.backend == "gpu" else None,
        "timeout_seconds": args.timeout,
    }
    metadata.write_text(json.dumps(saved, indent=2) + "\n", encoding="utf-8")
    return saved


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods", nargs="+", choices=["S4", "S6", "ASTRAL-Pro3"],
                        default=["S4", "S6", "ASTRAL-Pro3"])
    parser.add_argument("--replicates", nargs="+", type=int, default=list(range(1, 11)))
    parser.add_argument("--backend", choices=["gpu", "cpu"], default="gpu")
    parser.add_argument("--threads", type=int, default=28)
    parser.add_argument("--xmx", default="8g")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if any(rep < 1 or rep > 10 for rep in args.replicates):
        parser.error("replicates must be between 1 and 10")
    if args.backend == "gpu" and any(m != "ASTRAL-Pro3" for m in args.methods):
        probe = subprocess.run(["nvidia-smi", "-L"], capture_output=True,
                               text=True, check=False)
        if probe.returncode != 0:
            parser.error("GPU unavailable: " + (probe.stderr.strip() or probe.stdout.strip()))
    for rep in args.replicates:
        for method in args.methods:
            saved = run_one(method, rep, args)
            rss = saved.get("max_rss_kb")
            rss_mib = "n/a" if rss is None else f"{rss / 1024:.1f} MiB"
            print(f"R{rep} {method}: exit={saved['exit_code']} "
                  f"RSS={rss_mib} GPU={saved['gpu_peak_sampled_mib']} MiB "
                  f"wall={saved['wall_seconds']:.2f}s", flush=True)
    for method in args.methods:
        values = []
        for rep in args.replicates:
            stem = method.lower().replace("-", "_")
            path = RESULTS / args.backend / f"R{rep}" / f"{stem}.json"
            if not path.exists():
                continue
            data = json.loads(path.read_text(encoding="utf-8"))
            if data.get("exit_code") == 0 and data.get("max_rss_kb"):
                values.append(data["max_rss_kb"] / 1024)
        if values:
            sd = statistics.stdev(values) if len(values) > 1 else 0.0
            print(f"{method}: peak RSS {statistics.mean(values):.1f} ± {sd:.1f} "
                  f"MiB (n={len(values)})", flush=True)


if __name__ == "__main__":
    main()
