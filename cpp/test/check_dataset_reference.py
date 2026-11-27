#!/usr/bin/env python3
"""Run a real input sequentially through Java and native retained modes.

Keep commands, GNU time reports, inference logs, trees and exact parity checks.
No inference processes overlap. Compilation must be completed beforehand.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
import time
from parity import topology, score, cluster_set

ROOT = Path(__file__).resolve().parents[2]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--native", type=Path, required=True)
    p.add_argument("--java-classes", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--modes", nargs="+", default=["s1", "s2", "s3"])
    p.add_argument("--gpu", action="store_true")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for mode in a.modes:
        java_mode = {"s1": "S1", "s2": "S4", "s3": "S6"}[mode]
        row = {"mode": mode, "java_mode": java_mode}
        for label in ("java", "native"):
            tree = a.output / f"{mode}-{label}.tree"
            clusters = a.output / f"{mode}-{label}.clusters"
            if label == "java":
                command = ["java", "-Xmx8g", f"-Dstelarpro.home={ROOT}",
                           f"-Djava.library.path={ROOT / 'native'}", "-cp",
                           str(a.java_classes.resolve()), "stelarx.Main"]
            else:
                command = [str(a.native.resolve())]
            command += ["-i", str(a.input.resolve()), "-o", str(tree),
                        "--dump-clusters", str(clusters), "--search-space",
                        java_mode if label == "java" else mode,
                        "--gpu-strict" if a.gpu else "--cpu",
                        "--threads", str(a.threads)]
            timed = ["/usr/bin/time", "-v", "-o", str(a.output / f"{mode}-{label}.time"),
                     *command]
            (a.output / f"{mode}-{label}.command.json").write_text(json.dumps(timed, indent=2) + "\n")
            started = time.monotonic()
            with (a.output / f"{mode}-{label}.log").open("w") as f:
                result = subprocess.run(timed, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT,
                                        timeout=1200)
            if result.returncode:
                raise RuntimeError(f"{mode} {label} failed; inspect retained log")
            log = (a.output / f"{mode}-{label}.log").read_text()
            report = (a.output / f"{mode}-{label}.time").read_text()
            row[label] = {"seconds": time.monotonic() - started,
                          "peak_rss_kib": int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", report)[1]),
                          "score": score(log),
                          "weight_calls": list(map(int, re.findall(
                              r"Weight table: (\d+) splits scored" if label == "java" else
                              r"Weight calculation: (\d+) scored splits", log)))}
            if "no reachable splits to score" in log:
                row[label]["weight_calls"] = [0]
            print(mode, label, row[label], flush=True)
        different = topology((a.output / f"{mode}-java.tree").read_text()) != topology(
            (a.output / f"{mode}-native.tree").read_text())
        row["same_rooted_topology"] = not different
        row["same_candidate_clades"] = cluster_set(a.output / f"{mode}-java.clusters") == cluster_set(
            a.output / f"{mode}-native.clusters")
        rows.append(row)
        (a.output / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
        assert not different and row["same_candidate_clades"], row
        assert row["java"]["score"] == row["native"]["score"], row
        assert row["java"]["weight_calls"] == row["native"]["weight_calls"], row
        print(mode, "PASS: rooted topology, candidates, score, every weight call", flush=True)


if __name__ == "__main__":
    main()
