#!/usr/bin/env python3
"""Resume only the C++ s3 2,000-taxon R1 cell of the paired benchmark."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import benchmark_java_cpp as benchmark


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=benchmark.ROOT / "cpp/validation/paired-s2-s3-20261127")
    args = parser.parse_args()
    work = args.work.resolve()
    manifest = json.loads((work / "manifest.json").read_text())
    case = next(case for case in manifest["cases"] if case["dataset"] == "taxa2000_R1")
    assert benchmark.digest(Path(case["input"])) == case["input_sha256"]
    assert benchmark.digest(Path(case["truth"])) == case["truth_sha256"]
    native = benchmark.ROOT / "cpp/bin/stelar-pro-cpp"
    assert benchmark.digest(native) == manifest["cpp_binary_sha256"], "C++ executable differs from previous benchmark"
    runtime = SimpleNamespace(native=native, classes=work / "java-classes", work=work,
                              threads=manifest["threads"], xmx=manifest["xmx"],
                              timeout=7200, sample_interval=0.25)
    existing = json.loads((work / "results.json").read_text())
    row = benchmark.run_one(case, ("cpp-s3", "cpp", "s3"), runtime, manifest)
    rows = [old for old in existing if (old["dataset"], old["replicate"], old["method"])
            != (case["dataset"], case["replicate"], "cpp-s3")]
    rows.append(row)
    benchmark.report(manifest["cases"], rows, runtime)
    benchmark.save(work / "current.json", dict(status="partial_benchmark_complete",
                                                 completed_runs=len(rows),
                                                 large_cpp_s3_exit_code=row["exit_code"],
                                                 note="Fresh large Java S4/S6 reruns were stopped at user request; compare archived runs."))
    return row["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())
