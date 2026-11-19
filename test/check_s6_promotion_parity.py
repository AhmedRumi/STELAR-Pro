#!/usr/bin/env python3
"""Replay saved scope experiments against production S6 and an S4 control.

Inference is sequential. This is a correctness replay, not a timing benchmark.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from calculate_rf_rate import compare_newicks, read_tree


def signature(log):
    weights = re.findall(r"Weight table: (\d+) splits scored \[LONG\]", log)
    score = re.findall(r"Final triplet score = (\d+)", log)
    if not weights or not score:
        raise AssertionError("Missing exact weight/score result")
    return {"scored_splits_per_call": [int(n) for n in weights],
            "triplet_score": int(score[-1])}


def replay(archive, mode, rep, work):
    source = archive / rep / mode
    command = json.loads((source / "command.json").read_text())
    command = [arg for arg in command if arg != "--subsetlocal"
               and not arg.startswith("-Dstelarpro.experimental.scope=")]
    command[command.index("-cp") + 1] = str(REPO / "build")
    target = work / rep / mode
    target.mkdir(parents=True, exist_ok=True)
    tree = target / "species_tree.tre"
    command[command.index("--output") + 1] = str(tree)
    (target / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    with (target / "terminal.log").open("w") as stream:
        subprocess.run(command, cwd=REPO, stdout=stream,
                       stderr=subprocess.STDOUT, check=True, timeout=180)
    actual_log = (target / "terminal.log").read_text()
    expected = signature((source / "terminal.log").read_text())
    actual = signature(actual_log)
    if actual != expected:
        raise AssertionError(f"{mode} {rep}: {actual} != {expected}")
    rf = compare_newicks(read_tree(tree), read_tree(source / "species_tree.tre")).raw_rf
    if rf:
        raise AssertionError(f"{mode} {rep}: topology changed (RF {rf})")
    if mode.startswith("S6"):
        if "duplication-vertex scope" not in actual_log:
            raise AssertionError("Production S6 did not select duplication scope")
        if "shared S4 guide pool" in actual_log:
            raise AssertionError("Production S6 selected old global guide path")
    print(f"PASS {mode} {rep}: topology, score and every DP weight-call count match", flush=True)
    return {"replicate": rep, "reference": mode, "topology_rf": rf, **actual}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    rows = []
    for rep in (f"{i:02d}" for i in range(1, 11)):
        rows.append(replay(args.archive, "S6-subset-duplication", rep, args.work))
    rows.append(replay(args.archive, "S4-strict", "01", args.work))
    (args.work / "parity.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
