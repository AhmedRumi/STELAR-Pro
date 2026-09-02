#!/usr/bin/env python3
"""End-to-end inference invariants with independent final-score validation."""

from __future__ import annotations

import argparse
import os
import pathlib
import re
import subprocess
import tempfile

from test_stelar_pro_differential import annotate, oracle_score, parse_newick

ROOT = pathlib.Path(__file__).resolve().parents[1]
INCOMPLETE = ROOT / "test/input/test_incomplete.tre"
POLYTOMY = ROOT / "test/input/stelar_polytomy_incomplete_6taxa.tre"


def read_trees(path: pathlib.Path):
    return [parse_newick(line) for line in path.read_text().splitlines() if line.strip()]


def run_inference(genes_path: pathlib.Path, output_path: pathlib.Path,
                  compute: str,
                  extra: tuple[str, ...] = ()) -> tuple[str, int]:
    command = [
        "java", "-Djava.library.path=" + str(ROOT / "native"),
        "-cp", str(ROOT / "build"), "stelarx.Main", compute, "-q",
        "-i", str(genes_path), "-o", str(output_path),
        *extra,
    ]
    run = subprocess.run(command, cwd=ROOT, text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if run.returncode:
        raise AssertionError(f"inference failed ({extra}):\n{run.stdout}")
    if not output_path.is_file() or not output_path.read_text().strip():
        raise AssertionError("inference emitted no species tree")
    if compute == "--gpu-strict" and "[STELAR-Pro GPU] weight" not in run.stdout:
        raise AssertionError("strict CUDA inference used no GPU weight phase")
    matches = re.findall(r"(?:Final triplet score\s*=|Triplet score\s+)([0-9]+)", run.stdout)
    if not matches:
        raise AssertionError(f"inference emitted no final triplet score:\n{run.stdout}")
    return output_path.read_text().strip(), int(matches[-1])


def validate_run(genes_path: pathlib.Path, tree_text: str, reported: int) -> int:
    species = parse_newick(tree_text)
    annotate(species)
    if species.name is not None or len(species.children) != 2:
        raise AssertionError("inferred species tree does not retain an explicit binary root")
    genes = read_trees(genes_path)
    universe = frozenset().union(*(gene.taxa for gene in genes))
    if species.taxa != universe:
        raise AssertionError(f"species taxa {species.taxa} != gene universe {universe}")
    expected = oracle_score(species, genes)
    if expected != reported:
        raise AssertionError(f"independent final score={expected}, reported={reported}")
    return expected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpu", action="store_true", help="require strict CUDA")
    parser.add_argument("--no-build", action="store_true")
    args = parser.parse_args()
    if not args.no_build:
        subprocess.run([str(ROOT / "build.sh")], cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL)
    compute = "--gpu-strict" if args.gpu else "--cpu"
    executions = 0

    with tempfile.TemporaryDirectory(prefix="stelar-pro-inference-") as directory:
        work = pathlib.Path(directory)

        # Omitting search/scoring selectors is intentional: this is the product's
        # built-in S1 inference path, and its objective is independently enumerated.
        baseline_path = work / "incomplete-default.tre"
        baseline_tree, baseline_score = run_inference(
            INCOMPLETE, baseline_path, compute)
        validate_run(INCOMPLETE, baseline_tree, baseline_score)
        executions += 1

        # Native unresolved scoring and default deterministic preprocessing both
        # need end-to-end coverage on incomplete internal polytomies.
        for keep in (False, True):
            extra = ("--keep-polytomy-during-inference",) if keep else ()
            output = work / f"poly-{'keep' if keep else 'refine'}.tre"
            tree, score = run_inference(POLYTOMY, output, compute, extra)
            validate_run(POLYTOMY, tree, score)
            executions += 1

        # Threading, hash seed count, and reachability pruning are implementation
        # choices and must not change this deterministic fixture's result.
        if not args.gpu:
            variants = (
                ("threads-1", ("--threads", "1")),
                ("threads-4", ("--threads", str(min(4, os.cpu_count() or 1)))),
                ("seed-1", ("--seeds", "1")),
                ("seed-4", ("--seeds", "4")),
                ("unpruned", ("--no-prune-search-space",)),
            )
            for label, extra in variants:
                tree, score = run_inference(
                    INCOMPLETE, work / f"determinism-{label}.tre",
                    compute, extra)
                validate_run(INCOMPLETE, tree, score)
                if tree != baseline_tree or score != baseline_score:
                    raise AssertionError(f"determinism variant changed result: {label}")
                executions += 1

    print(f"STELAR-Pro inference invariants: PASS ({executions} end-to-end runs, "
          f"compute={'CUDA-strict' if args.gpu else 'CPU'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
