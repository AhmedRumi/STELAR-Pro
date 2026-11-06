#!/usr/bin/env python3
"""Independent rooted-triplet oracle for the built-in STELAR-Pro scorer.

The scorer roots and tags its input with ASTRAL-Pro3 before counting
triples; the oracle enumerates triples on that same rooted/tagged output
(see test_stelar_pro_differential.py for the observation rule).
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys
import tempfile

from test_stelar_pro_differential import oracle_score, parse_newick, rooted_tagged_genes

ROOT = pathlib.Path(__file__).resolve().parents[1]
GENES = ROOT / "test/input/test_5taxa.tre"
CANDIDATE = ROOT / "test/input/stelar_candidate_5taxa.tre"
UNROOTED = ROOT / "test/input/tc9_unrooted_simple.tre"
POLYTOMY_GENES = ROOT / "test/input/stelar_polytomy_5taxa.tre"
INCOMPLETE_POLYTOMY_GENES = ROOT / "test/input/stelar_polytomy_incomplete_6taxa.tre"
SIX_TAXON_CANDIDATE = ROOT / "test/input/stelar_candidate_6taxa.tre"


def stelar_score(genes_path: pathlib.Path = GENES,
                 candidate_path: pathlib.Path = CANDIDATE) -> int:
    command = [
        "java", "-cp", str(ROOT / "build"), "stelarx.Main",
        "-i", str(genes_path), "--score-species-tree", str(candidate_path),
        "--cpu", "-q",
    ]
    run = subprocess.run(command, cwd=ROOT, text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if run.returncode:
        raise AssertionError(f"built-in scorer failed:\n{run.stdout}")
    import re
    match = re.search(r"^TRIPLET_SCORE:\s*(\d+)\s*$", run.stdout, re.MULTILINE)
    if not match:
        raise AssertionError(f"built-in scorer emitted no triplet score:\n{run.stdout}")
    return int(match.group(1))


def check(genes_path: pathlib.Path, candidate_path: pathlib.Path,
          work: pathlib.Path) -> int:
    candidate = parse_newick(candidate_path.read_text())
    expected = oracle_score(candidate, rooted_tagged_genes(genes_path, work))
    score = stelar_score(genes_path, candidate_path)
    assert score == expected, (genes_path.name, score, expected)
    return score


def main() -> int:
    if os.environ.get("STELAR_PRO_SKIP_BUILD") != "1":
        subprocess.run([str(ROOT / "build.sh")], cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL)
    with tempfile.TemporaryDirectory(prefix="stelar-pro-triplets-") as directory:
        work = pathlib.Path(directory)
        score = check(GENES, CANDIDATE, work)
        polytomy_score = check(POLYTOMY_GENES, CANDIDATE, work)
        incomplete_polytomy_score = check(INCOMPLETE_POLYTOMY_GENES, SIX_TAXON_CANDIDATE, work)

        # Unrooted input is rooted by the same preprocessing, so it scores
        # exactly like its rooted/tagged form.  (The species tree may only
        # contain taxa present in the gene trees, so use a 4-taxon candidate.)
        unrooted_candidate = work / "unrooted-candidate.tre"
        unrooted_candidate.write_text("((A,B),(C,D));\n")
        unrooted_score = check(UNROOTED, unrooted_candidate, work)

    print(f"PASS: binary={score}; polytomy={polytomy_score}; "
          f"incomplete-polytomy={incomplete_polytomy_score}; unrooted={unrooted_score}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
