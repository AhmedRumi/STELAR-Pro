#!/usr/bin/env python3
"""Independent rooted-triplet oracle for the built-in STELAR-Pro scorer."""

from __future__ import annotations

import itertools
import os
import pathlib
import re
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
GENES = ROOT / "test/input/test_5taxa.tre"
CANDIDATE = ROOT / "test/input/stelar_candidate_5taxa.tre"
UNROOTED = ROOT / "test/input/tc9_unrooted_simple.tre"
POLYTOMY_GENES = ROOT / "test/input/stelar_polytomy_5taxa.tre"
INCOMPLETE_POLYTOMY_GENES = ROOT / "test/input/stelar_polytomy_incomplete_6taxa.tre"
SIX_TAXON_CANDIDATE = ROOT / "test/input/stelar_candidate_6taxa.tre"


class Node:
    def __init__(self, children=None, name=None):
        self.children = children or []
        self.name = name


def parse_newick(text: str) -> Node:
    stack: list[Node | str] = []
    token = ""
    for char in text.strip().rstrip(";"):
        if char == "(":
            stack.append(char)
        elif char in ",)":
            if token.strip():
                stack.append(Node(name=token.strip()))
                token = ""
            if char == ")":
                children = []
                while stack and stack[-1] != "(":
                    children.append(stack.pop())
                if not stack:
                    raise ValueError("unbalanced Newick")
                stack.pop()
                stack.append(Node(children=list(reversed(children))))
        else:
            token += char
    if len(stack) != 1 or not isinstance(stack[0], Node):
        raise ValueError("malformed Newick")
    return stack[0]


def leaves(node: Node) -> set[str]:
    if node.name is not None:
        return {node.name}
    return set().union(*(leaves(child) for child in node.children))


def rooted_pair(node: Node, triple: tuple[str, str, str]):
    if node.name is not None:
        return None
    child_sets = [leaves(child) for child in node.children]
    locations = [next(i for i, group in enumerate(child_sets) if taxon in group)
                 for taxon in triple]
    if locations[0] == locations[1] != locations[2]:
        return tuple(sorted(triple[:2]))
    if locations[0] == locations[2] != locations[1]:
        return tuple(sorted((triple[0], triple[2])))
    if locations[1] == locations[2] != locations[0]:
        return tuple(sorted(triple[1:]))
    if locations[0] == locations[1] == locations[2]:
        return rooted_pair(node.children[locations[0]], triple)
    return None


def oracle_score(species: Node, genes: list[Node]) -> int:
    taxa = sorted(leaves(species))
    total = 0
    for gene in genes:
        present = leaves(gene)
        for triple in itertools.combinations(taxa, 3):
            if set(triple) <= present:
                species_pair = rooted_pair(species, triple)
                gene_pair = rooted_pair(gene, triple)
                total += gene_pair is not None and species_pair == gene_pair
    return total


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
    match = re.search(r"^TRIPLET_SCORE:\s*(\d+)\s*$", run.stdout, re.MULTILINE)
    if not match:
        raise AssertionError(f"built-in scorer emitted no triplet score:\n{run.stdout}")
    return int(match.group(1))


def main() -> int:
    if os.environ.get("STELAR_PRO_SKIP_BUILD") != "1":
        subprocess.run([str(ROOT / "build.sh")], cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL)
    genes = [parse_newick(line) for line in GENES.read_text().splitlines() if line.strip()]
    candidate = parse_newick(CANDIDATE.read_text())
    expected = oracle_score(candidate, genes)
    assert expected == 21, expected
    score = stelar_score()
    assert score == expected, (score, expected)

    polytomy_genes = [parse_newick(line) for line in POLYTOMY_GENES.read_text().splitlines()
                      if line.strip()]
    polytomy_expected = oracle_score(candidate, polytomy_genes)
    polytomy_score = stelar_score(POLYTOMY_GENES)
    assert polytomy_score == polytomy_expected, (polytomy_score, polytomy_expected)

    incomplete_polytomy_genes = [
        parse_newick(line) for line in INCOMPLETE_POLYTOMY_GENES.read_text().splitlines()
        if line.strip()
    ]
    six_taxon_candidate = parse_newick(SIX_TAXON_CANDIDATE.read_text())
    incomplete_polytomy_expected = oracle_score(
        six_taxon_candidate, incomplete_polytomy_genes)
    incomplete_polytomy_score = stelar_score(
        INCOMPLETE_POLYTOMY_GENES, SIX_TAXON_CANDIDATE)
    assert incomplete_polytomy_score == incomplete_polytomy_expected, (
        incomplete_polytomy_score, incomplete_polytomy_expected)

    with tempfile.TemporaryDirectory(prefix="stelar-pro-expected-crash-") as directory:
        reject_env = os.environ.copy()
        reject_env["STELAR_PRO_CRASH_DIR"] = str(pathlib.Path(directory) / "crash_logs")
        rejected = subprocess.run(
            ["java", "-cp", str(ROOT / "build"), "stelarx.Main",
             "-i", str(UNROOTED), "--cpu", "-q"],
            cwd=ROOT, env=reject_env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert rejected.returncode != 0
    assert "never roots input trees arbitrarily" in rejected.stdout

    print(f"PASS: binary={score}; polytomy={polytomy_score}; "
          f"incomplete-polytomy={incomplete_polytomy_score}; unrooted rejected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
