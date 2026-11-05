#!/usr/bin/env python3
"""Deterministic randomized differential tests for STELAR-Pro fixed-tree scoring.

The expected value is computed by direct rooted-triple enumeration.  This test
does not reuse STELAR-Pro clusters, partitions, hashes, or weight formulas.

What STELAR-Pro scores.  Every run first roots and tags the input gene trees
with ASTRAL-Pro3 (``astral-pro3 -T``); the rooting given in the input file is
not kept, and ASTRAL-Pro3 treats every polytomy as a duplication event.  A
rooted triple then counts as an observation only when the gene-tree node that
resolves it is a speciation node (not tagged ``D``).  The oracle therefore
reads the same ``astral-pro3 -T`` output the scorer consumes and applies that
rule; everything downstream of the rooting (the triple enumeration and the
agreement test) stays independent of the Java implementation.
"""

from __future__ import annotations

import argparse
import itertools
import os
import pathlib
import random
import re
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]


class Node:
    def __init__(self, name: str | None = None, children: list["Node"] | None = None,
                 duplication: bool = False):
        self.name = name
        self.children = children or []
        self.duplication = duplication
        self.taxa: frozenset[str] = frozenset()


def leaf(name: str) -> Node:
    return Node(name=name)


def annotate(node: Node) -> frozenset[str]:
    if node.name is not None:
        node.taxa = frozenset((node.name,))
    else:
        node.taxa = frozenset().union(*(annotate(child) for child in node.children))
    return node.taxa


def newick(node: Node, mirror: bool = False) -> str:
    if node.name is not None:
        return node.name
    children = list(reversed(node.children)) if mirror else node.children
    return "(" + ",".join(newick(child, mirror) for child in children) + ")"


def parse_newick(text: str) -> Node:
    """Small independent parser.

    Accepts the generated undecorated subset and ASTRAL-Pro3's tagged output:
    branch lengths are ignored, bracket comments are ignored, and an internal
    node label ``D`` marks a duplication node.
    """
    text = re.sub(r"\[[^\]]*\]", "", text).strip().rstrip(";")
    stack: list[Node | str] = []
    token = ""
    i = 0
    while i < len(text):
        char = text[i]
        if char == "(":
            stack.append(char)
        elif char in ",)":
            if token.strip():
                stack.append(leaf(token.strip().split(":")[0]))
                token = ""
            if char == ")":
                children: list[Node] = []
                while stack and stack[-1] != "(":
                    item = stack.pop()
                    if not isinstance(item, Node):
                        raise ValueError("malformed Newick")
                    children.append(item)
                if not stack:
                    raise ValueError("unbalanced Newick")
                stack.pop()
                label = ""
                j = i + 1
                while j < len(text) and text[j] not in ",();":
                    label += text[j]
                    j += 1
                i = j - 1
                stack.append(Node(children=list(reversed(children)),
                                  duplication=label.split(":")[0] == "D"))
        else:
            token += char
        i += 1
    if token.strip():
        stack.append(leaf(token.strip().split(":")[0]))
    if len(stack) != 1 or not isinstance(stack[0], Node):
        raise ValueError("malformed Newick")
    annotate(stack[0])
    return stack[0]


def resolved_pair(node: Node, triple: tuple[str, str, str]) -> tuple[str, str] | None:
    """The pair of the triple that is together below the other, or None."""
    if node.name is not None:
        return None
    locations = [next(i for i, child in enumerate(node.children) if taxon in child.taxa)
                 for taxon in triple]
    if locations[0] == locations[1] == locations[2]:
        return resolved_pair(node.children[locations[0]], triple)
    if locations[0] == locations[1] != locations[2]:
        return tuple(sorted((triple[0], triple[1])))
    if locations[0] == locations[2] != locations[1]:
        return tuple(sorted((triple[0], triple[2])))
    if locations[1] == locations[2] != locations[0]:
        return tuple(sorted((triple[1], triple[2])))
    return None


def observed_pair(node: Node, triple: tuple[str, str, str]) -> tuple[str, str] | None:
    """Like resolved_pair, but a triple resolved at a duplication node is not an observation."""
    if node.name is not None:
        return None
    locations = [next(i for i, child in enumerate(node.children) if taxon in child.taxa)
                 for taxon in triple]
    if locations[0] == locations[1] == locations[2]:
        return observed_pair(node.children[locations[0]], triple)
    if node.duplication:
        return None
    return resolved_pair(node, triple)


def oracle_score(species: Node, genes: list[Node]) -> int:
    total = 0
    for gene in genes:
        for triple in itertools.combinations(sorted(species.taxa & gene.taxa), 3):
            observed = observed_pair(gene, triple)
            total += observed is not None and resolved_pair(species, triple) == observed
    return total


def astral_pro_executable() -> pathlib.Path:
    override = os.environ.get("STELAR_PRO_EXECUTABLE")
    return pathlib.Path(override) if override else ROOT / "ASTER-Linux/bin/astral-pro3"


def rooted_tagged_genes(genes_path: pathlib.Path, work: pathlib.Path) -> list[Node]:
    """The gene trees exactly as the scorer sees them: ASTRAL-Pro3 -T output."""
    tagged = work / (genes_path.name + ".tagged")
    if not tagged.exists():
        run = subprocess.run([str(astral_pro_executable()), "-T", "-o", str(tagged),
                              str(genes_path)], text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if run.returncode or not tagged.exists():
            raise AssertionError(f"astral-pro3 -T failed for {genes_path}:\n{run.stdout}")
    return [parse_newick(line) for line in tagged.read_text().splitlines() if line.strip()]


def read_genes(genes_path: pathlib.Path) -> list[Node]:
    return [parse_newick(line) for line in genes_path.read_text().splitlines() if line.strip()]


def random_subtree(names: list[str], rng: random.Random, poly_probability: float) -> Node:
    if len(names) == 1:
        return leaf(names[0])
    shuffled = names[:]
    rng.shuffle(shuffled)
    if len(names) >= 3 and rng.random() < poly_probability:
        arity = rng.randint(3, min(5, len(names)))
    else:
        arity = 2
    cuts = sorted(rng.sample(range(1, len(names)), arity - 1))
    groups: list[list[str]] = []
    start = 0
    for end in cuts + [len(names)]:
        groups.append(shuffled[start:end])
        start = end
    return Node(children=[random_subtree(group, rng, poly_probability) for group in groups])


def random_rooted_tree(names: list[str], rng: random.Random,
                       poly_probability: float = 0.0) -> Node:
    """Build an explicitly rooted tree: the supplied top-level root is binary."""
    shuffled = names[:]
    rng.shuffle(shuffled)
    cut = rng.randint(1, len(shuffled) - 1)
    root = Node(children=[random_subtree(shuffled[:cut], rng, poly_probability),
                          random_subtree(shuffled[cut:], rng, poly_probability)])
    annotate(root)
    return root


class Case:
    def __init__(self, name: str, species: Node, genes: list[Node]):
        self.name = name
        self.species = species
        self.genes = genes
        self.expected: int | None = None   # oracle on the rooted/tagged genes, set by prepare()

    def prepare(self, work: pathlib.Path) -> None:
        genes_path = work / f"{self.name}.genes.tre"
        species_path = work / f"{self.name}.species.tre"
        genes_path.write_text("".join(newick(gene) + ";\n" for gene in self.genes))
        species_path.write_text(newick(self.species) + ";\n")
        self.expected = oracle_score(self.species, rooted_tagged_genes(genes_path, work))


def fixed_cases() -> list[Case]:
    raw = [
        ("minimal-three", "(A,(B,C));", ["(A,(B,C));", "(B,(A,C));"]),
        ("balanced-duplicates", "((A,B),(C,D));",
         ["((A,B),(C,D));", "((A,B),(C,D));", "((A,C),(B,D));"]),
        ("complete-internal-polytomy", "((A,B),((C,D),(E,F)));",
         ["((A,B,C),(D,(E,F)));", "((A,(B,C)),(D,E,F));",
          "(((A,B),C),(D,(E,F)));"]),
        ("incomplete-polytomies", "(((A,B),C),((D,E),(F,G)));",
         ["((A,B,C),(D,E));", "((A,C),(F,D,E));", "((B,C),(D,(F,G)));",
          "((A,B),(E,G));"]),
        ("deep-and-shallow", "(A,(B,(C,(D,(E,(F,G))))));",
         ["(A,(B,(C,(D,(E,(F,G))))));", "(((A,B,C),D),(E,F,G));",
          "((A,(C,E)),(B,D,F,G));"]),
        # These two internal polytomies have the same unordered set of four
        # groups, but a different group is outside the polytomous node.  A
        # rooted partition key must therefore keep the complement distinguished.
        ("polytomy-complement-dedup", "(((a1,b),(a2,c)),d);",
         ["(((a1,a2),b,c),d);", "(((a1,a2),b,d),c);"]),
    ]
    cases = []
    for name, species_text, gene_texts in raw:
        cases.append(Case(name, parse_newick(species_text),
                          [parse_newick(text) for text in gene_texts]))
    return cases


def generated_cases(count: int, seed: int) -> list[Case]:
    rng = random.Random(seed)
    cases: list[Case] = []
    for index in range(count):
        n = rng.randint(4, 10)
        taxa = [f"T{i}" for i in range(n)]
        species = random_rooted_tree(taxa, rng)
        genes = [random_rooted_tree(taxa, rng, rng.uniform(0.0, 0.75))]
        for _ in range(rng.randint(3, 10)):
            kept = [taxon for taxon in taxa if rng.random() >= rng.uniform(0.05, 0.45)]
            if len(kept) < 3:
                kept = rng.sample(taxa, 3)
            genes.append(random_rooted_tree(kept, rng, rng.uniform(0.0, 0.85)))
        if index % 4 == 0:
            genes.extend((genes[0], genes[-1]))  # frequency/deduplication coverage
        cases.append(Case(f"random-{index:02d}-n{n}-g{len(genes)}", species, genes))

    # Metamorphic variants of the first random case.  They are scorer-vs-oracle
    # cases only: ASTRAL-Pro3 roots each gene tree against a species tree it
    # infers from the whole input, so reversing, mirroring or duplicating the
    # input can legitimately change the chosen roots and hence the score.  Both
    # the scorer and the oracle see the same rooted/tagged trees, so agreement
    # is still exact.
    if cases:
        base = cases[0]
        cases.append(Case("metamorphic-reversed-lines", base.species, list(reversed(base.genes))))
        mirrored_species = parse_newick(newick(base.species, mirror=True) + ";")
        mirrored_genes = [parse_newick(newick(gene, mirror=True) + ";") for gene in base.genes]
        cases.append(Case("metamorphic-mirrored-children", mirrored_species, mirrored_genes))
        cases.append(Case("metamorphic-doubled-genes", base.species, base.genes + base.genes))
    return cases


def run_score(case: Case, work: pathlib.Path, compute: str,
              numeric: str = "long", threads: int | None = None) -> int:
    genes_path = work / f"{case.name}.genes.tre"
    species_path = work / f"{case.name}.species.tre"
    if case.expected is None:
        case.prepare(work)

    command = [
        "java", "-Djava.library.path=" + str(ROOT / "native"),
        "-cp", str(ROOT / "build"), "stelarx.Main", compute, "-q",
        "-i", str(genes_path), "--score-species-tree", str(species_path),
    ]
    env = os.environ.copy()
    if numeric != "long":
        env["STELAR_PRO_WEIGHT_FORCE_DOUBLE"] = "1"
        command += ["--large-n-score-type", numeric]
    if threads is not None:
        command += ["--threads", str(threads)]
    run = subprocess.run(command, cwd=ROOT, env=env, text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if run.returncode:
        raise AssertionError(f"{case.name}/{compute}/{numeric} failed:\n{run.stdout}")
    match = re.search(r"^TRIPLET_SCORE:\s*([0-9]+(?:[.]0+)?)\s*$",
                      run.stdout, re.MULTILINE)
    if not match:
        raise AssertionError(f"{case.name} emitted no score:\n{run.stdout}")
    if compute == "--gpu-strict" and "[STELAR-Pro GPU] weight" not in run.stdout:
        raise AssertionError(f"{case.name} did not execute a GPU weight phase")
    return int(float(match.group(1)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=12,
                        help="number of seeded random cases (plus fixed/metamorphic cases)")
    parser.add_argument("--seed", type=lambda value: int(value, 0), default=0x5E1A7)
    parser.add_argument("--gpu", action="store_true", help="require strict CUDA execution")
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--skip-numeric", action="store_true",
                        help="skip forced Int128/Double path checks")
    args = parser.parse_args()
    if args.cases < 0:
        parser.error("--cases must be non-negative")
    if not args.no_build:
        subprocess.run([str(ROOT / "build.sh")], cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL)

    cases = fixed_cases() + generated_cases(args.cases, args.seed)
    compute = "--gpu-strict" if args.gpu else "--cpu"
    assertions = 0
    with tempfile.TemporaryDirectory(prefix="stelar-pro-differential-") as directory:
        work = pathlib.Path(directory)
        for case in cases:
            case.prepare(work)
            actual = run_score(case, work, compute)
            if actual != case.expected:
                raise AssertionError(
                    f"{case.name}: oracle={case.expected}, actual={actual}")
            assertions += 1

        # Exercise serial and parallel CPU implementations on representative binary,
        # polytomous, and incomplete inputs. GPU kernels are thread-count independent.
        if not args.gpu:
            for case in (cases[0], cases[2], cases[3], cases[-1]):
                for threads in (1, min(4, os.cpu_count() or 1)):
                    actual = run_score(case, work, compute, threads=threads)
                    if actual != case.expected:
                        raise AssertionError(f"{case.name}: threads={threads}, score={actual}")
                    assertions += 1

        # Force the normally large-N-only arithmetic paths on small exact oracles.
        if not args.skip_numeric:
            for case in (cases[2], cases[3], cases[5]):
                for numeric in ("int128", "double"):
                    actual = run_score(case, work, compute, numeric=numeric)
                    if actual != case.expected:
                        raise AssertionError(
                            f"{case.name}/{numeric}: "
                            f"oracle={case.expected}, actual={actual}")
                    assertions += 1

    print(f"STELAR-Pro differential oracle: PASS ({len(cases)} cases, "
          f"{assertions} scorer executions, seed=0x{args.seed:x}, "
          f"compute={'CUDA-strict' if args.gpu else 'CPU'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
