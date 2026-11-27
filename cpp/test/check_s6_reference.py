#!/usr/bin/env python3
"""Sequential Java S6/native s3 parity on multicopy fixtures above 100 taxa.

Fresh Java classes may be supplied with --java-classes; otherwise this script
compiles the retained source into its artifact directory. Inputs and logs are
kept with the JSON evidence. This is a correctness test, not a timing study.
"""
import argparse
import json
from pathlib import Path
import random
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def balanced(labels):
    if len(labels) == 1:
        return labels[0]
    mid = len(labels) // 2
    return f"({balanced(labels[:mid])},{balanced(labels[mid:])})"


def fixture(n, families=12):
    rng = random.Random(71893 + n)
    names = [f"T{i:04d}" for i in range(n)]
    genes = [balanced(names)]
    for family in range(1, families):
        observed = [x for x in names if rng.random() > 0.14]
        # Preserve much of the common species signal while perturbing some tips.
        for _ in range(5):
            a, b = rng.sample(range(len(observed)), 2)
            observed[a], observed[b] = observed[b], observed[a]
        if family % 3:
            size = min(35, max(6, n // 5))
            scope = observed[:size]
            one = [x for x in scope if rng.random() > 0.2]
            two = [x for x in scope if rng.random() > 0.2]
            if len(one) < 4:
                one = scope[:4]
            if len(two) < 4:
                two = scope[-4:]
            duplication = f"({balanced(one)},{balanced(two)})"
            genes.append(f"({duplication},{balanced(observed[size:])})")
        else:
            genes.append(balanced(observed))
    return "\n".join(x + ";" for x in genes) + "\n"


def clades(text):
    tokens = re.findall(r"[(),;]|[^(),;\s]+", text)
    stack = []
    found = set()
    for token in tokens:
        if token == "(":
            stack.append([])
        elif token == ")":
            taxa = frozenset().union(*stack.pop())
            found.add(taxa)
            if stack:
                stack[-1].append(taxa)
        elif token not in (",", ";"):
            stack[-1].append(frozenset([token.split(":", 1)[0]]))
    return found


def run(command, folder, name):
    (folder / f"{name}.command.json").write_text(json.dumps(command, indent=2) + "\n")
    completed = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, timeout=240)
    (folder / f"{name}.log").write_text(completed.stdout)
    if completed.returncode:
        raise RuntimeError(f"{name} failed:\n{completed.stdout[-4000:]}")
    return completed.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--java-classes", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--taxa", type=int, nargs="+", default=[101, 137, 200])
    parser.add_argument("--rounds", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--cross-tree", action="store_true")
    parser.add_argument("--gpu", action="store_true", help="native GPU versus Java CPU")
    parser.add_argument("--seed", type=int, default=0xDEADBEEFCAFE)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    classes = args.java_classes or args.work / "java-classes"
    if args.java_classes is None:
        classes.mkdir(exist_ok=True)
        run(["javac", "-d", str(classes), *map(str, sorted((ROOT / "src").rglob("*.java")))],
            args.work, "javac")
    results = []
    for n in args.taxa:
        input_file = args.work / f"genes-{n}.trees"
        input_file.write_text(fixture(n))
        for rounds in args.rounds:
            for threads in args.threads:
                folder = args.work / f"n{n}-r{rounds}-t{threads}"
                folder.mkdir(exist_ok=True)
                options = ["--cpu", "--threads", str(threads), "--s6-rounds", str(rounds),
                           "--s6-seed", str(args.seed), "--input", str(input_file)]
                if args.cross_tree:
                    options.append("--cross-tree-splits")
                java_tree = folder / "java.tree"
                native_tree = folder / "native.tree"
                java_log = run(["java", "-Xmx2g", f"-Dstelarpro.home={ROOT}",
                                f"-Djava.library.path={ROOT / 'native'}", "-cp", str(classes),
                                "stelarx.Main", *options, "--search-space", "S6",
                                "--output", str(java_tree)], folder, "java")
                native_log = run([str(args.native), *options, "--search-space", "s3",
                                  *(["--gpu-strict"] if args.gpu else []),
                                  "--output", str(native_tree)], folder, "native")
                java_score = int(re.findall(r"Final triplet score = (\d+)", java_log)[-1])
                native_score = int(re.findall(r"Final triplet score = (\d+)", native_log)[-1])
                java_calls = list(map(int, re.findall(r"Weight table: (\d+) splits scored", java_log)))
                native_calls = list(map(int, re.findall(r"Weight calculation: (\d+) scored splits", native_log)))
                missing = clades(java_tree.read_text()) ^ clades(native_tree.read_text())
                row = dict(taxa=n, rounds=rounds, threads=threads, seed=args.seed, java_score=java_score,
                           native_score=native_score, java_weight_calls=java_calls,
                           native_weight_calls=native_calls, different_rooted_clades=len(missing))
                results.append(row)
                (args.work / "parity.json").write_text(json.dumps(results, indent=2) + "\n")
                assert java_score == native_score, row
                assert not missing, row
                assert java_calls == native_calls, row
                print(f"PASS {n} taxa / {rounds} rounds / {threads} threads: "
                      f"score {java_score}, every weight call, rooted topology", flush=True)


if __name__ == "__main__":
    main()
