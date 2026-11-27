#!/usr/bin/env python3
"""End-to-end parity for the retained native modes and unchanged Java sources.

Run from anywhere:
  python3 cpp/test/parity.py --native /path/to/stelar-pro-cpp --output /tmp/parity
Add --gpu /path/to/stelar-pro-cpp to replay every accepted inference on CUDA.
All input fixtures, commands, outputs and logs are retained in --output.
"""
import argparse
import json
import os
import pathlib
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]


def invoke(command, log, timeout=120):
    started = time.monotonic()
    with log.open("w") as stream:
        result = subprocess.run([str(word) for word in command], cwd=ROOT,
                                stdout=stream, stderr=subprocess.STDOUT,
                                timeout=timeout)
    return result.returncode, log.read_text(), time.monotonic() - started


def topology(text):
    # Inference emits topology-only Newick. Ignore child orientation while
    # preserving the root, and handle quoted species names without dependencies.
    tokens = re.findall(r"'(?:[^']|'')*'|[^\s(),;:]+|[(),;:]", text.strip())
    stack, clades, labels = [], set(), []
    for token in tokens:
        if token == "(":
            stack.append([])
        elif token == ")":
            children = stack.pop()
            united = frozenset().union(*children)
            clades.add(united)
            if stack:
                stack[-1].append(united)
        elif token not in (",", ";"):
            label = token[1:-1].replace("''", "'") if token.startswith("'") else token
            labels.append(label)
            if stack:
                stack[-1].append(frozenset((label,)))
    if len(labels) != len(set(labels)):
        raise AssertionError("inference output contains duplicate species")
    return frozenset(labels), frozenset(clades)


def score(log, kind="Final"):
    expression = (r"Final triplet score\s*=\s*([0-9]+)" if kind == "Final" else
                  r"optimization-objective triplet score\s*=\s*([0-9]+)")
    matches = re.findall(expression, log)
    if not matches and kind == "Final":
        matches = re.findall(r"TRIPLET_SCORE:\s*([0-9]+)", log)
    if not matches:
        raise AssertionError(f"missing {kind} triplet score")
    return int(matches[-1])


def cluster_set(path):
    result = set()
    for line in path.read_text().splitlines():
        if line.startswith("{"):
            result.add(frozenset(line[1:-1].split(",")))
        elif line and not line.startswith("cluster_id"):
            words = line.split("\t")
            result.add(frozenset(words[2].split(",")))
    return result


def evaluate(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    classes = output / "java-classes"
    compile_code, compile_log, _ = invoke(
        ["javac", "-sourcepath", "src", "-d", classes, "src/stelarx/Main.java"],
        output / "javac.log")
    if compile_code:
        raise RuntimeError("fresh Java reference compilation failed: " + compile_log)
    binary = args.native.resolve()
    backend = ROOT / "ASTER-Linux/bin/astral-pro3"
    java = ["java", f"-Dstelarpro.home={ROOT}", "-cp", classes, "stelarx.Main"]
    standard = ["--threads", "4", "--astral-pro-executable", backend]
    source = output / "inputs"
    source.mkdir(exist_ok=True)
    fixtures = {
        "binary": "((A,B),((C,D),(E,F)));\n((A,B),((C,D),(E,F)));\n((A,C),((B,D),(E,F)));\n",
        "incomplete": "((A,B),(C,D));\n((A,B),(E,F));\n((A,E),(C,(D,F)));\n((A,B),C);\n",
        "multicopy": (ROOT / "cpp/test/fixtures/s4_scopes.trees").read_text(),
        "polytomy": "(A,(B,C,D,E,F));\n((A,B),(C,D,E,F));\n((A,C),((B,D),(E,F)));\n",
    }
    for name, text in fixtures.items():
        (source / f"{name}.trees").write_text(text)
    taxa = source / "taxa.txt"
    taxa.write_text("# requested species, duplicate and absent entry\nA\nF\nD\nC\nA\nAbsent\n")
    species = source / "species.tree"
    species.write_text("((A,B),((C,D),(E,F)));\n")

    cases = [
        ("binary-s1", "binary", "s1", "S1", []),
        ("binary-s2", "binary", "s2", "S4", []),
        ("multicopy-s1", "multicopy", "s1", "S1", []),
        ("multicopy-s2", "multicopy", "s2", "S4", []),
        ("incomplete-s1", "incomplete", "s1", "S1", []),
        ("incomplete-s2", "incomplete", "s2", "S4", []),
        ("exhaustive", "multicopy", "s2", "S4", ["--exhaustive"]),
        ("cross-tree", "multicopy", "s2", "S4", ["--cross-tree-splits"]),
        ("cross-exhaustive", "multicopy", "s2", "S4", ["--exhaustive", "--cross-tree-splits"]),
        ("cross-alias", "multicopy", "s2", "S4", ["--search-mode", "full"]),
        ("cross-override-off", "multicopy", "s2", "S4", ["--cross-tree-splits", "--no-cross-tree-splits"]),
        ("s3-small-fallback", "multicopy", "s3", "S6", ["--s6-rounds", "2", "--s6-seed", "7"]),
        ("no-prune-s1", "multicopy", "s1", "S1", ["--no-prune-search-space"]),
        ("no-prune-s2", "multicopy", "s2", "S4", ["--no-prune-search-space"]),
        ("taxa-s1", "binary", "s1", "S1", ["--taxa-file", taxa]),
        ("taxa-s2", "binary", "s2", "S4", ["--taxa-file", taxa]),
        ("polytomy-s1", "polytomy", "s1", "S1", []),
        ("polytomy-s2", "polytomy", "s2", "S4", []),
        ("polytomy-keep", "polytomy", "s1", "S1", ["--keep-polytomy-during-inference"]),
        ("hash-seeds", "multicopy", "s2", "S4", ["--seeds", "4"]),
        ("completion-distance", "incomplete", "s2", "S4", ["--completion-method", "distance"]),
        ("score-type-int128", "binary", "s1", "S1", ["--large-n-score-type", "int128"]),
        ("score-type-double", "binary", "s2", "S4", ["--large-n-score-type", "double"]),
    ]
    records = []
    for name, fixture, native_mode, java_mode, extra in cases:
        path = output / name
        path.mkdir(exist_ok=True)
        inp = source / f"{fixture}.trees"
        java_tree = path / "java.tree"
        java_clusters = path / "java.clusters"
        cmd = java + ["-i", inp, "-o", java_tree, "--cpu", "--search-space", java_mode,
                      "--dump-clusters", java_clusters] + standard + extra
        java_code, java_log, java_time = invoke(cmd, path / "java.log")
        record = {"case": name, "java_seconds": java_time, "results": {}}
        for label, executable, compute in [("cpu", binary, "--cpu")] + (
                [("gpu", args.gpu.resolve(), "--gpu-strict")] if args.gpu else []):
            native_tree = path / f"{label}.tree"
            native_clusters = path / f"{label}.clusters"
            native_cmd = [executable, "-i", inp, "-o", native_tree, compute,
                          "--search-space", native_mode, "--dump-clusters", native_clusters] + standard + extra
            code, log, elapsed = invoke(native_cmd, path / f"{label}.log")
            try:
                if java_code or code:
                    raise AssertionError(f"exit Java={java_code}, native={code}")
                if topology(java_tree.read_text()) != topology(native_tree.read_text()):
                    raise AssertionError("rooted output topology differs")
                if score(java_log) != score(log):
                    raise AssertionError(f"final triplet score Java={score(java_log)}, native={score(log)}")
                if score(java_log, "Optimization") != score(log, "Optimization"):
                    raise AssertionError("optimization triplet score differs")
                if cluster_set(java_clusters) != cluster_set(native_clusters):
                    raise AssertionError("candidate clade sets differ")
                js = re.findall(r"Weight table:\s*(\d+) splits scored", java_log)
                stats = json.loads(pathlib.Path(str(native_tree) + ".stats.json").read_text())
                if js and int(js[-1]) != stats["scored_splits"]:
                    raise AssertionError(f"scored splits Java={js[-1]}, native={stats['scored_splits']}")
                result = {"status": "PASS", "score": score(log),
                          "scored_splits": stats["scored_splits"], "seconds": elapsed}
            except (AssertionError, FileNotFoundError, ValueError) as error:
                result = {"status": "FAIL", "error": str(error), "seconds": elapsed}
            record["results"][label] = result
            print(f"{name} {label}: {result}", flush=True)
        records.append(record)

    for name, fixture, extra in [
        ("score-binary", "binary", ["--score-species-tree", species]),
        ("score-polytomy", "polytomy", ["--score-species-tree", species]),
        ("score-taxa", "binary", ["--score-species-tree", species, "--taxa-file", taxa]),
    ]:
        path = output / name
        path.mkdir(exist_ok=True)
        inp = source / f"{fixture}.trees"
        jc, jl, _ = invoke(java + ["-i", inp, "--cpu"] + standard + extra, path / "java.log")
        nc, nl, _ = invoke([binary, "-i", inp, "--cpu"] + standard + extra, path / "cpu.log")
        ok = jc == nc == 0 and score(jl) == score(nl)
        records.append({"case": name, "results": {"cpu": {"status": "PASS" if ok else "FAIL",
                        "java_exit": jc, "native_exit": nc,
                        "java_score": score(jl) if jc == 0 else None,
                        "native_score": score(nl) if nc == 0 else None}}})
        print(name, records[-1]["results"], flush=True)

    for op in ("union", "intersection"):
        path = output / f"extract-{op}"
        path.mkdir(exist_ok=True)
        inp = source / "incomplete.trees"
        common = ["-i", inp, "--extract-taxa", "--taxa-set", op]
        jc, _, _ = invoke(java + common + ["-o", path / "java.txt"], path / "java.log")
        nc, _, _ = invoke([binary] + common + ["-o", path / "cpu.txt"], path / "cpu.log")
        ok = jc == nc == 0 and (path / "java.txt").read_text() == (path / "cpu.txt").read_text()
        records.append({"case": f"extract-{op}", "results": {"cpu": {"status": "PASS" if ok else "FAIL"}}})

    rejects = [
        ("reject-s1-cross", "s1", "S1", ["--cross-tree-splits"]),
        ("reject-exhaustive-s1", "s1", "S1", ["--exhaustive"]),
        ("reject-autocomplete", "s2", "S4", ["--autocomplete-incomplete-gene-trees"]),
        ("reject-anchor", "s2", "S4", ["--anchor-outgroup"]),
        ("reject-consensus", "s2", "S4", ["--consensus-experimental"]),
    ]
    for name, native_mode, java_mode, extra in rejects:
        path = output / name
        path.mkdir(exist_ok=True)
        common = ["-i", source / "binary.trees", "--cpu"] + standard
        # Search presets reset older enrichment options. Put the deliberately
        # unsupported configuration after the preset to exercise rejection.
        jc, _, _ = invoke(java + common + ["--search-space", java_mode] + extra, path / "java.log")
        nc, _, _ = invoke([binary] + common + ["--search-space", native_mode] + extra, path / "cpu.log")
        records.append({"case": name, "results": {"cpu": {"status": "PASS" if jc and nc else "FAIL",
                        "java_exit": jc, "native_exit": nc}}})

    (output / "results.json").write_text(json.dumps(records, indent=2) + "\n")
    failures = [(r["case"], engine, result) for r in records for engine, result in r["results"].items()
                if result["status"] != "PASS"]
    print(f"{len(records)} cases; {len(failures)} failures; evidence: {output / 'results.json'}")
    return bool(failures)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--native", type=pathlib.Path, required=True)
    parser.add_argument("--gpu", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    sys.exit(evaluate(parser.parse_args()))
