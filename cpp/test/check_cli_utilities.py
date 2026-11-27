#!/usr/bin/env python3
"""Utility/diagnostic flags, output safety and score-only CUDA regression checks."""
import argparse
import json
from pathlib import Path
import subprocess
from parity import score

ROOT = Path(__file__).resolve().parents[2]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native", type=Path, required=True)
    p.add_argument("--java-classes", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--gpu", action="store_true")
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    source = a.work / "genes.trees"
    original = "((A,B),((C,D),(E,F)));\n((A,C),((B,D),(E,F)));\n"
    source.write_text(original)
    native = [str(a.native.resolve()), "-i", str(source), "--threads", "3"]
    java = ["java", f"-Dstelarpro.home={ROOT}", "-cp", str(a.java_classes.resolve()),
            "stelarx.Main", "-i", str(source), "--cpu", "--threads", "3"]
    records = []

    def run(command, name, expected=0):
        result = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, timeout=120)
        (a.work / f"{name}.stdout").write_text(result.stdout)
        (a.work / f"{name}.stderr").write_text(result.stderr)
        assert (result.returncode == 0) == (expected == 0), (name, result.stderr)
        return result.stdout

    for flag in ("--verify-distance-matrix", "--verify-similarity-matrix", "--verify-upgma"):
        name = flag[2:]
        j = run(java + ["--search-space", "S1", flag], name + ".java")
        c = run(native + ["--cpu", "--search-space", "s1", flag], name + ".native")
        assert j == c, (name, j, c)
        records.append({"case": name, "status": "PASS", "java_stdout_equal": True})

    tagged_j, tagged_c = a.work / "java.tagged", a.work / "native.tagged"
    run(java + ["--tag-only", "-o", str(tagged_j)], "tag.java")
    run(native + ["--cpu", "--tag-only", "-o", str(tagged_c)], "tag.native")
    assert tagged_j.read_text() == tagged_c.read_text()
    records.append({"case": "tag-only", "status": "PASS", "java_output_equal": True})
    mapping = a.work / "mapping.txt"
    mapping.write_text("A A\nB B\nC C\nD D\nE E\nF F\n")
    run(java + ["--tag-only", "--gene-species-map", str(mapping), "-o", str(tagged_j)], "mapping.java")
    run(native + ["--cpu", "--tag-only", "--gene-species-map", str(mapping), "-o", str(tagged_c)],
        "mapping.native")
    assert tagged_j.read_text() == tagged_c.read_text()
    run(java + ["--gene-species-map", str(mapping)], "reject-map.java", expected=2)
    run(native + ["--cpu", "--gene-species-map", str(mapping)], "reject-map.native", expected=2)
    records.append({"case": "tag-only-map-and-inference-map-rejection", "status": "PASS"})
    legacy_j, legacy_c = a.work / "java.consensus", a.work / "native.consensus"
    run(java + ["--verify-greedy-consensus", "-o", str(legacy_j)], "legacy.java")
    run(native + ["--cpu", "--verify-greedy-consensus", "-o", str(legacy_c)], "legacy.native")
    assert legacy_j.read_text() == legacy_c.read_text()
    records.append({"case": "legacy-consensus-reference-verifier-bridge", "status": "PASS"})

    # Never truncate an input, including hard-link aliases and derived artifacts.
    alias = a.work / "input-alias"
    if not alias.exists():
        alias.hardlink_to(source)
    for target in (source, alias):
        run(native + ["--cpu", "-o", str(target)], "protect-" + target.name, expected=2)
        assert source.read_text() == original
    base = a.work / "derived"
    derived = Path(str(base) + ".stats.json")
    derived.write_text(original)
    run([str(a.native.resolve()), "-i", str(derived), "-o", str(base)], "protect-derived", expected=2)
    assert derived.read_text() == original
    records.append({"case": "input-and-derived-output-alias-protection", "status": "PASS"})

    for mode in ("s1", "s2", "s3"):
        out = a.work / f"{mode}.tree"
        run(native + ["--cpu", "--search-space", mode, "--verify-weights", "--log-file",
                      str(a.work / f"{mode}.log"), "-o", str(out)], mode + ".weights")
        assert Path(str(out) + ".weights.tsv").stat().st_size > 0
        assert "Final triplet score" in (a.work / f"{mode}.log").read_text()
        records.append({"case": mode + "-weight-and-log-dumps", "status": "PASS"})
    species = a.work / "species.tree"
    species.write_text("((A,B),((C,D),(E,F)));\n")
    j = run(java + ["-c", str(species)], "score.java")
    c = run(native + ["--gpu-strict" if a.gpu else "--cpu", "-c", str(species)], "score.native")
    assert score(j) == score(c)
    records.append({"case": "score-only-gpu" if a.gpu else "score-only-cpu", "status": "PASS"})
    (a.work / "results.json").write_text(json.dumps(records, indent=2) + "\n")
    print(f"PASS: {len(records)} CLI utility/diagnostic groups; evidence {a.work}")


if __name__ == "__main__":
    main()
