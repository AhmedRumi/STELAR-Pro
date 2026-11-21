#!/usr/bin/env python3
"""Integration check of S6's 100-species switch in all three launchers."""

from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def balanced(labels):
    nodes = list(labels)
    while len(nodes) > 1:
        nodes = [f"({nodes[i]},{nodes[i + 1]})" if i + 1 < len(nodes)
                 else nodes[i] for i in range(0, len(nodes), 2)]
    return nodes[0] + ";\n"


def run(name, tree, outdir, search="S6", options=(), suffix=""):
    output = outdir / f"{name}-{search}{suffix}.tre"
    args = ["--cpu", "--threads", "2", "--s6-rounds", "1",
            "--input", str(tree), "--output", str(output), *options]
    if name == "production":
        command = ["java", "-Xmx2g", f"-Dstelarpro.home={ROOT}",
                   f"-Djava.library.path={ROOT / 'native'}", "-cp",
                   str(ROOT / "build"), "stelarx.Main", *args,
                   "--search-space", search]
    elif name == "without-subset-local":
        command = ["bash", str(ROOT / "temporary-builds" /
                   "s6-without-subset-local" / "run.sh"), *args, "--xmx", "2g"]
    else:
        command = ["bash", str(ROOT / "temporary-builds" /
                   "s6-subset-local-strict" / "run.sh"), *args, "--xmx", "2g"]
    result = subprocess.run(command, cwd=ROOT, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=180, check=False)
    if result.returncode:
        raise AssertionError(f"{name} {search} exited {result.returncode}:\n"
                             + result.stdout[-3000:])
    return output.read_text().strip(), result.stdout


def main():
    with tempfile.TemporaryDirectory(prefix="stelar-s6-small-") as work:
        folder = Path(work)
        for n in (5, 100):
            fixture = folder / f"gene-{n}.tre"
            fixture.write_text(balanced(f"T{i:03d}" for i in range(n)))
            reference, ref_log = run("production", fixture, folder, "S4")
            assert "Phase 1b S4" in ref_log
            for name in ("production", "without-subset-local",
                         "subset-local-strict"):
                result, log = run(name, fixture, folder)
                assert result == reference, (name, n, "S4 topology mismatch")
                assert f"S6: {n} species (<= 100); using S4" in log, (name, n)
                assert "Phase 1b S4" in log and "S6 candidates:" not in log, (name, n)
                if n == 5 and name == "production":
                    _, cross_log = run(name, fixture, folder,
                                       options=("--cross-tree-splits",), suffix="-cross")
                    assert "Phase 5b Cross-tree transitions" in cross_log
        fixture = folder / "gene-101.tre"
        fixture.write_text(balanced(f"T{i:03d}" for i in range(101)))
        for name in ("production", "without-subset-local", "subset-local-strict"):
            result, log = run(name, fixture, folder)
            assert result.endswith(";"), (name, "missing tree")
            assert "S6 candidates:" in log and "S6: 101 species" not in log, name
            expected_upgma = 1 if name == "without-subset-local" else 0
            assert f"{expected_upgma} UPGMA tree(s)" in log, (name, "wrong S6 variant")
            if name == "production":
                local_log = log
                _, cross_log = run(name, fixture, folder,
                                   options=("--cross-tree-splits",), suffix="-cross")
                assert "Phase 5b Cross-tree transitions" not in local_log
                assert "Phase 5b Cross-tree transitions" in cross_log
                candidate_lines = lambda value: [line.split("S6 candidates:", 1)[1]
                    for line in value.splitlines() if "S6 candidates:" in line]
                assert candidate_lines(local_log) == candidate_lines(cross_log)
        print("S6 small-species switch: PASS (5, 100, 101 species; all three variants)")


if __name__ == "__main__":
    main()
