#!/usr/bin/env python3
"""Compare raw S4 construction stages with the unchanged Java implementation.

Usage: check_s4_reference.py /path/to/test_s4 /path/to/compiled/java/classes
The Java helper is compiled with:
  javac -sourcepath src -d /tmp/stelar-java-reference cpp/test/JavaS4Reference.java
"""
import argparse
import pathlib
import random
import subprocess
import tempfile


def canonical_lines(path):
    # Java and iostream spell equivalent hexadecimal doubles differently.
    lines = []
    for line in path.read_text().splitlines():
        if line.startswith("P "):
            words = line.split()
            words[-1] = float.fromhex(words[-1])
            lines.append(tuple(words))
        else:
            lines.append(line)
    return lines


def random_family(rng, names):
    leaves = [name for name in names if rng.random() < 0.7]
    leaves += rng.sample(names, rng.randrange(0, len(names) // 2 + 1))
    if len(leaves) < 4:
        leaves = names[:4]
    # Assign speciation versus duplication by actual species-set overlap.
    forest = [(leaf, {leaf}) for leaf in leaves]
    while len(forest) > 1:
        a = forest.pop(rng.randrange(len(forest)))
        b = forest.pop(rng.randrange(len(forest)))
        tag = "D" if a[1] & b[1] else "S"
        forest.append((f"({a[0]},{b[0]}){tag}", a[1] | b[1]))
    return forest[0][0] + ";"


def compare(binary, java_classes, input_path, work, name):
    baseline = None
    for threads in (1, 4):
        java = work / f"{name}.java.{threads}.txt"
        native = work / f"{name}.native.{threads}.txt"
        subprocess.run(["java", "-cp", str(java_classes), "JavaS4Reference",
                        str(input_path), str(java), str(threads)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([str(binary), str(input_path), str(native), str(threads)], check=True)
        j, c = canonical_lines(java), canonical_lines(native)
        if j != c:
            mismatch = next((i for i, pair in enumerate(zip(j, c)) if pair[0] != pair[1]),
                            min(len(j), len(c)))
            raise AssertionError(f"{name}, threads={threads}, line={mismatch + 1}: "
                                 f"Java={j[mismatch:mismatch+1]} C++={c[mismatch:mismatch+1]}")
        if baseline is not None and baseline != c:
            raise AssertionError("native S4 changes with thread count")
        baseline = c
    print(f"{name}: Java/C++ raw scopes, exact similarities, anchor ranks, UPGMA and completion PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=pathlib.Path)
    parser.add_argument("java_classes", type=pathlib.Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    java_classes = args.java_classes.resolve()
    with tempfile.TemporaryDirectory(prefix="stelar-cpp-s4-parity-") as temporary:
        work = pathlib.Path(temporary)
        compare(binary, java_classes, pathlib.Path(__file__).parent / "fixtures/s4_scopes.trees",
                work, "nested-scopes")
        rng = random.Random(137)
        for size in (8, 17, 70):
            names = [f"species{x}" for x in range(size)]
            families = [random_family(rng, names) for _ in range(12)]
            # Every registry species is observed and the >64 case crosses the
            # packed scope boundary regardless of random family sampling.
            scaffold = names[0]
            for name in names[1:]:
                scaffold = f"({scaffold},{name})S"
            families.insert(0, scaffold + ";")
            source = work / f"random{size}.trees"
            source.write_text("\n".join(families) + "\n")
            compare(binary, java_classes, source, work, f"random-{size}")


if __name__ == "__main__":
    main()
