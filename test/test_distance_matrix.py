#!/usr/bin/env python3
"""
test_distance_matrix.py
========================
20 randomized tests for the STELAR-Pro distance matrix computation.

For each test case:
  1. Generate n random taxa and k random binary gene trees (some incomplete).
  2. Compute the expected distance matrix by brute-force path-finding in Python.
  3. Run STELAR-Pro with --verify-distance-matrix and parse its output.
  4. Assert all distances match within tolerance.

Usage:
  python3 test/test_distance_matrix.py [--stelar-pro-root PATH] [--mode cpu|gpu]
  --mode gpu  will be supported once the GPU distance-matrix kernel is implemented.
"""

import argparse
import math
import os
import random
import re
import subprocess
import sys
import tempfile

# ── Tree representation ───────────────────────────────────────────────────────

class Node:
    __slots__ = ("left", "right", "taxon", "parent")
    def __init__(self, left=None, right=None, taxon=None):
        self.left   = left
        self.right  = right
        self.taxon  = taxon
        self.parent = None
    def is_leaf(self):
        return self.left is None

def random_binary_tree(taxa):
    """Return a random balanced binary Node tree over the given taxa list."""
    taxa = list(taxa)
    if len(taxa) == 1:
        return Node(taxon=taxa[0])
    random.shuffle(taxa)
    k = random.randint(1, len(taxa) - 1)
    left  = random_binary_tree(taxa[:k])
    right = random_binary_tree(taxa[k:])
    node = Node(left=left, right=right)
    left.parent  = node
    right.parent = node
    return node

def to_newick(node):
    if node.is_leaf():
        return node.taxon
    return f"({to_newick(node.left)},{to_newick(node.right)})"

# ── Brute-force distance matrix ───────────────────────────────────────────────

def leaf_depths(node, depth=0):
    """Returns dict: taxon -> depth (edges from root)."""
    if node.is_leaf():
        return {node.taxon: depth}
    result = {}
    result.update(leaf_depths(node.left,  depth + 1))
    result.update(leaf_depths(node.right, depth + 1))
    return result

def paths_from_root(node):
    """Returns dict: taxon -> list of Node objects from root (inclusive) to leaf."""
    if node.is_leaf():
        return {node.taxon: [node]}
    result = {}
    for child in (node.left, node.right):
        for taxon, path in paths_from_root(child).items():
            result[taxon] = [node] + path
    return result

def pairwise_distances(tree_root):
    """
    Returns dict: (a, b) -> topological distance for all ordered pairs in the tree.
    Uses LCA via longest common prefix of root-to-leaf paths.
    """
    paths = paths_from_root(tree_root)
    taxa  = list(paths.keys())
    dists = {}
    for i, a in enumerate(taxa):
        for b in taxa[i + 1:]:
            pa, pb = paths[a], paths[b]
            # LCA depth = index of last node common to both paths
            lca_depth = 0
            for j in range(min(len(pa), len(pb))):
                if pa[j] is pb[j]:
                    lca_depth = j
                else:
                    break
            d = (len(pa) - 1) + (len(pb) - 1) - 2 * lca_depth
            dists[(a, b)] = d
            dists[(b, a)] = d
    return dists

def brute_force_distance_matrix(trees_with_taxa, all_taxa):
    """
    trees_with_taxa: list of (Node root, set of taxa in this tree)
    all_taxa:        ordered list of all taxa (index = taxon id)
    Returns: n×n list-of-lists of average distances (inf if never co-occur).
    """
    n = len(all_taxa)
    idx = {t: i for i, t in enumerate(all_taxa)}
    dist_sum    = [[0.0] * n for _ in range(n)]
    cooccur     = [[0]   * n for _ in range(n)]

    for (root, _present) in trees_with_taxa:
        dists = pairwise_distances(root)
        for (a, b), d in dists.items():
            ia, ib = idx[a], idx[b]
            dist_sum[ia][ib]  += d
            cooccur[ia][ib]   += 1

    result = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                result[i][j] = 0.0
            elif cooccur[i][j] > 0:
                result[i][j] = dist_sum[i][j] / cooccur[i][j]
            else:
                result[i][j] = math.inf
    return result

# ── Run STELAR-Pro and parse output ─────────────────────────────────────────────

def parse_rooted_newick(text):
    """
    Parse one rooted binary Newick line (ASTRAL-Pro3 -T output: branch lengths
    and internal labels are ignored) into the Node class above.
    """
    text = re.sub(r"\[[^\]]*\]", "", text).strip().rstrip(";")
    stack = []
    token = ""
    i = 0
    while i < len(text):
        c = text[i]
        if c == "(":
            stack.append("(")
        elif c in ",)":
            if token.strip():
                stack.append(Node(taxon=token.strip().split(":")[0]))
                token = ""
            if c == ")":
                children = []
                while stack and stack[-1] != "(":
                    children.append(stack.pop())
                stack.pop()
                children.reverse()
                if len(children) != 2:
                    raise ValueError("rooted/tagged tree is not binary")
                node = Node(left=children[0], right=children[1])
                children[0].parent = node
                children[1].parent = node
                j = i + 1
                while j < len(text) and text[j] not in ",();":
                    j += 1          # skip the internal label / branch length
                i = j - 1
                stack.append(node)
        else:
            token += c
        i += 1
    if token.strip():
        stack.append(Node(taxon=token.strip().split(":")[0]))
    if len(stack) != 1:
        raise ValueError("malformed Newick")
    return stack[0]


def preprocessed_trees(newick_lines, stelarx_root):
    """
    The gene trees as STELAR-Pro's distance matrix sees them: every run first
    roots and tags its input with ASTRAL-Pro3 (stelar-pro -T), which replaces
    the rooting given in the file.  Returns [(Node root, set of taxa)].
    """
    build_dir = os.path.join(stelarx_root, "build")
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tre", delete=False) as f:
        f.write("\n".join(newick_lines) + "\n")
        input_path = f.name
    output_path = input_path + ".tagged"
    try:
        cmd = ["java", "-cp", build_dir, "stelarx.Main", "-T", "-q",
               "-i", input_path, "-o", output_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            raise RuntimeError(
                f"STELAR-Pro -T exited {result.returncode}\nstderr:\n{result.stderr[:2000]}")
        trees = []
        with open(output_path) as f:
            for line in f:
                if line.strip():
                    root = parse_rooted_newick(line)
                    trees.append((root, set(leaf_depths(root).keys())))
        return trees
    finally:
        os.unlink(input_path)
        if os.path.exists(output_path):
            os.unlink(output_path)


def run_stelar_pro_distance_matrix(newick_lines, stelarx_root, mode):
    """
    Write trees to a temp file, invoke STELAR-Pro --verify-distance-matrix,
    parse and return (taxa_order, n×n list-of-lists).
    """
    build_dir  = os.path.join(stelarx_root, "build")
    native_dir = os.path.join(stelarx_root, "native")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tre", delete=False) as f:
        f.write("\n".join(newick_lines) + "\n")
        input_path = f.name

    try:
        cmd = ["java",
               f"-Djava.library.path={native_dir}",
               "-cp", build_dir,
               "stelarx.Main",
               "-i", input_path,
               "--verify-distance-matrix", "-q"]
        if mode == "gpu":
            cmd.append("--gpu-strict")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            raise RuntimeError(
                f"STELAR-Pro exited {result.returncode}\nstderr:\n{result.stderr[:2000]}")
        return parse_dm_output(result.stdout)
    finally:
        os.unlink(input_path)

def parse_dm_output(output):
    """
    Parse the DISTANCE_MATRIX block from STELAR-Pro stdout.
    Returns (taxa_order: list[str], matrix: list[list[float]]).
    """
    lines = [l.strip() for l in output.splitlines() if l.strip()]
    taxa_order = []
    matrix     = []
    in_dm      = False
    n          = 0

    for line in lines:
        if line == "DISTANCE_MATRIX":
            in_dm = True
            continue
        if not in_dm:
            continue
        if line.startswith("n="):
            n = int(line[2:])
        elif line.startswith("taxa="):
            taxa_order = line[5:].split(",")
        elif line.startswith("row"):
            vals = line.split("=", 1)[1].split(",")
            row = []
            for v in vals:
                row.append(math.inf if v == "inf" else float(v))
            matrix.append(row)

    if not taxa_order or not matrix:
        raise RuntimeError(f"Failed to parse DISTANCE_MATRIX from output:\n{output[:1000]}")
    return taxa_order, matrix

# ── Test runner ───────────────────────────────────────────────────────────────

TOLERANCE = 1e-4

def run_test(seed, stelarx_root, mode, verbose=False):
    random.seed(seed)
    n_taxa    = random.randint(5, 50)
    k_trees   = random.randint(5, 20)
    miss_rate = random.uniform(0.0, 0.4)   # fraction of taxa to drop per tree

    all_taxa = [f"t{i}" for i in range(n_taxa)]

    trees_with_taxa = []
    newick_lines    = []

    for _ in range(k_trees):
        # Drop a random subset of taxa from this tree
        n_drop   = int(n_taxa * miss_rate)
        present  = all_taxa[:]
        if n_drop > 0 and len(present) - n_drop >= 2:
            drop  = random.sample(present, n_drop)
            present = [t for t in present if t not in drop]
        root = random_binary_tree(present)
        trees_with_taxa.append((root, set(present)))
        newick_lines.append(to_newick(root) + ";")

    # Python brute-force answer on the rooted/tagged trees STELAR-Pro actually
    # measures (its preprocessing re-roots every input tree; a root on a
    # different edge changes the number of edges between leaves).
    expected = brute_force_distance_matrix(
        preprocessed_trees(newick_lines, stelarx_root), all_taxa)

    # STELAR-Pro answer
    taxa_order, actual_raw = run_stelar_pro_distance_matrix(newick_lines, stelarx_root, mode)

    # Reorder actual to match all_taxa order (Java assigns IDs by first-occurrence)
    idx_java = {t: i for i, t in enumerate(taxa_order)}
    n = len(all_taxa)
    actual = [[0.0] * n for _ in range(n)]
    for i, ta in enumerate(all_taxa):
        for j, tb in enumerate(all_taxa):
            ji, jj = idx_java[ta], idx_java[tb]
            actual[i][j] = actual_raw[ji][jj]

    # Compare
    max_err = 0.0
    for i in range(n):
        for j in range(n):
            e = expected[i][j]
            a = actual[i][j]
            if math.isinf(e) and math.isinf(a):
                continue
            if math.isinf(e) or math.isinf(a):
                err = math.inf
            else:
                err = abs(e - a)
            if err > max_err:
                max_err = err

    passed = max_err <= TOLERANCE
    if verbose or not passed:
        status = "PASS" if passed else "FAIL"
        print(f"  seed={seed:5d}  n={n_taxa:3d}  k={k_trees:3d}  "
              f"miss_rate={miss_rate:.2f}  max_err={max_err:.2e}  [{status}]")
    return passed

def main():
    ap = argparse.ArgumentParser(description="Test STELAR-Pro distance matrix")
    ap.add_argument("--stelar-pro-root", default=".", help="Path to STELAR-Pro root")
    ap.add_argument("--mode", choices=["cpu", "gpu"], default="cpu",
                    help="Compute mode (gpu requires compiled GPU kernel)")
    ap.add_argument("--seeds", type=int, nargs="+",
                    help="Specific seeds to test (default: 1..20)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    seeds = args.seeds if args.seeds else list(range(1, 21))
    print(f"Running {len(seeds)} distance-matrix tests (mode={args.mode}) ...")

    passed = failed = 0
    for seed in seeds:
        try:
            ok = run_test(seed, args.stelar_pro_root, args.mode, verbose=args.verbose)
            if ok:
                passed += 1
            else:
                failed += 1
        except Exception as exc:
            print(f"  seed={seed}  ERROR: {exc}")
            failed += 1

    total = passed + failed
    print(f"\n{passed}/{total} tests passed", end="")
    if failed:
        print(f"  ({failed} FAILED)")
        sys.exit(1)
    else:
        print("  ✓")

if __name__ == "__main__":
    main()
