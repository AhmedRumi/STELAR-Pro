#!/usr/bin/env python3
"""Known-answer and end-to-end tests for the RF-rate scripts."""

from __future__ import annotations

import csv
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from calculate_rf_rate import (  # noqa: E402
    NewickError,
    compare_newicks,
    parse_newick,
    profile_newick,
)


def reference_splits(newick: str) -> set[frozenset[frozenset[str]]]:
    """Independent edge-removal split oracle for unique-leaf trees."""
    root = parse_newick(newick)
    adjacency: dict[int, set[int]] = {}
    labels: dict[int, str] = {}
    next_id = 0

    def build(node, parent: int | None = None) -> int:
        nonlocal next_id
        node_id = next_id
        next_id += 1
        adjacency[node_id] = set()
        if node.is_leaf:
            assert node.label is not None
            labels[node_id] = node.label
        if parent is not None:
            adjacency[node_id].add(parent)
            adjacency[parent].add(node_id)
        for child in node.children:
            build(child, node_id)
        return node_id

    build(root)
    taxa = frozenset(labels.values())
    splits: set[frozenset[frozenset[str]]] = set()
    for left, neighbors in adjacency.items():
        for right in neighbors:
            if left >= right:
                continue
            stack = [left]
            visited = {right}
            side: set[str] = set()
            while stack:
                node_id = stack.pop()
                if node_id in visited:
                    continue
                visited.add(node_id)
                if node_id in labels:
                    side.add(labels[node_id])
                stack.extend(adjacency[node_id] - visited)
            other = taxa - side
            if len(side) < 2 or len(other) < 2:
                continue
            side_frozen = frozenset(side)
            splits.add(frozenset((side_frozen, other)))
    return splits


def random_binary_newick(labels: list[str], rng: random.Random) -> str:
    components = labels[:]
    while len(components) > 1:
        left = components.pop(rng.randrange(len(components)))
        right = components.pop(rng.randrange(len(components)))
        components.append(f"({left},{right})")
    return components[0] + ";"


class RFRateTest(unittest.TestCase):
    def test_identical_quartet(self):
        result = compare_newicks("((A,B),(C,D));", "((A,B),(C,D));")
        self.assertEqual(result.raw_rf, 0)
        self.assertEqual(result.rf_rate, 0.0)

    def test_maximally_different_quartet(self):
        result = compare_newicks("((A,C),(B,D));", "((A,B),(C,D));")
        self.assertEqual((result.raw_rf, result.max_rf), (2, 2))
        self.assertEqual(result.rf_rate, 1.0)

    def test_one_of_two_five_taxon_splits_differs(self):
        result = compare_newicks("(((A,B),D),(C,E));", "(((A,B),C),(D,E));")
        self.assertEqual((result.shared_splits, result.raw_rf, result.max_rf), (1, 2, 4))
        self.assertEqual(result.rf_rate, 0.5)

    def test_unrooted_result_is_root_placement_invariant(self):
        result = compare_newicks("(A,(B,(C,D)));", "((A,B),(C,D));")
        self.assertEqual(result.rf_rate, 0.0)

    def test_labels_lengths_comments_and_internal_tags(self):
        output = "[&R](('A one':1,B:2)D:3,(C:1,D:1)90:1);"
        truth = "(('A one':9,B:9):9,(C:9,D:9):9);"
        result = compare_newicks(output, truth)
        self.assertEqual(result.rf_rate, 0.0)

    def test_duplicate_taxa_fail_in_both_species_trees(self):
        with self.assertRaisesRegex(ValueError, "output species tree.*duplicate"):
            compare_newicks("(((A,A),B),(C,D));", "((A,B),(C,D));")
        with self.assertRaisesRegex(ValueError, "true species tree.*duplicate"):
            compare_newicks("((A,B),(C,D));", "(((A,A),B),(C,D));")

    def test_mismatched_taxa_fail(self):
        with self.assertRaisesRegex(ValueError, "different taxa"):
            compare_newicks("((A,B),(C,X));", "((A,B),(C,D));")

    def test_malformed_newick_fails(self):
        with self.assertRaises(NewickError):
            parse_newick("((A,B),(C,D);")

    def test_random_unique_trees_match_edge_removal_oracle(self):
        rng = random.Random(29091)
        for taxa_count in range(4, 15):
            labels = [f"T{index}" for index in range(taxa_count)]
            for _ in range(10):
                newick = random_binary_newick(labels, rng)
                profile = profile_newick(newick)
                taxa = frozenset(profile.taxa)
                decoded = set()
                for mask in profile.splits:
                    side = frozenset(
                        profile.taxa[index]
                        for index in range(taxa_count)
                        if mask & (1 << index)
                    )
                    decoded.add(frozenset((side, taxa - side)))
                self.assertEqual(decoded, reference_splits(newick))

    def test_dataset_cli(self):
        with tempfile.TemporaryDirectory() as temporary:
            dataset = Path(temporary) / "dataset"
            for replicate in ("R1", "R2"):
                (dataset / replicate).mkdir(parents=True)
            truth = "((A,B),(C,D));\n"
            (dataset / "R1" / "s_tree.trees").write_text(truth)
            (dataset / "R2" / "s_tree.trees").write_text(truth)
            (dataset / "R1" / "out-stelar-pro.tre").write_text(
                "((A,B),(C,D));\n"
            )
            (dataset / "R2" / "out-stelar-pro.tre").write_text(
                "((A,C),(B,D));\n"
            )
            output = Path(temporary) / "results"
            command = [
                sys.executable,
                str(ROOT / "scripts" / "calculate_dataset_rf_rates.py"),
                str(dataset),
                "--results-dir",
                str(output),
            ]
            completed = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (output / "rf-rates-summary.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            self.assertEqual(int(rows[0]["replicates"]), 2)
            self.assertAlmostEqual(float(rows[0]["mean_rf_rate"]), 0.5)
            with (output / "rf-rates.csv").open() as handle:
                rate_rows = list(csv.DictReader(handle))
            self.assertEqual([row["replicate"] for row in rate_rows], ["R1", "R2"])
            self.assertEqual([float(row["rf_rate"]) for row in rate_rows], [0.0, 1.0])


if __name__ == "__main__":
    unittest.main()
