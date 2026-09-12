#!/usr/bin/env python3
"""GDL normalization tests, independent of native executables."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from calculate_rf_rate import compare_newicks, profile_newick


class GDLInputTests(unittest.TestCase):
    def prepare(self, work, text, expected=2, reference="((0,1),(2,3));\n"):
        source = work / "all_gt.trees"
        source.write_text(text)
        truth = work / "s_tree.trees"
        truth.write_text(reference)
        result = subprocess.run([
            sys.executable, str(ROOT / "scripts/prepare_gdl_input.py"),
            "--input", str(source), "--reference", str(truth),
            "--output", str(work / "normalized.trees"),
            "--metadata", str(work / "metadata.json"), "--expected-genes", str(expected),
        ], capture_output=True, text=True)
        self.assertEqual(source.read_text(), text)
        return result

    def test_copies_lengths_comments_internal_labels_and_multiline(self):
        text = "[note ; , (] (('0_0_0':0.1,0_1_0:2.0)dup,(1_0_0:3,\n2_0_0:4)99,3_0_0:5);\n((0,1),(2,3));\n"
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            result = self.prepare(work, text)
            self.assertEqual(result.returncode, 0, result.stderr)
            normalized = (work / "normalized.trees").read_text()
            expected = text.replace("'0_0_0'", "'0'").replace("0_1_0", "0").replace("1_0_0", "1").replace("2_0_0", "2").replace("3_0_0", "3").replace("3,\n2", "3, 2")
            self.assertEqual(normalized, expected)
            metadata = json.loads((work / "metadata.json").read_text())
            self.assertEqual(metadata["gene_trees"], 2)
            self.assertEqual(metadata["renamed_leaves"], 5)
            self.assertEqual(metadata["observed_species"], 4)
            self.assertEqual(metadata["input_sha256"], hashlib.sha256(text.encode()).hexdigest())
            self.assertEqual(metadata["normalized_sha256"], hashlib.sha256(normalized.encode()).hexdigest())

    def test_reject_unknown_or_ambiguous_labels(self):
        for label in ("9_0_0", "0_extra_copy", "0_0", "other"):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as directory:
                result = self.prepare(Path(directory), f"(({label},1_0_0),(2_0_0,3_0_0));", expected=1)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("cannot map leaf", result.stderr)

    def test_reject_partial_or_malformed_input(self):
        for text in ("((0,1),(2,3));", "((0,1),(2,3);", "((0,1),(2,3))"):
            with self.subTest(text=text), tempfile.TemporaryDirectory() as directory:
                self.assertNotEqual(self.prepare(Path(directory), text).returncode, 0)

    def test_multiple_trees_same_line(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.prepare(Path(directory), "((0,1),(2,3)); ((0,2),(1,3));")
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_loss_does_not_invent_missing_species(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            result = self.prepare(work, "((0_0_0,1_0_0),2_0_0);", expected=1)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(json.loads((work / "metadata.json").read_text())["missing_species"], ["3"])
            with self.assertRaisesRegex(ValueError, "different taxa"):
                compare_newicks("((0,1),2);", "((0,1),(2,3));")

    def test_deep_species_tree(self):
        tree = "0"
        for species in range(1, 1500):
            tree = f"({tree},{species})"
        self.assertEqual(len(profile_newick(tree + ";").taxa), 1500)


if __name__ == "__main__":
    unittest.main()
