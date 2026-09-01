#!/usr/bin/env python3
"""Calculate standard unrooted RF rate between inferred and true species trees."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


@dataclass(slots=True)
class Node:
    children: list["Node"]
    label: str | None = None

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickError(ValueError):
    """Raised when a Newick tree is malformed."""


class NewickParser:
    """Small topology-only parser supporting lengths, labels, and comments."""

    def __init__(self, text: str):
        self.text = text
        self.pos = 0

    def parse(self) -> Node:
        self._skip_ignored()
        if self.pos >= len(self.text):
            raise NewickError("tree is empty")
        root = self._subtree()
        self._skip_ignored()
        if self._peek() == ";":
            self.pos += 1
            self._skip_ignored()
        if self.pos != len(self.text):
            raise self._error("unexpected content after the tree")
        return root

    def _subtree(self) -> Node:
        self._skip_ignored()
        if self._peek() == "(":
            self.pos += 1
            children = [self._subtree()]
            while True:
                self._skip_ignored()
                token = self._peek()
                if token == ",":
                    self.pos += 1
                    children.append(self._subtree())
                elif token == ")":
                    self.pos += 1
                    break
                else:
                    raise self._error("expected ',' or ')'")
            if len(children) < 2:
                raise self._error("an internal node needs at least two children")
            label = self._label(required=False)
            self._branch_length()
            return Node(children, label)

        label = self._label(required=True)
        self._branch_length()
        return Node([], label)

    def _label(self, required: bool) -> str | None:
        self._skip_ignored()
        if self.pos >= len(self.text) or self._peek() in "(),:;":
            if required:
                raise self._error("expected a leaf label")
            return None

        if self._peek() == "'":
            self.pos += 1
            chars: list[str] = []
            while self.pos < len(self.text):
                char = self.text[self.pos]
                self.pos += 1
                if char == "'":
                    if self._peek() == "'":
                        chars.append("'")
                        self.pos += 1
                        continue
                    return "".join(chars)
                chars.append(char)
            raise self._error("unterminated quoted label")

        start = self.pos
        while self.pos < len(self.text):
            char = self.text[self.pos]
            if char.isspace() or char in "(),:;[]":
                break
            self.pos += 1
        if self.pos == start:
            if required:
                raise self._error("expected a leaf label")
            return None
        return self.text[start:self.pos]

    def _branch_length(self) -> None:
        self._skip_ignored()
        if self._peek() != ":":
            return
        self.pos += 1
        start = self.pos
        while self.pos < len(self.text):
            char = self.text[self.pos]
            if char in ",);":
                break
            if char == "[":
                self._comment()
            else:
                self.pos += 1
        if not self.text[start:self.pos].strip():
            raise self._error("branch length is empty")

    def _skip_ignored(self) -> None:
        while self.pos < len(self.text):
            if self.text[self.pos].isspace():
                self.pos += 1
            elif self.text[self.pos] == "[":
                self._comment()
            else:
                break

    def _comment(self) -> None:
        depth = 0
        while self.pos < len(self.text):
            char = self.text[self.pos]
            self.pos += 1
            if char == "[":
                depth += 1
            elif char == "]":
                depth -= 1
                if depth == 0:
                    return
        raise self._error("unterminated Newick comment")

    def _peek(self) -> str:
        return self.text[self.pos] if self.pos < len(self.text) else ""

    def _error(self, message: str) -> NewickError:
        return NewickError(f"{message} at character {self.pos + 1}")


@dataclass(frozen=True, slots=True)
class TreeProfile:
    taxa: tuple[str, ...]
    splits: frozenset[int]


@dataclass(frozen=True, slots=True)
class RFResult:
    taxa: int
    output_splits: int
    true_splits: int
    shared_splits: int
    output_only_splits: int
    true_only_splits: int
    raw_rf: int
    max_rf: int
    rf_rate: float


def parse_newick(newick: str) -> Node:
    return NewickParser(newick.strip()).parse()


def _leaf_labels(root: Node) -> list[str]:
    labels: list[str] = []
    stack = [root]
    while stack:
        node = stack.pop()
        if node.is_leaf:
            if node.label is None:
                raise NewickError("leaf has no label")
            labels.append(node.label)
        else:
            stack.extend(node.children)
    return labels


def profile_newick(
    newick: str,
    expected_taxa: tuple[str, ...] | None = None,
    tree_name: str = "tree",
) -> TreeProfile:
    """Map a single-copy tree's edges to unique nontrivial bipartitions."""
    root = parse_newick(newick)
    labels = _leaf_labels(root)
    if not labels:
        raise ValueError("tree contains no leaves")
    seen: set[str] = set()
    duplicate_labels: set[str] = set()
    for label in labels:
        if label in seen:
            duplicate_labels.add(label)
        seen.add(label)
    if duplicate_labels:
        raise ValueError(
            f"{tree_name} must contain each taxon exactly once; duplicate taxa: "
            + ", ".join(sorted(duplicate_labels)[:10])
        )

    observed = set(labels)
    taxa = tuple(sorted(observed)) if expected_taxa is None else expected_taxa
    expected = set(taxa)
    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        details = []
        if missing:
            details.append("missing=" + ",".join(missing[:10]))
        if extra:
            details.append("extra=" + ",".join(extra[:10]))
        raise ValueError(
            "output and true species trees have different taxa ("
            + "; ".join(details) + ")"
        )

    taxon_index = {label: index for index, label in enumerate(taxa)}
    all_mask = (1 << len(taxa)) - 1
    splits: set[int] = set()

    def visit(node: Node, is_root: bool = False) -> int:
        if node.is_leaf:
            assert node.label is not None
            return 1 << taxon_index[node.label]

        mask = 0
        for child in node.children:
            mask |= visit(child)

        if not is_root:
            side_taxa = mask.bit_count()
            if 2 <= side_taxa <= len(taxa) - 2:
                split_mask = mask
                complement = all_mask ^ split_mask
                if side_taxa > len(taxa) - side_taxa:
                    split_mask = complement
                elif side_taxa == len(taxa) - side_taxa:
                    split_mask = min(split_mask, complement)
                splits.add(split_mask)
        return mask

    visit(root, is_root=True)
    return TreeProfile(taxa=taxa, splits=frozenset(splits))


def compare_profiles(output: TreeProfile, truth: TreeProfile) -> RFResult:
    if output.taxa != truth.taxa:
        raise ValueError("output and true profiles use different taxon orders")
    shared = len(output.splits & truth.splits)
    output_only = len(output.splits - truth.splits)
    true_only = len(truth.splits - output.splits)
    raw_rf = output_only + true_only
    max_rf = max(0, 2 * (len(truth.taxa) - 3))
    rate = raw_rf / max_rf if max_rf else 0.0
    return RFResult(
        taxa=len(truth.taxa),
        output_splits=len(output.splits),
        true_splits=len(truth.splits),
        shared_splits=shared,
        output_only_splits=output_only,
        true_only_splits=true_only,
        raw_rf=raw_rf,
        max_rf=max_rf,
        rf_rate=rate,
    )


def compare_newicks(output_newick: str, true_newick: str) -> RFResult:
    truth = profile_newick(true_newick, tree_name="true species tree")
    output = profile_newick(
        output_newick,
        expected_taxa=truth.taxa,
        tree_name="output species tree",
    )
    return compare_profiles(output, truth)


def read_tree(path: Path, tree_index: int = 1) -> str:
    if tree_index < 1:
        raise ValueError("tree index must be at least 1")
    try:
        with path.open(encoding="utf-8") as handle:
            current = 0
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                current += 1
                if current == tree_index:
                    return line
    except FileNotFoundError as error:
        raise ValueError(f"tree file does not exist: {path}") from error
    raise ValueError(f"tree {tree_index} was not found in {path}")


def _print_text(result: RFResult) -> None:
    print(f"RF rate: {result.rf_rate:.10f}")
    print(f"RF distance: {result.raw_rf}/{result.max_rf}")
    print(f"Taxa: {result.taxa}")
    print(
        "Splits: "
        f"output={result.output_splits}, true={result.true_splits}, "
        f"shared={result.shared_splits}"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate the normalized unrooted RF rate between an inferred "
            "output species tree and the true species tree."
        )
    )
    parser.add_argument("output_tree", type=Path, help="inferred species-tree Newick file")
    parser.add_argument("true_tree", type=Path, help="true species-tree Newick file")
    parser.add_argument("--json", action="store_true", help="emit one JSON object")
    args = parser.parse_args(argv)

    try:
        output_newick = read_tree(args.output_tree)
        true_newick = read_tree(args.true_tree)
        result = compare_newicks(output_newick, true_newick)
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(asdict(result), sort_keys=True))
    else:
        _print_text(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
