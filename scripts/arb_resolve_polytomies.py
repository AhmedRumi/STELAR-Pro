#!/usr/bin/env python3
"""Arbitrarily refine every Newick polytomy without changing leaf labels."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from uniquify_leaves import restore_leaf_names, uniquify_tree_with_mapping

sys.setrecursionlimit(10_000_000)


def _scan_parts(text: str, separator: str) -> list[str]:
    """Split outside parentheses, comments, and quoted labels."""
    parts: list[str] = []
    start = 0
    depth = 0
    comment_depth = 0
    quoted = False
    i = 0
    while i < len(text):
        char = text[i]
        if quoted:
            if char == "'":
                if i + 1 < len(text) and text[i + 1] == "'":
                    i += 1
                else:
                    quoted = False
        elif comment_depth:
            if char == "[":
                comment_depth += 1
            elif char == "]":
                comment_depth -= 1
        elif char == "'":
            quoted = True
        elif char == "[":
            comment_depth = 1
        elif char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                raise ValueError("unbalanced ')' in Newick input")
        elif char == separator and depth == 0:
            part = text[start:i].strip()
            if part:
                parts.append(part)
            start = i + 1
        i += 1
    if quoted or comment_depth or depth:
        raise ValueError("unbalanced Newick input")
    tail = text[start:].strip()
    if tail:
        parts.append(tail)
    return parts


def _outer_close(subtree: str) -> int:
    depth = 0
    comment_depth = 0
    quoted = False
    i = 0
    while i < len(subtree):
        char = subtree[i]
        if quoted:
            if char == "'":
                if i + 1 < len(subtree) and subtree[i + 1] == "'":
                    i += 1
                else:
                    quoted = False
        elif comment_depth:
            if char == "[":
                comment_depth += 1
            elif char == "]":
                comment_depth -= 1
        elif char == "'":
            quoted = True
        elif char == "[":
            comment_depth = 1
        elif char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0:
                return i
            if depth < 0:
                break
        i += 1
    raise ValueError("malformed internal node in Newick input")


def resolve_subtree(subtree: str) -> str:
    subtree = subtree.strip()
    if not subtree.startswith("("):
        return subtree

    close = _outer_close(subtree)
    children = _scan_parts(subtree[1:close], ",")
    if len(children) < 2:
        raise ValueError("an internal Newick node has fewer than two children")
    children = [resolve_subtree(child) for child in children]

    # Match DendroPy's deterministic rule: join the first pair and append it.
    while len(children) > 2:
        first = children.pop(0)
        second = children.pop(0)
        children.append(f"({first},{second}):0.0")
    return f"({children[0]},{children[1]}){subtree[close + 1:]}"


def resolve_file(input_path: Path, output_path: Path) -> int:
    text = input_path.read_text(encoding="utf-8").replace("[&R]", "")
    trees = _scan_parts(text, ";")
    if not trees:
        raise ValueError(f"input tree file is empty: {input_path}")
    resolved = []
    for tree in trees:
        # DendroPy requires unique leaves. Restore the original species labels
        # immediately after resolution so ASTRAL-Pro still sees gene copies.
        unique_tree, restoration = uniquify_tree_with_mapping(tree)
        resolved_tree = resolve_subtree(unique_tree)
        resolved.append(restore_leaf_names(resolved_tree, restoration) + ";\n")
    output_path.write_text("".join(resolved), encoding="utf-8")
    return len(trees)


def main() -> int:
    parser = argparse.ArgumentParser(description="Arbitrarily resolve Newick polytomies")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path, nargs="?")
    args = parser.parse_args()
    output = args.output or Path(str(args.input) + ".resolved")
    try:
        resolve_file(args.input, output)
    except (OSError, ValueError) as error:
        parser.exit(2, f"polytomy resolution failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
