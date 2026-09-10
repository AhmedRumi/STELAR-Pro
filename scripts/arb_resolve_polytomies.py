#!/usr/bin/env python3
"""Arbitrarily refine every Newick polytomy without changing leaf labels.

The resolver is a single left-to-right scan with an explicit stack, so a tree
of ``n`` characters is processed in O(n) time and memory regardless of its
depth or of the degree of its polytomies. It emits exactly the same text as
the earlier recursive implementation (kept in ``backup-code/``): children of
a polytomy are folded pairwise from the front, DendroPy style, each inserted
node carries the branch length ``:0.0``, and every label, branch length, and
comment of the input is preserved verbatim.
"""

from __future__ import annotations

import argparse
import sys
from collections import deque
from pathlib import Path

from uniquify_leaves import restore_leaf_names, uniquify_tree_with_mapping


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


# A resolved tree is either a leaf (``str``, copied verbatim from the input)
# or an internal node ``(left, right, suffix)`` where ``suffix`` is the text
# that followed the node's closing parenthesis (label, branch length, comment).


def _fold(children: list) -> tuple:
    """Fold a child list into a binary node, DendroPy style.

    The first two children are joined and the pair is appended to the end of
    the list until only two items remain. A deque keeps this O(k).
    """
    if len(children) == 2:
        return children[0], children[1]
    queue = deque(children)
    while len(queue) > 2:
        first = queue.popleft()
        second = queue.popleft()
        queue.append((first, second, ":0.0"))
    return queue[0], queue[1]


class _Frame:
    """Parser state for one open internal node (or the virtual root)."""

    __slots__ = ("children", "part_start", "seen", "opaque", "node", "node_end")

    def __init__(self, part_start: int) -> None:
        self.children: list = []
        self.part_start = part_start   # index where the current child part starts
        self.seen = False              # any content seen yet in the current part
        self.opaque = 0                # depth of parentheses that are plain text
        self.node = None               # (left, right) of a closed child node
        self.node_end = -1             # index just after that child's ')'

    def finish_part(self, text: str, end: int) -> None:
        if self.node is not None:
            # Keep whatever followed ')' up to the separator, minus trailing
            # whitespace, exactly as the recursive version did.
            suffix = text[self.node_end:end].rstrip()
            self.children.append((self.node[0], self.node[1], suffix))
        else:
            leaf = text[self.part_start:end].strip()
            if leaf:
                self.children.append(leaf)
        self.node = None
        self.node_end = -1
        self.seen = False
        self.opaque = 0

    def start_part(self, start: int) -> None:
        self.part_start = start


def _parse_resolved(text: str):
    """Parse ``text`` (which starts with '(') into a resolved item tree."""
    n = len(text)
    root = _Frame(0)
    frames = [root]
    i = 0
    while i < n:
        char = text[i]
        if char.isspace():
            i += 1
            continue
        frame = frames[-1]
        if char == "'":
            frame.seen = True
            i += 1
            while True:
                if i >= n:
                    raise ValueError("unbalanced quoted label in Newick input")
                if text[i] == "'":
                    if i + 1 < n and text[i + 1] == "'":
                        i += 2
                        continue
                    i += 1
                    break
                i += 1
            continue
        if char == "[":
            frame.seen = True
            depth = 1
            i += 1
            while i < n and depth:
                if text[i] == "[":
                    depth += 1
                elif text[i] == "]":
                    depth -= 1
                i += 1
            if depth:
                raise ValueError("unbalanced comment in Newick input")
            continue
        if char == "(":
            if frame.seen:
                # A '(' that is not the first token of its part is plain text
                # (the recursive version treated such parts as opaque leaves).
                frame.opaque += 1
            else:
                frame.seen = True
                frames.append(_Frame(i + 1))
        elif char == ")":
            if frame.opaque:
                frame.opaque -= 1
            else:
                frame.finish_part(text, i)
                if len(frame.children) < 2:
                    raise ValueError(
                        "an internal Newick node has fewer than two children")
                frames.pop()
                parent = frames[-1]
                parent.node = _fold(frame.children)
                parent.node_end = i + 1
                parent.seen = True
                if parent is root:
                    # The outermost node is closed. Whatever follows is its
                    # suffix and is passed through untouched, exactly as the
                    # recursive version did; resolve_file() already validated
                    # the balance of the whole file.
                    break
        elif char == "," and frame.opaque == 0:
            frame.finish_part(text, i)
            frame.start_part(i + 1)
        else:
            frame.seen = True
        i += 1
    if len(frames) != 1 or root.node is None:
        raise ValueError("malformed internal node in Newick input")
    root.finish_part(text, n)
    return root.children[0]


def _serialize(item) -> str:
    """Write a resolved item tree to Newick with one pass and one join."""
    out: list[str] = []
    stack = [item]
    while stack:
        current = stack.pop()
        if type(current) is str:
            out.append(current)
        else:
            left, right, suffix = current
            stack.append(suffix)
            stack.append(")")
            stack.append(right)
            stack.append(",")
            stack.append(left)
            stack.append("(")
    return "".join(out)


def resolve_subtree(subtree: str) -> str:
    subtree = subtree.strip()
    if not subtree.startswith("("):
        return subtree
    return _serialize(_parse_resolved(subtree))


def resolve_file(input_path: Path, output_path: Path) -> int:
    text = input_path.read_text(encoding="utf-8").replace("[&R]", "")
    trees = _scan_parts(text, ";")
    if not trees:
        raise ValueError(f"input tree file is empty: {input_path}")
    resolved = []
    for tree in trees:
        # Duplicate leaves are made unique for the topology pass and the
        # original species labels are restored immediately afterwards so
        # ASTRAL-Pro still sees gene copies.
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
