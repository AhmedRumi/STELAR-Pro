#!/usr/bin/env python3
"""Fixture for the GDL runner; never used by production code."""
import os
from pathlib import Path
import sys

args = sys.argv[1:]
source = Path(args[args.index("-i") + 1])
target = Path(args[args.index("-o") + 1])
with open(os.environ["GDL_CALLS"], "a") as stream:
    stream.write(" ".join(args) + "\n")
text = source.read_text()
assert "_0_" not in text, "runner failed to normalize copies"
assert text.count(";") == 2
assert "(0,0)" in text, "runner collapsed distinct copies"
if os.environ.get("GDL_FAIL") == "1":
    sys.exit(9)
if os.environ.get("GDL_EMPTY") != "1":
    target.write_text("((0,1),(2,3));\n")
