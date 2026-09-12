#!/usr/bin/env python3
"""Stream GDL copy-labelled Newick trees into species-labelled method input.

Only exact reference species labels or SimPhy species_locus_copy labels are
accepted. Copies remain separate leaves; only their species identity changes.
Original data is never modified. Metadata records SHA-256 reproducibility IDs.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

from calculate_rf_rate import profile_newick
from uniquify_leaves import _leaf_label_spans


def trees(stream, digest):
    """Bounded-memory splitter supporting multiline trees, quotes and comments."""
    parts = []
    quoted = False
    comments = depth = 0
    for line in stream:
        digest.update(line.encode("utf-8"))
        start = i = 0
        while i < len(line):
            char = line[i]
            if comments:
                if char == "[":
                    comments += 1
                elif char == "]":
                    comments -= 1
            elif quoted:
                if char == "'":
                    if i + 1 < len(line) and line[i + 1] == "'":
                        i += 1
                    else:
                        quoted = False
            elif char == "'":
                quoted = True
            elif char == "[":
                comments = 1
            elif char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
                if depth < 0:
                    raise ValueError("unbalanced Newick parentheses")
            elif char == ";":
                if depth:
                    raise ValueError("tree ends before closing its parentheses")
                parts.append(line[start:i + 1])
                yield "".join(parts).strip()
                parts.clear()
                start = i + 1
            i += 1
        parts.append(line[start:])
    if quoted or comments or depth or "".join(parts).strip():
        raise ValueError("unterminated Newick tree, quote, or comment")


def prepare(args):
    source, reference, output = map(Path, (args.input, args.reference, args.output))
    metadata_path = Path(args.metadata)
    if output.resolve() in (source.resolve(), reference.resolve()) or metadata_path.resolve() in (
        source.resolve(), reference.resolve(), output.resolve()
    ):
        raise ValueError("normalized output/metadata must not overwrite source data or each other")
    reference_bytes = reference.read_bytes()
    reference_text = reference_bytes.decode("utf-8")
    taxa = set(profile_newick(reference_text, tree_name="reference species tree").taxa)
    if any("\n" in species or "\r" in species for species in taxa):
        raise ValueError("species labels must not contain embedded line breaks")
    observed = set()
    raw_digest, normalized_digest = hashlib.sha256(), hashlib.sha256()
    count = leaves = changed = 0
    with source.open(newline="") as stream, output.open("w") as target:
        for count, tree in enumerate(trees(stream, raw_digest), 1):
            fragments = []
            end = 0
            spans = _leaf_label_spans(tree)
            if not spans:
                raise ValueError(f"tree {count}: no leaves")
            for start, stop, token in spans:
                label = token[1:-1].replace("''", "'") if token.startswith("'") else token
                if label in taxa:
                    species = label
                else:
                    match = re.fullmatch(r"([0-9]+)_[0-9]+_[0-9]+", label)
                    if not match or match[1] not in taxa:
                        raise ValueError(f"tree {count}: cannot map leaf {label!r} to a reference species")
                    species = match[1]
                changed += species != label
                observed.add(species)
                leaves += 1
                replacement = "'" + species.replace("'", "''") + "'" if token.startswith("'") else species
                fragments.extend((tree[end:start], replacement))
                end = stop
            fragments.append(tree[end:])
            # STELAR's reader expects one gene tree per line. Whitespace folding
            # changes no topology, lengths, labels or comment semantics.
            normalized = "".join(fragments).replace("\r", " ").replace("\n", " ") + "\n"
            target.write(normalized)
            normalized_digest.update(normalized.encode())
    if not count:
        raise ValueError("combined gene-tree file is empty")
    if args.expected_genes is not None and count != args.expected_genes:
        raise ValueError(f"expected {args.expected_genes} gene trees, found {count}; refusing a partial dataset")
    metadata = dict(
        policy="exact reference species label, otherwise SimPhy species_locus_copy -> species",
        original_input=str(source.resolve()), reference_tree=str(reference.resolve()),
        input_sha256=raw_digest.hexdigest(), reference_sha256=hashlib.sha256(reference_bytes).hexdigest(),
        normalized_sha256=normalized_digest.hexdigest(), gene_trees=count, gene_leaves=leaves,
        renamed_leaves=changed, observed_species=len(observed), reference_species=len(taxa),
        missing_species=sorted(taxa - observed),
    )
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Prepared {count} gene trees; {len(observed)}/{len(taxa)} reference species; {changed} copy labels normalized.")
    if taxa - observed:
        print("Warning: species absent from the gene-tree union; strict RF may be unavailable.", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--expected-genes", type=int)
    args = parser.parse_args()
    try:
        prepare(args)
    except (OSError, ValueError) as exc:
        print(f"GDL input preparation failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
