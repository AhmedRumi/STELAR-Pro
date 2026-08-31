#!/usr/bin/env python3
"""Uniquify duplicate Newick leaf labels by adding numeric suffixes.

This script processes Newick format gene trees where multi-copy genes may have
the same name (e.g., ((speciesA,speciesB),(speciesB,speciesC))).

Logic:
1. Validate leaf labels containing underscores. A label with one underscore is
   treated as species_copy. Labels with more than one underscore are ambiguous
   and rejected for now.
2. Check for duplicate plain species labels in each tree.
3. Add _1, _2, _3, etc. suffixes to duplicate plain species labels.

Usage:
    python3 uniquify_leaves.py <input_file> <output_file>

The resolver imports the reversible API below and restores these temporary
labels before ASTRAL-Pro rooting/tagging.
"""

import sys
from collections import Counter

DELIMITERS = set('(),;:[] \t\n\r')


def _leaf_label_spans(newick_str):
    """Return ``(start, end, label)`` for each leaf label token."""
    spans = []
    expect_leaf = True
    i = 0
    while i < len(newick_str):
        char = newick_str[i]
        if char in "(,":
            expect_leaf = True
            i += 1
            continue
        if char == ")":
            expect_leaf = False
            i += 1
            continue
        if char == "[":
            depth = 1
            i += 1
            while i < len(newick_str) and depth:
                if newick_str[i] == "[":
                    depth += 1
                elif newick_str[i] == "]":
                    depth -= 1
                i += 1
            continue
        if char.isspace() or char in ";:":
            i += 1
            continue
        if not expect_leaf:
            i += 1
            continue

        start = i
        if char == "'":
            i += 1
            while i < len(newick_str):
                if newick_str[i] == "'":
                    if i + 1 < len(newick_str) and newick_str[i + 1] == "'":
                        i += 2
                        continue
                    i += 1
                    break
                i += 1
        else:
            while i < len(newick_str) and newick_str[i] not in DELIMITERS:
                i += 1
        spans.append((start, i, newick_str[start:i]))
        expect_leaf = False
    return spans


def _add_suffix(label, suffix):
    """Add a suffix inside a quoted label, or after an unquoted label."""
    if len(label) >= 2 and label[0] == "'" and label[-1] == "'":
        return label[:-1] + suffix + "'"
    return label + suffix


def previous_non_whitespace(text, index):
    """Return the previous non-whitespace character before index, or None."""
    index -= 1
    while index >= 0 and text[index] in ' \t\n\r':
        index -= 1
    return text[index] if index >= 0 else None

def extract_leaf_names(newick_str):
    """Extract all leaf names from a Newick string."""
    return [label for _, _, label in _leaf_label_spans(newick_str)]

def underscored_leaf_names(trees):
    """Return leaf names containing underscores across all trees."""
    underscored = []
    for tree_index, tree in enumerate(trees, 1):
        leaves = extract_leaf_names(tree)
        for leaf in leaves:
            if '_' in leaf:
                underscored.append((tree_index, leaf))
    return underscored

def ambiguous_copy_labels(trees):
    """Return underscore labels that are not valid species_copy labels."""
    ambiguous = []
    for tree_index, tree in enumerate(trees, 1):
        leaves = extract_leaf_names(tree)
        for leaf in leaves:
            if '_' not in leaf:
                continue

            parts = leaf.split('_')
            if len(parts) != 2 or not parts[0] or not parts[1]:
                ambiguous.append((tree_index, leaf))

    return ambiguous

def duplicate_copy_labels(trees):
    """Return exact duplicate labels that already use species_copy form."""
    duplicates = []
    for tree_index, tree in enumerate(trees, 1):
        leaves = extract_leaf_names(tree)
        counts = Counter(leaves)
        for leaf, count in counts.items():
            if '_' in leaf and count > 1:
                duplicates.append((tree_index, leaf, count))

    return duplicates

def uniquify_tree_with_mapping(newick_str):
    """Return a uniquified tree and the generated-label restoration map."""
    spans = _leaf_label_spans(newick_str)
    leaves = [label for _, _, label in spans]
    leaf_counts = Counter(leaves)
    duplicates = {leaf for leaf, count in leaf_counts.items() if count > 1}
    if not duplicates:
        return newick_str, {}

    result = []
    replacement_counter = {leaf: 0 for leaf in duplicates}
    used_labels = set(leaves)
    restoration = {}
    cursor = 0
    for start, end, label in spans:
        result.append(newick_str[cursor:start])
        replacement = label
        if label in duplicates:
            while True:
                replacement_counter[label] += 1
                replacement = _add_suffix(
                    label, f"_{replacement_counter[label]}")
                if replacement not in used_labels:
                    used_labels.add(replacement)
                    restoration[replacement] = label
                    break
        result.append(replacement)
        cursor = end
    result.append(newick_str[cursor:])
    return "".join(result), restoration


def uniquify_tree(newick_str):
    """Add suffixes to duplicate leaf names in a single tree."""
    return uniquify_tree_with_mapping(newick_str)[0]


def restore_leaf_names(newick_str, restoration):
    """Restore generated leaf labels after topology-only processing."""
    if not restoration:
        return newick_str
    result = []
    cursor = 0
    for start, end, label in _leaf_label_spans(newick_str):
        result.append(newick_str[cursor:start])
        result.append(restoration.get(label, label))
        cursor = end
    result.append(newick_str[cursor:])
    return "".join(result)

def process_trees(input_file, output_file):
    """Process all trees in the input file."""
    with open(input_file, 'r') as f:
        trees = [line.strip() for line in f if line.strip()]
    
    # Labels with one underscore are accepted as species_copy. More underscores
    # are ambiguous because species names themselves cannot contain underscores.
    ambiguous = ambiguous_copy_labels(trees)
    if ambiguous:
        print("Error: ambiguous underscores in input leaf labels.", file=sys.stderr)
        print('wQFM-GDL treats "_" as the species/copy separator and supports labels as species_copy.', file=sys.stderr)
        print('Species labels and copy IDs must not themselves contain "_".', file=sys.stderr)
        print("Examples:", file=sys.stderr)
        for tree_index, leaf in ambiguous[:5]:
            print(f"  tree {tree_index}: {leaf}", file=sys.stderr)
        if len(ambiguous) > 5:
            print(f"  ... and {len(ambiguous) - 5} more labels", file=sys.stderr)
        sys.exit(1)

    duplicate_copies = duplicate_copy_labels(trees)
    if duplicate_copies:
        print("Error: duplicate copy-labeled leaf labels found.", file=sys.stderr)
        print("Input labels in species_copy form must be unique within each gene tree.", file=sys.stderr)
        print("Examples:", file=sys.stderr)
        for tree_index, leaf, count in duplicate_copies[:5]:
            print(f"  tree {tree_index}: {leaf} appears {count} times", file=sys.stderr)
        if len(duplicate_copies) > 5:
            print(f"  ... and {len(duplicate_copies) - 5} more labels", file=sys.stderr)
        sys.exit(1)
    
    print("Checking for duplicate plain species labels...")
    
    # Process each tree to uniquify duplicates
    modified_trees = []
    trees_with_duplicates = 0
    duplicate_species = 0
    
    for i, tree in enumerate(trees, 1):
        leaves = extract_leaf_names(tree)
        leaf_counts = Counter(leaves)
        duplicates_in_tree = sum(1 for count in leaf_counts.values() if count > 1)
        
        if duplicates_in_tree > 0:
            trees_with_duplicates += 1
            duplicate_species += duplicates_in_tree
            modified_tree = uniquify_tree(tree)
            modified_trees.append(modified_tree)
            if i <= 5:  # Show first few examples
                print(f"  Tree {i}: Found {duplicates_in_tree} species with duplicates")
        else:
            modified_trees.append(tree)
    
    print(f"Trees with duplicates: {trees_with_duplicates}")
    print(f"Duplicate-label groups: {duplicate_species}")
    print(f"Writing uniquified trees to {output_file}")
    
    # Write output
    with open(output_file, 'w') as f:
        for tree in modified_trees:
            f.write(tree + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 uniquify_leaves.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_trees(input_file, output_file)
    print("Done!")
