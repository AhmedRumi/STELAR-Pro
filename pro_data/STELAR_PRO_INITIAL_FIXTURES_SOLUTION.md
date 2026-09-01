# STELAR-Pro initial multicopy-fixture solution

This report independently checks the current STELAR-Pro constrained S1 pipeline
on the two no-loss multicopy inputs added on 2026-09-01.

## Calculation rules

For every rooted/tagged gene tree:

1. Visit every internal node after rooting/tagging.
2. Ignore a node tagged `D`.
3. At an untagged speciation node, replace each child's leaf-copy multiset by
   its distinct species set. For example,
   `[1,2,1,2] | [3,4] -> {1,2} | {3,4}`.
4. Treat the two child sets as an unordered candidate bipartition. Repeated
   occurrences are retained in the frequency but share one candidate identity.

For candidate `A | B` and gene bipartition `M0 | M1`, the independent weight
calculation used the STELAR rooted-triplet formula

```text
2 w(A|B, M0|M1)
  = |A∩M0||B∩M1|(|A∩M0| + |B∩M1| - 2)
  + |A∩M1||B∩M0|(|A∩M1| + |B∩M0| - 2).
```

Weights below sum this quantity over every selected speciation-node occurrence
and divide the final sum by two. The independent oracle parses Newick itself and
uses Python sets; it does not use STELAR-Pro hashing, `PartitionTable`, the Java
tree parser, or the Java weight/DP implementation.

## Overall results

| Measurement | 7 taxa / 5 trees | 12 taxa / 2 trees |
|---|---:|---:|
| Leaf copies per gene tree | 11 | 24 |
| Duplication nodes skipped | 10 total (2/tree) | 12 total (6/tree) |
| Selected speciation occurrences | 40 | 34 |
| Unique candidate bipartitions | 27 | 15 |
| Parent clusters with splits | 23 | 14 |
| Clusters in `X` | 29 | 25 |
| Maximum candidate weight | 62 | 128 |
| Sum of all candidate weights | 365 | 616 |
| Independently calculated optimum | 91 | 312 |
| STELAR-Pro score | 91 | 312 |

Both input files are already binary. The arbitrary resolver returned each file
byte-for-byte unchanged before rooting/tagging.

## Dataset 1: seven taxa, five gene trees

Input: `test_gene_trees_7taxa_dup_no_loss.tre`

### Rooted/tagged trees

Branch lengths are omitted here because they do not enter the calculation.

```text
T1: (((5,6),(5,6))D,(7,(((1,2),(1,2))D,(3,4))));
T2: (2,(((1,3),(1,3))D,(((4,5),(4,5))D,(6,7))));
T3: (6,(((2,3),(2,3))D,((1,((4,5),(4,5))D),7)));
T4: ((((1,7),(1,7))D,(2,6)),(3,((4,5),(4,5))D));
T5: (2,(1,(((3,7),(3,7))D,(((4,6),(4,6))D,5))));
```

### Per-gene-tree candidate calculation

Node numbers are postorder internal-node numbers. `x2` means two distinct
speciation roots produce the same duplicate-collapsed child bipartition.

#### Tree 1

Selected, 8 occurrences:

- N01,N02: `{5} | {6}` (`x2`)
- N04,N05: `{1} | {2}` (`x2`)
- N07: `{3} | {4}`
- N08: `[1,2,1,2] | [3,4] -> {1,2} | {3,4}`
- N09: `[7] | [1,2,1,2,3,4] -> {7} | {1,2,3,4}`
- N10, root: `[5,6,5,6] | [7,1,2,1,2,3,4] -> {5,6} | {1,2,3,4,7}`

Skipped duplication roots:

- N03: `{5,6} | {5,6}`
- N06: `{1,2} | {1,2}`

#### Tree 2

Selected, 8 occurrences:

- N01,N02: `{1} | {3}` (`x2`)
- N04,N05: `{4} | {5}` (`x2`)
- N07: `{6} | {7}`
- N08: `[4,5,4,5] | [6,7] -> {4,5} | {6,7}`
- N09: `[1,3,1,3] | [4,5,4,5,6,7] -> {1,3} | {4,5,6,7}`
- N10, root: `{2} | {1,3,4,5,6,7}`

Skipped duplication roots:

- N03: `{1,3} | {1,3}`
- N06: `{4,5} | {4,5}`

#### Tree 3

Selected, 8 occurrences:

- N01,N02: `{2} | {3}` (`x2`)
- N04,N05: `{4} | {5}` (`x2`)
- N07: `[1] | [4,5,4,5] -> {1} | {4,5}`
- N08: `[1,4,5,4,5] | [7] -> {1,4,5} | {7}`
- N09: `[2,3,2,3] | [1,4,5,4,5,7] -> {2,3} | {1,4,5,7}`
- N10, root: `{6} | {1,2,3,4,5,7}`

Skipped duplication roots:

- N03: `{2,3} | {2,3}`
- N06: `{4,5} | {4,5}`

#### Tree 4

Selected, 8 occurrences:

- N01,N02: `{1} | {7}` (`x2`)
- N04: `{2} | {6}`
- N05: `[1,7,1,7] | [2,6] -> {1,7} | {2,6}`
- N06,N07: `{4} | {5}` (`x2`)
- N09: `[3] | [4,5,4,5] -> {3} | {4,5}`
- N10, root: `{1,2,6,7} | {3,4,5}`

Skipped duplication roots:

- N03: `{1,7} | {1,7}`
- N08: `{4,5} | {4,5}`

#### Tree 5

Selected, 8 occurrences:

- N01,N02: `{3} | {7}` (`x2`)
- N04,N05: `{4} | {6}` (`x2`)
- N07: `[4,6,4,6] | [5] -> {4,6} | {5}`
- N08: `[3,7,3,7] | [4,6,4,6,5] -> {3,7} | {4,5,6}`
- N09: `{1} | {3,4,5,6,7}`
- N10, root: `{2} | {1,3,4,5,6,7}`

Skipped duplication roots:

- N03: `{3,7} | {3,7}`
- N06: `{4,6} | {4,6}`

### All 27 unique candidates and independent weights

`f` is the selected-occurrence frequency. `T1..T5` are the independently
calculated per-gene-tree weight contributions.

| ID | Unordered child bipartition | f | T1 | T2 | T3 | T4 | T5 | Weight |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| C01 | `{1,2,3,4,5,7} | {6}` | 1 | 10 | 2 | 15 | 4 | 1 | 32 |
| C02 | `{1,2,3,4,7} | {5,6}` | 1 | 25 | 5 | 11 | 5 | 6 | 52 |
| C03 | `{1,2,3,4} | {7}` | 1 | 6 | 1 | 2 | 1 | 0 | 10 |
| C04 | `{1,2,6,7} | {3,4,5}` | 1 | 9 | 9 | 7 | 30 | 7 | 62 |
| C05 | `{1,2} | {3,4}` | 1 | 4 | 1 | 0 | 4 | 2 | 11 |
| C06 | `{1,3,4,5,6,7} | {2}` | 2 | 2 | 15 | 6 | 4 | 15 | 42 |
| C07 | `{1,3} | {4,5,6,7}` | 1 | 5 | 16 | 5 | 3 | 9 | 38 |
| C08 | `{1,4,5,7} | {2,3}` | 1 | 2 | 9 | 16 | 4 | 7 | 38 |
| C09 | `{1,4,5} | {7}` | 1 | 1 | 1 | 3 | 1 | 1 | 7 |
| C10 | `{1,7} | {2,6}` | 1 | 1 | 1 | 2 | 4 | 1 | 9 |
| C11 | `{1} | {2}` | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| C12 | `{1} | {3,4,5,6,7}` | 1 | 2 | 6 | 1 | 3 | 10 | 22 |
| C13 | `{1} | {3}` | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| C14 | `{1} | {4,5}` | 1 | 0 | 1 | 1 | 1 | 1 | 4 |
| C15 | `{1} | {7}` | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| C16 | `{2} | {3}` | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| C17 | `{2} | {6}` | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| C18 | `{3,7} | {4,5,6}` | 1 | 4 | 4 | 3 | 2 | 9 | 22 |
| C19 | `{3} | {4,5}` | 1 | 0 | 1 | 1 | 1 | 1 | 4 |
| C20 | `{3} | {4}` | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| C21 | `{3} | {7}` | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| C22 | `{4,5} | {6,7}` | 1 | 0 | 4 | 2 | 4 | 1 | 11 |
| C23 | `{4,6} | {5}` | 1 | 0 | 0 | 0 | 0 | 1 | 1 |
| C24 | `{4} | {5}` | 6 | 0 | 0 | 0 | 0 | 0 | 0 |
| C25 | `{4} | {6}` | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| C26 | `{5} | {6}` | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| C27 | `{6} | {7}` | 1 | 0 | 0 | 0 | 0 | 0 | 0 |

### Inference check

STELAR-Pro returned:

```text
(2,((1,3),((5,4),(6,7))));
```

Its internal candidate weights are
`42 + 38 + 0 + 11 + 0 + 0 = 91`. The independent DP obtains the same
topology modulo child order and the same optimum `91`.

## Dataset 2: twelve taxa, two gene trees

Input: `test_gene_trees_12taxa_2trees_dup_no_loss.tre`

### Rooted/tagged trees

```text
T1: (((9,10),(9,10))D,(((11,12),(11,12))D,((((1,2),(1,2))D,((3,4),(3,4))D),(((5,6),(5,6))D,((7,8),(7,8))D))));
T2: (((7,8),(7,8))D,(((5,6),(5,6))D,((((1,2),(1,2))D,((3,4),(3,4))D),(((9,10),(9,10))D,((11,12),(11,12))D))));
```

### Per-gene-tree candidate calculation

#### Tree 1

Selected, 17 occurrences:

- N01,N02: `{9} | {10}` (`x2`)
- N04,N05: `{11} | {12}` (`x2`)
- N07,N08: `{1} | {2}` (`x2`)
- N10,N11: `{3} | {4}` (`x2`)
- N13: `{1,2} | {3,4}`
- N14,N15: `{5} | {6}` (`x2`)
- N17,N18: `{7} | {8}` (`x2`)
- N20: `{5,6} | {7,8}`
- N21: `{1,2,3,4} | {5,6,7,8}`
- N22: `{11,12} | {1,2,3,4,5,6,7,8}`
- N23, root: `{9,10} | {1,2,3,4,5,6,7,8,11,12}`

Skipped duplication roots: `{9,10}|{9,10}`, `{11,12}|{11,12}`,
`{1,2}|{1,2}`, `{3,4}|{3,4}`, `{5,6}|{5,6}`, and `{7,8}|{7,8}`.

#### Tree 2

Selected, 17 occurrences:

- N01,N02: `{7} | {8}` (`x2`)
- N04,N05: `{5} | {6}` (`x2`)
- N07,N08: `{1} | {2}` (`x2`)
- N10,N11: `{3} | {4}` (`x2`)
- N13: `{1,2} | {3,4}`
- N14,N15: `{9} | {10}` (`x2`)
- N17,N18: `{11} | {12}` (`x2`)
- N20: `{9,10} | {11,12}`
- N21: `{1,2,3,4} | {9,10,11,12}`
- N22: `{5,6} | {1,2,3,4,9,10,11,12}`
- N23, root: `{7,8} | {1,2,3,4,5,6,9,10,11,12}`

Skipped duplication roots: `{7,8}|{7,8}`, `{5,6}|{5,6}`,
`{1,2}|{1,2}`, `{3,4}|{3,4}`, `{9,10}|{9,10}`, and
`{11,12}|{11,12}`.

### All 15 unique candidates and independent weights

| ID | Unordered child bipartition | f | Tree 1 | Tree 2 | Weight |
|---|---|---:|---:|---:|---:|
| C01 | `{1,2,3,4,5,6,7,8,11,12} | {9,10}` | 1 | 100 | 28 | 128 |
| C02 | `{1,2,3,4,5,6,7,8} | {11,12}` | 1 | 64 | 24 | 88 |
| C03 | `{1,2,3,4,5,6,9,10,11,12} | {7,8}` | 1 | 28 | 100 | 128 |
| C04 | `{1,2,3,4,9,10,11,12} | {5,6}` | 1 | 24 | 64 | 88 |
| C05 | `{1,2,3,4} | {5,6,7,8}` | 1 | 48 | 32 | 80 |
| C06 | `{1,2,3,4} | {9,10,11,12}` | 1 | 32 | 48 | 80 |
| C07 | `{1,2} | {3,4}` | 2 | 4 | 4 | 8 |
| C08 | `{9} | {10}` | 4 | 0 | 0 | 0 |
| C09 | `{9,10} | {11,12}` | 1 | 4 | 4 | 8 |
| C10 | `{11} | {12}` | 4 | 0 | 0 | 0 |
| C11 | `{1} | {2}` | 4 | 0 | 0 | 0 |
| C12 | `{3} | {4}` | 4 | 0 | 0 | 0 |
| C13 | `{5,6} | {7,8}` | 1 | 4 | 4 | 8 |
| C14 | `{5} | {6}` | 4 | 0 | 0 | 0 |
| C15 | `{7} | {8}` | 4 | 0 | 0 | 0 |

### Inference check

STELAR-Pro returned:

```text
((9,10),((12,11),(((5,6),(8,7)),((4,3),(1,2)))));
```

Its internal candidate weights are
`128 + 0 + 88 + 0 + 80 + 8 + 0 + 0 + 8 + 0 + 0 = 312`.
The independent DP also obtains optimum `312`.

There are tied optima because the two input trees place the `{5,6}/{7,8}` and
`{9,10}/{11,12}` blocks symmetrically. The independent oracle's deterministic
tie choice may therefore have a different topology, but the STELAR-Pro output
contains only listed candidates and independently evaluates to the optimum.

## Exact implementation/oracle comparison

For each dataset, a read-only Java inspector extracted the production
`PartitionTable` as canonical `child-set-1 | child-set-2, frequency` rows. A
separate Python Newick/set oracle generated the same rows without using project
parsing, hashing, partition, weight, or DP code.

```text
7 taxa:  implementation rows = 27, oracle rows = 27, exact diff = empty
12 taxa: implementation rows = 15, oracle rows = 15, exact diff = empty
```

The built-in Phase 4 and Phase 5 verifiers also reported `ALL ASSERTIONS PASSED`
for both datasets. Therefore the selected candidate identities, frequencies,
duplicate collapse, duplication filtering, candidate weights, reachable S1 DP
space, inferred objective, and output-tree objective all match the independent
calculations on these fixtures.
