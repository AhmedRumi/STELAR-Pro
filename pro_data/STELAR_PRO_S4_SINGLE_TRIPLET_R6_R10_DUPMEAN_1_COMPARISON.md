# S4 single-triplet similarity versus S3, S5, ASTRAL-Pro3 (1-lb, R6–R10)

Dataset: `simphy-200taxa-1000gt-10rep-dupmean-1-lb-2e-10-noloss`;
200 taxa and 1,000 gene trees per replicate. RF is standard unrooted RF
against the true species tree, normalized by 394; lower is better. Times are
full-process wall seconds.

| Method | Mean RF rate (R6–R10) | Mean runtime |
|:--|--:|--:|
| STELAR-Pro S3 | 3.25% | 46.95 s |
| STELAR-Pro S4, one triplet matrix | 3.25% | 53.39 s |
| STELAR-Pro S5, overall optimization | 3.45% | 145.16 s |
| ASTRAL-Pro3 | **2.23%** | 73.72 s |

| Replicate | S3 RF | New S4 RF | S5 RF | ASTRAL-Pro3 RF | S3 time | New S4 time | S5 time | ASTRAL-Pro3 time |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| R6 | 6 | 6 | 6 | 4 | 28.04 | 34.37 | 109.95 | 67.41 |
| R7 | 20 | 20 | 24 | 10 | 67.46 | 75.73 | 184.34 | 73.80 |
| R8 | 12 | 12 | 12 | 8 | 43.47 | 48.48 | 136.42 | 74.49 |
| R9 | 18 | 18 | 18 | 14 | 73.23 | 82.24 | 192.30 | 73.31 |
| R10 | 8 | 8 | 8 | 8 | 22.54 | 26.14 | 102.78 | 79.62 |

RF table entries are raw counts out of 394. Every new S4 tree has pairwise
unrooted RF 0 against its S3 counterpart (and the prior S4 output). Thus the
corrected guide changed the matrix route but did not change unrooted output
topologies or accuracy in these five replicates. The prior dual-matrix S4
runs averaged 49.83 s versus 53.39 s here; these separate-run timings do
not establish a speed regression.

S4 now builds **one** rooted-triplet sister-frequency matrix from the
single-copy DISCO trees and reuses it for incomplete-tree completion, UPGMA,
and greedy consensus. Its distance view is `1 − similarity`. S2/S3 still
use their quartet similarity matrix. The completion algorithm still uses
its four-point insertion rule, but reads S4's triplet matrix; this change
does not replace that insertion rule.

Sanity checks: Java build, CPU/CUDA backend check, CLI tests (including an
S4 single-matrix route assertion), and hand-calculated CPU/GPU/packed
triplet-matrix parity all passed. A strict-GPU incomplete-tree S4 run passed.
Each R6–R10 log records one GPU rooted-triplet similarity build; all five
benchmark runs completed below five minutes. A separate broad regression
script stopped before reaching S4 on an unrelated S1 score-only assertion
(observed 27, expected 21); it is not counted as a passing suite.

New S4 trees, times, and logs: `benchmark-stelar-pro-s4-single-triplet-gpu-dupmean-1/R6`
through `R10`. The S3 and updated S5 results are from their existing
1-lb GPU benchmarks; ASTRAL-Pro3 results are from
`benchmark-stelar-pro-vs-astral-pro3-gpu/replicate-results.csv`.
STELAR-Pro used strict CUDA and 28 CPU threads; the reused ASTRAL-Pro3
run used 28 CPU threads and seed 42. Runtimes are indicative because they
were recorded in separate runs with different compute backends.
