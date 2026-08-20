# Validation strategy

The migration is checked at three levels.

## Exact small-tree oracle

`test/test_stelarx_triplets.py` independently parses rooted Newick, enumerates
all taxon triples, determines the rooted pair by LCA, and compares that total to
STELAR-X fixed-tree scoring. It does not reuse production weight code.

The five-taxon fixture in `test/input/test_5taxa.tre` and candidate
`test/input/stelar_candidate_5taxa.tre` have an independently enumerated score
of 21.

## Scoring-path parity

The fixture is scored with I1, I2, I3, and I4. All CPU and CUDA implementations
must return the same exact score. CUDA uses the same doubled formulas for binary
and rooted-polytomous nodes.

`test/input/stelar_polytomy_5taxa.tre` separately exercises the generalized
rooted-polytomy formula against the same enumeration oracle in every
intersection mode. Unresolved triples are excluded rather than treated as
agreements.

`test/input/stelar_polytomy_incomplete_6taxa.tre` adds multiple internal
polytomies and missing taxa. Its independently enumerated score is 22, checking
that taxa outside a rooted node—and taxa absent from a gene tree—do not become
synthetic child groups.

## Rootedness and root preservation

- A conventional unrooted Newick root with three children must be rejected.
- Rooted input with two top-level children must parse without rerooting.
- Completion tests verify that restricting a completed tree to its original
  taxa recovers the original rooted clades.

Run the full focused suite with `test/run_stelarx_tests.sh`. It builds the Java
code, runs the independent binary and polytomy oracle, checks that completion
preserves every original rooted triple, exercises the cluster/partition/DP
verifiers, runs S1/S2/S3, checks product naming, and probes the native JNI
libraries. The narrower oracle can be run alone with
`python3 test/test_stelarx_triplets.py`.

## Representative resource comparison

The migrated working tree and the original `HEAD` revision were compiled
separately and run with the same JDK, one CPU thread, 64 MiB initial heap, 2 GiB
heap limit, S1/I2, and the 200-tree `all_gt_bs_rooted_37.tre` fixture.
`/usr/bin/time` reported:

| build | wall time | maximum RSS |
|---|---:|---:|
| original ASTRAL-X | 0.78 s | 164,732 KiB |
| migrated STELAR-X | 0.73 s | 129,240 KiB |

This representative run is about 6% faster and uses about 22% less peak RSS.
The structural reason is also checked directly: STELAR-X retains descendant
clades only and emits rooted child transitions only, eliminating the
complement orientations and Type-2 rotations required by the unrooted search.
The asymptotic bounds of the retained hash, DP, similarity, UPGMA, consensus,
and cross-tree implementations are unchanged.

## CUDA hardware validation

`test/run_stelarx_gpu_tests.sh` runs with `--gpu-strict`, which prevents an
unavailable CUDA device from being silently replaced by the normal CPU mode.
It was run successfully on an NVIDIA GeForce RTX 3050 6GB Laptop GPU (compute
capability 8.6), driver 580.173.02, driver API 13.0, and bundled CUDA runtime
12.0.

- I1, I2, I3, and I4 each executed their CUDA weight kernel and returned the
  independent binary score 21, rooted-polytomy score 11, and incomplete
  rooted-polytomy score 22.
- On 200 rooted gene trees with 37 taxa, every CUDA method returned 1,390,544,
  exactly matching the CPU result for `true_37.tre`.
- Strict S1, S2, and S3 inference all returned 122 on the incomplete-tree
  fixture. S2/S3 additionally exercised the GPU similarity matrix and GPU
  cross-tree DP construction.
- End-to-end S2 inference preserving native input polytomies returned 27 on CPU
  and strict CUDA and produced byte-identical species-tree topology.

All rooted-polytomy CUDA paths use the linear-time formula over actual children:
the prefix-sum, smaller-side, bitset, and tree-walk kernels, including their
long/double and Int128 variants.

## Comprehensive suite

The broader regression entry point is `test/run_stelarx_comprehensive_tests.sh`.
It adds seeded independent-oracle scoring and inference, forced numeric modes,
matrix/UPGMA oracles, parser and path-safety failures, low-level arithmetic and
large-matrix boundaries, portable packaging, CUDA batching controls, and strict
randomized hardware parity. See `03-comprehensive-testing.md` for the coverage
map and fast/full invocations.
