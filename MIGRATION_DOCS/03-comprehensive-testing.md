# Comprehensive validation

The primary validation entry point is:

```bash
test/run_stelarx_comprehensive_tests.sh --require-gpu
```

It builds once and runs these independent layers:

1. Exact `Int128` arithmetic against `BigInteger` (100,016 checks), numeric-mode
   boundary/override checks, threading failure propagation, packed-matrix index
   boundaries, heap preflight, and wide Euler/RMQ boundaries.
2. Rooted internal-polytomy lifecycle checks: strict binary root retention,
   native versus deterministic refinement, serial/parallel determinism,
   dense/packed completion parity, and preservation of every originally
   displayed or unresolved rooted triple.
3. Seeded differential scoring against direct rooted-triple enumeration. The
   cases mix complete/incomplete trees, binary/deep/caterpillar shapes, internal
   degrees 3–5, duplicates, child/line permutations, and 4–10 taxa. Every case
   runs through I1, I2, I3, and I4. Representative polytomous cases additionally
   force the normally large-N-only Double and Int128 paths.
4. End-to-end S1/S2/S3 inference through every intersection method. The final
   output taxon set/root and reported objective are independently checked, as
   are thread-count, hash-seed-count, pruning, polytomy-refinement, and native
   polytomy invariants.
5. Independent randomized similarity, topological-distance, and UPGMA oracles;
   packed/dense matrix and RMQ parity; focused cluster/partition/DP verifiers.
6. CLI and parser rejection cases: empty/malformed files, duplicate taxa,
   non-binary supplied roots, incomplete/unknown/multiple species trees,
   invalid methods/presets/numeric modes, `--unrooted`, and all material
   read/write path collisions. Decorated Newick branch lengths and internal
   labels are scored through all four methods.
7. A self-contained CPU portable image build plus its packaged launcher,
   inference, log-file, and version smoke tests.
8. Strict CUDA tests (no fallback): I1–I4 on fixed and randomized binary,
   incomplete, and polytomous cases; forced Long/Double/Int128 kernels; automatic,
   fixed-size, fixed-count, and disabled batching; 37-taxon parity; S1/S2/S3;
   CPU/GPU topology parity; and randomized GPU similarity/distance matrices.

For a shorter edit-time check, use:

```bash
test/run_stelarx_comprehensive_tests.sh --quick --cpu-only
```

Random generation is deterministic (`0x5e1a7` by default), and a failing case
prints its case name, method, numeric mode, expected score, and captured program
output. The number of generated CPU cases can be changed with the differential
test's `--cases`; the strict GPU wrapper accepts `STELARX_GPU_RANDOM_CASES`.

On the RTX 3050 development host, the complete `--require-gpu` run took 259
seconds (4 minutes 19 seconds), including a clean Java build and portable-image
construction. The `--quick --cpu-only` tier took 71 seconds.

Three production defects were exposed while introducing this suite and are now
permanently covered: unsigned-low-limb multiplication in `Int128.mulScalar`,
acceptance of empty inputs, and unsafe read/write path aliases.

## Accuracy and scalability regression benchmark

`test/test_stelarx_scalability.py` adds measured performance checks without
placing inherently noisy wall-time assertions in the ordinary correctness
suite. It:

- independently enumerates all 1,554,000 gene-tree/triple observations in the
  37-taxon, 200-gene fixture;
- runs S1/S2/S3 × I1/I2/I3/I4 on CPU and strict CUDA, requiring byte-identical
  topology and exact score parity across implementations;
- independently re-scores each CPU preset's inferred tree;
- replicates the gene collection by 1×, 4×, and 16× (200, 800, and 3,200 trees),
  requiring the score to scale exactly linearly for every CPU/CUDA intersection
  method; and
- records external wall time and peak resident memory with `/usr/bin/time`.

The optional `--reference-dir` comparison uses three-run medians and fails on a
time ratio over 1.35× or an RSS ratio over 1.25×. Against clean pre-migration
commit `88af054` on the RTX 3050 host, S1/I2 CPU measured:

| Build | Median wall time | Median peak RSS |
|---|---:|---:|
| Pre-migration ASTRAL-X reference | 1.110 s | 360.5 MiB |
| Current STELAR-X | 1.100 s | 204.5 MiB |
| Current/reference | 0.991× | 0.567× |

All 24 S×I×CPU/CUDA inference runs returned score 1,391,309 with identical
topology. The independent fixed-tree oracle returned 1,390,544. Across fixed
scoring, 1×/4×/16× returned exactly 1,390,544 / 5,562,176 / 22,248,704 for all
eight CPU/CUDA method paths. At 16× observations the slowest measured run was
1.17 s and the largest peak RSS was 246.7 MiB, versus 0.91 s and 170.7 MiB at
1×. Thus this measured range shows no time or memory regression and no accuracy
loss; the retained asymptotic algorithms and separate 46k/50k/75k/100k boundary
tests cover sizes that cannot be exercised economically in routine end-to-end
runs.
