# STELAR-X

STELAR-X is a scalable rooted species-tree summary method. It maximizes agreement
with the rooted triplets displayed by rooted gene trees and provides compact
multi-seed hashes, range intersections, cross-tree recombination, configurable
search spaces and intersections, similarity/UPGMA guidance, parallel CPU paths,
and CUDA acceleration for binary and polytomous rooted inputs.

The checkout directory name and location are arbitrary; scripts resolve the
project root from their own location. The program and artifacts are named
STELAR-X regardless of the checkout name.
The source, Java package, JNI symbols, native libraries, launchers, and portable
artifacts all use the `stelarx` name.

## Build and run

JDK 21 or newer is required. CUDA is optional.

```bash
./build.sh
./stelarx -i rooted_gene_trees.tre -o species_tree.tre --cpu
```

Input contains one rooted Newick tree per non-empty line. STELAR-X uses the
supplied top-level root exactly. A top-level node with anything other than two
children is rejected because ordinary Newick has no independent rootedness flag
and STELAR-X never invents an arbitrary root.

SimPhy datasets default to `$PHYLOGENY_DATA_DIR/simphy/data`. The simulation,
testing, bulk-transfer, and statistics scripts create that directory when it is
missing. Their explicit SimPhy data-directory options still override the
environment-based default.

Use `--search-space S1`, `S2`, or `S3` and `--intersection-method I1` through
`I4`. The default is S1/I2. S2 completes missing taxa with root-preserving
nearest-anchor insertion, adds a UPGMA guide, and enables hash-based cross-tree
recombination. S3 additionally enables consensus-guided enrichment.

Score a supplied rooted species tree with:

```bash
./stelarx -i rooted_gene_trees.tre \
  --score-species-tree rooted_species_tree.tre --cpu
```

The machine-readable result is `TRIPLET_SCORE: N`.

## Crash reports

Unexpected Java failures and JVM fatal-error logs are written under
`crash_logs/`, which STELAR-X creates automatically instead of placing logs in
the repository root. Set `STELARX_CRASH_DIR=/path/to/directory` to override the
location when using the repository launchers.

## CUDA

Build native libraries with `./build_native.sh`. The native libraries are
`libstelarx_weight`, `libstelarx_dp`, `libstelarx_dist`, and `libstelarx_sim`
(with the platform's shared-library suffix). Every intersection method supports
rooted-polytomy weights on both CPU and CUDA.

## Migration details

See [MIGRATION_DOCS/00-overview.md](MIGRATION_DOCS/00-overview.md) for the design,
invariants, retained components, and validation evidence. The canonical names
for every implementation surface are listed in
[MIGRATION_DOCS/04-stelarx-identity.md](MIGRATION_DOCS/04-stelarx-identity.md).

Run the migration-focused validation suite with:

```bash
test/run_stelarx_tests.sh
```

Run the complete layered suite (CPU, independent randomized oracles, malformed
inputs, end-to-end inference, packaging, and CUDA automatically when usable):

```bash
test/run_stelarx_comprehensive_tests.sh
```

Useful variants are `--require-gpu` (CUDA absence is a failure), `--cpu-only`,
`--quick`, and `--skip-packaging`. On a CUDA host, the standalone no-fallback
hardware layer remains available with:

```bash
test/run_stelarx_gpu_tests.sh
```

For accuracy plus wall-time/peak-RSS scaling measurements across the complete
S1–S3 × I1–I4 matrix, run:

```bash
python3 test/test_stelarx_scalability.py --require-gpu
```

Pass `--reference-dir PATH` to compare median S1/I2 CPU resources against a
separately built reference checkout with guarded time and memory ratios.
