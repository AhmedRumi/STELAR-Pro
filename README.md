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

JDK 21 or newer and Python 3 are required. CUDA is optional. The polytomy
resolver is self-contained and requires no Python packages.

```bash
./build.sh
./stelarx -i unrooted_gene_trees.tre -o species_tree.tre --cpu
```

Input contains one unrooted Newick gene tree per non-empty line. STELAR-Pro
temporarily uniquifies repeated leaves, arbitrarily resolves polytomies, restores
the repeated species labels, then roots and tags every gene tree before inference.
Tag-only mode skips polytomy resolution and uniquification.

## STELAR-Pro rooting and tagging

STELAR-Pro accepts unrooted, multi-copy gene trees. Normal inference resolves
polytomies before rooting and tagging. To root/tag without resolving and exit:

```bash
./stelarx -T -i unrooted_multicopy_gene_trees.tre \
  -o rooted_tagged_gene_trees.tre
```

ASTRAL-Pro3 writes `D` on duplication nodes and leaves speciation nodes unlabeled.
Use `--gene-species-map FILE` when gene-copy labels require an explicit two-column
gene-to-species mapping. `--astral-pro-executable FILE` overrides the bundled
`ASTER-Linux/bin/astral-pro3`. Tag-only mode suppresses backend messages and emits
only brief STELAR-Pro status lines. Normal inference currently stops explicitly
before weight calculation: the S1 candidate DP is duplicate-aware, while the
per-tree multi-copy index maps and scoring path are the next implementation stage.

SimPhy datasets default to `$PHYLOGENY_DATA_DIR/simphy/data`. The simulation,
testing, bulk-transfer, and statistics scripts create that directory when it is
missing. Their explicit SimPhy data-directory options still override the
environment-based default.

To remove that complete directory—including simulated datasets and every
inferred result beneath it—preview or run the dedicated cleanup command:

```bash
./clear-bulk-simulated.sh --dry-run
./clear-bulk-simulated.sh --yes
```

The current STELAR-Pro implementation uses the S1 search path and I1
smaller-side traversal. S2/S3 and I2/I3/I4 remain legacy STELAR-X code paths and
are rejected until their duplicate-aware STELAR-Pro versions are implemented.

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

Run the focused tests for the implemented STELAR-Pro stages with:

```bash
test/run_stelar_pro_tests.sh
```

The `run_stelarx_*` suites remain useful for low-level regression checks, but
their STELAR-X topology/score expectations are not STELAR-Pro correctness
oracles. The original migration-focused suite is:

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
