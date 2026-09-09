# STELAR-Pro

STELAR-Pro is a scalable rooted species-tree summary method. It maximizes agreement
with the rooted triplets displayed by rooted gene trees and provides compact
multi-seed hashes, range intersections, cross-tree recombination, configurable
search spaces and intersections, similarity/UPGMA guidance, parallel CPU paths,
and CUDA acceleration for binary and polytomous rooted inputs. S1 and the
duplicate-aware intersection implementation are built-in defaults.

The checkout directory name and location are arbitrary; scripts resolve the
project root from their own location. The program and artifacts are named
STELAR-Pro regardless of the checkout name.
The source, Java package, JNI symbols, native libraries, launchers, and portable
artifacts all use the `stelar-pro` name.

## Build and run

JDK 21 or newer and Python 3 are required. CUDA is optional. The polytomy
resolver is self-contained and requires no Python packages.

```bash
./build.sh
./stelar-pro -i unrooted_gene_trees.tre -o species_tree.tre --cpu
```

Input contains one unrooted Newick gene tree per non-empty line. STELAR-Pro
temporarily uniquifies repeated leaves, arbitrarily resolves polytomies, restores
the repeated species labels, then roots and tags every gene tree before inference.
Tag-only mode skips polytomy resolution and uniquification.

## STELAR-Pro rooting and tagging

STELAR-Pro accepts unrooted, multi-copy gene trees. Normal inference resolves
polytomies before rooting and tagging. To root/tag without resolving and exit:

```bash
./stelar-pro -T -i unrooted_multicopy_gene_trees.tre \
  -o rooted_tagged_gene_trees.tre
```

ASTRAL-Pro3 writes `D` on duplication nodes and leaves speciation nodes unlabeled.
Use `--gene-species-map FILE` when gene-copy labels require an explicit two-column
gene-to-species mapping. `--astral-pro-executable FILE` overrides the bundled
`ASTER-Linux/bin/astral-pro3`, which is compiled for the current machine on first
use (see [Machine-specific backends](#machine-specific-backends)). Tag-only mode
suppresses backend messages and emits only brief STELAR-Pro status lines; if the
backend fails, its own output is included in the error. S1 subtree/partition hashing, candidate DP,
and CPU/CUDA intersection indexing are duplicate-aware. Each tree stores a
sorted position vector for every species, so repeated copies count once.

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

The current STELAR-Pro implementation uses S1 and its built-in smaller-side
intersection path by default; neither needs a command-line option. S2 and S3
are reserved names and are rejected until their STELAR-Pro implementations are
ready. The old intersection-selector options have been removed.

Score a supplied rooted species tree with:

```bash
./stelar-pro -i rooted_gene_trees.tre \
  --score-species-tree rooted_species_tree.tre --cpu
```

The machine-readable result is `TRIPLET_SCORE: N`.

## Crash reports

Unexpected Java failures and JVM fatal-error logs are written under
`crash_logs/`, which STELAR-Pro creates automatically instead of placing logs in
the repository root. Set `STELAR_PRO_CRASH_DIR=/path/to/directory` to override the
location when using the repository launchers.

## Machine-specific backends

Two compiled components are **not tracked by git** because a build from one
machine generally does not run on another:

| Component | Built with | Why it does not travel |
|---|---|---|
| `ASTER-Linux/bin/astral-pro3` (rooting/tagging, required) | `g++ -march=native` | glibc symbol versions (`version 'GLIBC_2.38' not found`, exit status 1 before `main`), CPU instruction set (SIGILL) |
| `native/libstelar_pro_{weight,dp,dist,sim}.so` (CUDA, optional) | `nvcc` | glibc symbol versions, GPU compute capability below the artifact's minimum |

`./run.sh` (and therefore `./stelar-pro` and the monitor wrapper) calls
`./ensure_backends.sh` before every run. It is a no-op when both components were
built on this machine from the current sources; otherwise it rebuilds what is
missing, foreign, or stale and records a build stamp (`ASTER-Linux/bin/.build-stamp`,
`native/.build-stamp`) describing the glibc, CPU, and GPUs it was built for.

```bash
./ensure_backends.sh            # build/refresh both for this machine (what run.sh does)
./ensure_backends.sh --force    # rebuild everything
./ensure_backends.sh --cpu-only # skip the CUDA libraries
./stelar-pro --diagnose         # shows whether astral-pro3 and CUDA are usable here
```

Requirements: `g++` and `make` for ASTRAL-Pro3 (`sudo apt install build-essential`);
`nvcc` for the CUDA libraries, which are skipped with a CPU fallback when nvcc or
an NVIDIA GPU is absent. `CUDA_ARCH` overrides the nvcc target (default `native`,
i.e. exactly this machine's GPU; `all-major` for a build meant to run on other
GPUs, as `build_portable.sh` uses). `--no-build` skips both the Java build and the
backend check. To use an ASTRAL-Pro3 built elsewhere, pass
`--astral-pro-executable FILE` or set `STELAR_PRO_EXECUTABLE`.

The CUDA libraries can also be built directly with `./build_native.sh`; they are
`libstelar_pro_weight`, `libstelar_pro_dp`, `libstelar_pro_dist`, and
`libstelar_pro_sim` (with the platform's shared-library suffix). The built-in
intersection path supports rooted-polytomy weights on both CPU and CUDA.

## Migration details

See [MIGRATION_DOCS/00-overview.md](MIGRATION_DOCS/00-overview.md) for the design,
invariants, retained components, and validation evidence. The canonical names
for every implementation surface are listed in
[MIGRATION_DOCS/04-stelar-pro-identity.md](MIGRATION_DOCS/04-stelar-pro-identity.md).

Run the focused tests for the implemented STELAR-Pro stages with:

```bash
test/run_stelar_pro_tests.sh
```

The `run_stelar_pro_*` suites remain useful for low-level regression checks, but
their STELAR-Pro topology/score expectations are not STELAR-Pro correctness
oracles. The original migration-focused suite is:

```bash
test/run_stelar_pro_tests.sh
```

Run the complete layered suite (CPU, independent randomized oracles, malformed
inputs, end-to-end inference, packaging, and CUDA automatically when usable):

```bash
test/run_stelar_pro_comprehensive_tests.sh
```

Useful variants are `--require-gpu` (CUDA absence is a failure), `--cpu-only`,
`--quick`, and `--skip-packaging`. On a CUDA host, the standalone no-fallback
hardware layer remains available with:

```bash
test/run_stelar_pro_gpu_tests.sh
```

For accuracy plus wall-time/peak-RSS scaling measurements of the current
default implementation, run:

```bash
python3 test/test_stelar_pro_scalability.py --require-gpu
```

Pass `--reference-dir PATH` to compare median default-path CPU resources against a
separately built reference checkout with guarded time and memory ratios.
