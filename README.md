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

## Run ASTRAL-Pro3 as a comparison method

The bundled native ASTRAL-Pro3 can be run on any Newick gene-tree file through
the monitored wrapper:

```bash
./run-astral-pro3-with-monitor.sh \
  -i path/to/all_gt.tre -o results/out-astral-pro3.tre \
  --threads 16 --no-notify
```

Use `--mapping gene_to_species.map` when copy labels differ from species labels,
or pass additional native options with `--opts`. The wrapper records wall time,
maximum CPU/GPU memory, optional RF rate, and exit status next to the output in
`out-astral-pro3_stats.csv`.

ASTRAL-Pro3 is also a first-class method in both dataset sweep scripts:

```bash
./run-bulk-simulated.sh --method astral-pro3 \
  --taxa-list "100,200" --genes-list "100,1000" \
  --num-replicates 10 --opts "--thread 16 --seed 42"

./run-bulk-standard.sh --method "stelar-pro;astral-pro3" \
  --folder "37-taxon" --astral-pro3-opts "--thread 16 --seed 42"
```

Simulated outputs are stored below each replicate in
`astral-pro3-outputs/<setting>/`; `collect-stats-simulated.sh` merges them with
STELAR-Pro results. The aliases `astral-pro` and `apro3` are accepted for method
selection, while `astral` continues to mean the separate Java ASTRAL baseline.
`--search-space S1` is specific to STELAR-Pro and must not be supplied to
ASTRAL-Pro3. Use ASTRAL-Pro3's `--round`/`--subsample` options (or `-R`) when a
larger ASTRAL-Pro3 search is desired.

Both methods also maintain a lightweight safety copy at
`$PHYLOGENY_DATA_DIR/outputs/gdl-simulation/<method>-outputs/<dataset>/<R>/<setting>`.
It contains inference results and reproducibility command records, but never
the large simulated gene trees or SimPhy databases. Use
`sync-simulated-outputs.sh` to backfill existing runs. The complete layout and
safety rules are documented in
[DOCS/simulated-outputs-mirror.md](DOCS/simulated-outputs-mirror.md).

## Run actual GDL simulated datasets

Use the separate, inference-only sweep for existing data under
`$PHYLOGENY_DATA_DIR/gdl-simulation/data/taxa<T>_gt<G>_dup<D>_loss<L>_pop<P>/R*/`:

```bash
./run-bulk-gdl-simulated.sh --method stelar-pro \
  --taxa-list "10,100" --gt-list "10,1000" \
  --dup-list "1,2" --loss-list "0,0.5" --pop-list "50000000" \
  --num-replicates 5 --opts-list "--search-space S1" --no-notify

./run-bulk-gdl-simulated.sh --method astral-pro3 \
  --num-replicates 1 --opts "--thread 16" --no-notify
```

Only `true-genetrees/all_gt.trees` and `species-tree/s_tree.trees` are read;
`separate/` is never needed. Gene-copy labels such as `10_0_0` are normalized to
species `10` in temporary input for both methods, without changing source data.
Editable parameter lists are at the top of the new bulk script. Missing inputs
are skipped by default; `--strict-missing` aborts before inference, and
`--dry-run` previews without writes. Notifications are enabled unless disabled.
See [GDL sweep details](DOCS/gdl-simulated-runs.md) for running both methods,
output safety, stats collection, and the local sample command.

For the original SimPhy workflow, the dedicated cleanup removes `simphy/data`,
including its simulated inputs and inference results. It leaves actual GDL
datasets and the outputs safety mirror untouched:

```bash
./clear-bulk-simulated.sh --dry-run
./clear-bulk-simulated.sh --yes
```

STELAR-Pro uses S1 and its built-in smaller-side intersection path by default.
S2 and S4 retain the original rooted/tagged multicopy trees for scoring, but
their default candidate search uses only lineage-scoped completed DISCO
single-copy trees and a UPGMA guide—not S1 candidates. A component is completed
only with taxa present below its root before DISCO cuts. Use
`--exhaustive` with S2 or S4 to restore S1 candidates and completion
of every DISCO tree to all observed species; this does not turn on cross-tree
transitions.
Cross-tree splits are disabled by default in S2, S4, S5, S6, and S7; enable them with
`--cross-tree-splits`. For S5/S6/S7 this also controls their subset-backbone DP.
Normal S2/S4 inference completes DISCO trees in bounded
CPU-parallel batches, merges their candidates and local transitions, and drops
their node graphs; only compact indexes for novel candidate exemplars remain.
S3 retains its existing global DISCO completion plus S1 candidates and builds
greedy consensus trees from the completed DISCO trees, adding their clades and
polytomy refinements. S4 uses S2's bounded candidate pipeline, but
builds a rooted-triplet sister-frequency matrix for completion and UPGMA;
the former S4 consensus stage is dormant.
Its distance view is derived as `1 − similarity`. Select
them with `--search-space S2`, `S3`, or `S4`. S5 is independent of those
DISCO/completion/consensus paths: it adds four ASTRAL-Pro-style quartet-greedy
species-tree proposals to S1, registers both rooted orientations of each
greedy-tree edge, and uses STELAR-Pro's rooted-triplet DP to choose candidates.
Greedy placement/NNI uses parallel CPU workers; the final weight calculation
remains STELAR-Pro's CPU/GPU triplet scorer. The original rooted/tagged trees
are the only scoring evidence in all modes. Select S5 with `--search-space S5`.
Each DP pass adds four greedy trees by default (`--s5-rounds N` changes this
batch size); `--s5-seed N` makes the pass sequence reproducible. With fewer than 100
taxa, S5 uses direct full-taxon greedy placement. At 100 taxa or more, each
initial round uses ASTRAL-Pro's two-step search: build randomized subset guide trees,
infer their backbone using STELAR-Pro's rooted-triplet DP on the original
tagged gene trees restricted to that subset, then insert the remaining taxa
in groups and run NNI. The default is one search/DP pass. Set
`--max-search-passes N` to permit additional passes: each induces a new random
backbone from the preceding STELAR tree, adds another batch of greedy proposals,
and reruns the DP. All S1 and earlier S5
candidates remain available. Search stops when the rooted tree no longer
changes, following the ASTRAL-Pro2 supplement's tree-stability rule; the
bundled ASTRAL-Pro3 binary instead uses score improvement and a different cap.
Every backbone and full-tree DP remains rooted-triplet based.

S6 is independent of S2–S5: select `--search-space S6`. It keeps every S1
candidate and adds four rooted species-tree proposals per DP pass by default
(`--s6-rounds N`, `--s6-seed N`). Unlike S5, **all greedy placement and rooted
NNI decisions use STELAR-Pro's speciation-driven triplet score** from the
original tagged multicopy gene trees. S6 stores speciation-child taxon lists
and incremental colour counts, not a dense triplet tensor; duplicate copies
are collapsed within each child species set, and duplication-rooted nodes
contribute no observations.
For 100 or more taxa, S6 builds randomized subset guides, infers a rooted
STELAR-DP backbone, groups the remaining taxa by their preferred backbone
edge, then places and refines them by triplet score. With
`--max-search-passes N` greater than one, later passes restrict the previous
rooted STELAR tree to a fresh subset, accumulate all S1 and S6 candidates,
and stop when the rooted DP topology is unchanged or the pass limit is reached.
S6 does not use ASTRAL-Pro's quartet placement counter or
DISCO/completion/consensus; ASTRAL-Pro3 is used only for the common input
rooting/tagging preprocessing. CPU workers score independent guides and NNI
choices; CUDA builds bounded O(n²) slices for taxa being inserted. The final
STELAR weight calculation remains CUDA-capable. Guide-slice reuse is capped at
32 MiB, and one-pass runs release their search-only indexes and native pages
before final DP.
The subset guides score only gene families induced on their selected taxa;
their backbone DP takes candidates from those guides and weights from the
induced original genes, as in S5. Native CPU colour states accelerate the
triplet placement and NNI queries, with a Java fallback when unavailable.

S7 combines S5's A-Pro quartet-greedy proposals with S2's DISCO-completed
single-copy proposals and UPGMA guide (`--search-space S7`). It retains S1
clusters, both rooted orientations of greedy-tree edges, and local transitions
from completed trees. Normal inference streams completion in bounded batches;
the completed node graphs are released after their candidates enter X and DP.
Only the original rooted/tagged multicopy gene trees contribute triplet weights.
S7 has its own `--s7-rounds N` and `--s7-seed N` settings. It runs one DP pass;
cross-tree transitions are optional through `--cross-tree-splits`.

The 1-lb, ten-replicate comparison can be reproduced with
`python3 scripts/benchmark_s6_dupmean1.py`; it writes a resumable per-replicate
run directory and `pro_data/STELAR_PRO_S6_VS_S5_ASTRAL_PRO3_DUPMEAN_1_GPU_SUMMARY.md`.
The old intersection-selector options have been removed.

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
