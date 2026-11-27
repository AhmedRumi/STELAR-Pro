# STELAR-Pro C++17

Separate native implementation of the three retained production search modes.
The Java sources, Java launchers, previous temporary builds and original JNI
sources are unchanged. CPU and CUDA inference are implemented here, not by
launching Java.

| C++ option | Java reference | Candidate construction |
|---|---|---|
| `--search-space s1` | S1 | Original speciation-rooted candidates and local transitions |
| `--search-space s2` | S4 | Strict-lineage DISCO completion, triplet similarity and UPGMA; no S1 or cross-tree transitions by default |
| `--search-space s3` | S6 | Subset-local, duplication-scope DISCO backbones, triplet placements, group UPGMA and triplet NNI; four guides by default |

`s3` dispatches to the complete `s2` pipeline at **100 or fewer species**.
The overall CLI default remains `s1`. Java S2/S3/S5 and experimental completion
variants were not brought into production under other names.

## Build and run

From the repository root, using GCC/Clang with C++17 and POSIX threads:

```bash
make -C cpp -j4 cpu
make -C cpp -j4 gpu CUDA_ARCH=sm_89

cpp/bin/stelar-pro-cpp-cpu -i genes.trees -o species.tree --search-space s2 --cpu -t 28
cpp/bin/stelar-pro-cpp -i genes.trees -o species.tree --search-space s3 --gpu-strict -t 28
```

The CUDA build defaults to `/usr/local/cuda/bin/nvcc`; set `NVCC` to override it.
Without `CUDA_ARCH`, Make uses `all-major` for a portable multi-architecture
build. `sm_89` above is the GPU tested on this machine, not a requirement for all
users. Clean/rebuild when changing architecture or compiler flags.

CMake 3.24 or newer is also supported:

```bash
cmake -S cpp -B cpp/build/cmake-cpu -DSTELAR_CUDA=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build/cmake-cpu -j4
ctest --test-dir cpp/build/cmake-cpu --output-on-failure

cmake -S cpp -B cpp/build/cmake-gpu -DSTELAR_CUDA=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build cpp/build/cmake-gpu -j4
ctest --test-dir cpp/build/cmake-gpu --output-on-failure
```

Both Make executables have been built locally. The CPU executable has no CUDA
dependency. `--auto` chooses CUDA when available; `--gpu` permits CPU fallback;
`--gpu-strict` fails on unavailable CUDA or a failing CUDA operation. Some small
matrices and NNI deliberately remain on CPU even in strict GPU mode.

### Existing preprocessing dependencies

As in Java, normal inference calls the repository's Python polytomy/copy-label
preprocessor and the ASTRAL-Pro3 C++ executable for rooting/tagging. Keep the
repository alongside the executable; the build records its source location.
Set `STELAR_PRO_HOME` after relocating the repository. Backend overrides are
`--astral-pro-executable`, `STELAR_PRO_EXECUTABLE` and `STELAR_PRO_PYTHON`.
No JDK is needed for normal native inference.

**Input labels denote species.** Do not feed simulated copy identifiers such
as `61_12_0` as distinct species if the intended species is `61`. Use the same
species-mapped input as the previous experiments. As in Java, `--gene-species-map`
is available for tag-only output; direct mapping ingestion for inference/scoring
is rejected. Native parsing intentionally does not guess which
underscores belong to a species name. Branch lengths are ignored by the
triplet objective. Quoted labels and nested Newick comments are supported.

## Options and compatibility

Run `cpp/bin/stelar-pro-cpp --help` for the full accepted list. Principal flags:

| Group | Options |
|---|---|
| Compute | `--cpu`, `--gpu`, `--auto`, `--gpu-strict`, `-t` / `--threads` |
| Independent cross-tree override | `--cross-tree-splits`, `--no-cross-tree-splits`, `--search-mode local\|full` |
| Exhaustive s2 | `--exhaustive`: add S1 and complete DISCO components globally; usable with or without cross-tree splits |
| s3 | `--s3-rounds`, `--s3-seed`; old `--s6-rounds` / `--s6-seed` aliases retained |
| Scoring/restriction | `-c` / `--score-species-tree`, `--taxa-file`, `--prune-search-space`, `--no-prune-search-space` |
| Utilities | `-T` / `--tag-only`, `--extract-taxa`, `--taxa-set union\|intersection`, `--diagnose` |
| GPU memory/batching | `--gpu-batch-size`, `--gpu-batches`, `--no-gpu-batch`, `--gpu-vram-occupancy-factor`, `--gpu-vram-control-factor`, `--gpu-treewalk-vram-cap-mb`, `--gpu-sim-vram-cap-mb`, `--gpu-dist-tile-size` |
| GPU reporting/cross-tree output | `--gpu-progress-interval`, `--gpu-dp-state-space-construction-output-cap`, `--gpu-dp-state-space-progress-time-interval`, `--gpu-dp-state-space-progress-max-steps` |
| Output/debugging | `--log-file`, `--dump-clusters`, `--dump-completed-gene-trees`, `--verify-parse`, `--verify-hash`, `--verify-clusters`, `--verify-partitions`, `--verify-dp`, `--verify-weights`, matrix/UPGMA verifiers, verbosity/quiet flags |

As in the selected Java modes, cross-tree transitions are not enabled for s1;
`--exhaustive` belongs to s2. Search overrides are independent of the placement
of `--search-space` on the command line; the last explicit cross-tree override
wins. Only one outer search pass is supported. Removed S6 aliases
`--subsetlocal` and `--s6-no-upgma-candidate` remain removed.

Compatibility details that should not be mistaken for new algorithms:

- `--completion-method distance` and Step-B settings remain accepted legacy
  selectors; they do not replace the selected modes' active triplet pipeline.
- `--seeds` is retained. Native set equality checks exact species masks;
  randomized hash collisions cannot merge different candidate sets.
- Native accumulation is exact: CPU uses unsigned 128-bit scores; CUDA chooses
  safe 64-bit or two-limb 128-bit accumulation. `--large-n-score-type int128`
  forces the wide CUDA path. `double` is accepted, with an explicit notice,
  but does **not** introduce approximate floating-point scores in this port.
- `--rooted`, `--keep-polytomy-during-inference`, scoring and output restriction
  flags are retained. Unsupported unrooted/anchored, autocomplete and consensus
  enrichment combinations still produce errors instead of changing the model.
- `--xms`, `--xmx`, `--no-build` and `--no-notify` are accepted launcher
  compatibility options; they do not configure a native allocator or JVM.
- Native hash/cluster/partition/DP/weight dumps use exact-set TSV schemas, not
  Java object/hash serialization. Weights are written to `<output>.weights.tsv`.
  s2 completed-tree dumps stream in component order without retaining graphs.
- S1 distance/similarity/UPGMA diagnostics match the Java CPU text format.
  s2 diagnostics expose its DISCO triplet matrix; s3's matrix diagnostic shows
  observed DISCO evidence, **not** a global matrix used by its subset-local
  search. Distance-only diagnostics run on CPU.
- `--verify-greedy-consensus` preserves the original, consensus-only verifier
  by compiling/running the unchanged Java reference in a private temporary
  directory. It needs a JDK and is explicitly announced. It does not add
  consensus to native inference; no Java or JNI is used by s1/s2/s3 inference.
- `--tagged-input` is an additional diagnostic-only shortcut for already
  rooted/event-tagged trees; ordinary datasets should use normal preprocessing.

Inference writes a tree and `<output>.stats.json` with candidates, scored
splits, original partitions, exact scores, elapsed time and peak process RSS.
RSS in this JSON is the native process's high-water mark, not a sum of all
subprocesses. Use `/usr/bin/time -v` for whole-command resource comparisons.
Input/output aliases, hard links and derived-artifact collisions are rejected.

## Memory and parallelism

- Trees are flat node arrays (28-byte nodes here), child/sibling indices and
  leaf orderings. Raw DISCO components do not allocate a global-sized position
  map or CSR index per component. Scopes are trimmed missing-species bitsets.
- Exact species sets occupy a shared contiguous word arena with compact
  records and open-addressed interning. Candidate transitions are 12-byte ID
  triples with a compact deduplication table and CSR dependencies.
- Original rooted/tagged genes alone supply scoring observations. Identical
  observed partitions retain their frequency; identical candidates score once.
- Independent families, pair-owned matrix rows, completion batches and CPU
  split scores use a persistent bounded worker pool. Nested work uses a serial
  inner path; exceptions propagate to the caller without pool starvation.
- s2 completes bounded component batches, interns their candidates and destroys
  their node arrays immediately. No completed-tree exemplar forest survives.
- s3 reuses its compact **uncompleted** DISCO pool across rounds. Backbone
  temporaries die before placements. Partition matrices accumulate directly in
  one component traversal, without materializing every restricted group tree.
  All guide-support state dies before final inference.
- Similarity uses a packed double triangle. Completion ranks currently use
  32-bit IDs (4n² bytes), not Java's narrower rank representation. This is an
  explicit remaining memory opportunity at very large n.
- CUDA owns all buffers through RAII. Original species bitsets are shared across
  batches; weights use one cooperative block per split with parallel partition
  reduction. Exact partial scores are staged when observations exceed the
  memory budget. No split-by-partition result matrix is allocated.
- Focused s3 support uses selected-backbone slices, resident or bounded staged
  observations, and reused scratch. No n³ triplet tensor is allocated. CUDA
  cross-tree lookup confirms exact subset/residual masks after hash filtering.
  CUDA allocations, launches and copies are checked; buffers are released by
  their owning operation/state destructor.

Weight complexity still depends on candidate count, observed partitions and
bitset width. Removing JNI alone cannot guarantee a speedup on every dataset.
See [VALIDATION.md](VALIDATION.md) for measured evidence and current limits.

## Tests

```bash
make -C cpp -j4 test
make -C cpp -j4 test-gpu CUDA_ARCH=sm_89
python3 cpp/test/parity.py --native cpp/bin/stelar-pro-cpp-cpu \
  --gpu cpp/bin/stelar-pro-cpp --output cpp/validation/cli
python3 cpp/test/check_cli_utilities.py --native cpp/bin/stelar-pro-cpp \
  --java-classes cpp/validation/cli/java-classes --work cpp/validation/utilities --gpu
python3 cpp/test/check_s6_reference.py --native cpp/bin/stelar-pro-cpp \
  --java-classes cpp/validation/cli/java-classes --work cpp/validation/s3 --gpu
```

Java comparison helpers compile into test output directories, never over the
production Java build. Test fixtures, commands, logs and JSON results are
retained by the integration harnesses. Executables and generated validation
artifacts are local build products and excluded from Git.
