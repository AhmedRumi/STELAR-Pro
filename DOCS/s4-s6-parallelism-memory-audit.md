# S4 and S6: production defaults, parallelism, and memory audit

Reviewed against the source on 2026-11-19. The changes in this task promote the
tested S6 variant and preserve the previous S6. The performance opportunities
below are proposals; they have not been silently substituted into the algorithm.

## Current choices and reproducibility

| Choice | Candidate construction | Final DP candidates | Scoring observations |
|---|---|---|---|
| `--search-space S4` | Strict lineage-scoped DISCO completion, one global triplet similarity matrix, UPGMA | Scoped completed DISCO clades/transitions and UPGMA; no S1 by default | Original rooted/tagged multicopy genes |
| `--search-space S6` | Subset-local DISCO completion with duplication scope, triplet-DP backbones, triplet placements, local UPGMA partition trees, triplet NNI | Four rooted greedy guides by default; no standalone global UPGMA or S1 candidates | Original rooted/tagged multicopy genes; backbone DP uses their restrictions |
| Temporary **S6 without subset-local** | Previous strict-lineage global DISCO guide pool, shared global triplet matrix and UPGMA | Previous greedy guides plus global UPGMA by default | Same original-gene model |

S6 needs no `--subsetlocal` flag now. That flag remains an accepted S6-only
compatibility alias. `--s6-no-upgma-candidate` is also an accepted alias because
the final S6 candidate list already contains only its greedy guides. UPGMA is
still used inside guide construction. Cross-tree splits remain off by default.
S4 retains its existing `--exhaustive` option; its normal strict scope is unchanged.
The top-level search-space default remains S1; select S4 or S6 explicitly.

For a component detached from duplication vertex v, new S6 saves
`observed species below v minus species remaining in the component`.
It computes the union bottom-up before applying any cuts. The residual top
component retains the original gene's observed species scope. The masks are
packed and trimmed; no per-node scope matrix is retained. Each younger
duplication uses its own observed union. There is no inherited top-down
completion scope and no restoration of species absent from the whole gene.

The previous S6's global guide preparation was removed from production
`S6S2BBGuideData`. Its jar, relevant pre-change source snapshots and launcher
are under [temporary-builds/s6-without-subset-local](../temporary-builds/s6-without-subset-local/README.md).
That launcher uses a frozen Java jar and the repository's native backends.

Validation: the new default was replayed on estimated-gene replicates 01–10 of
`sim200_dup1_loss1`. All ten match the saved `S6-subset-duplication` experiments
in unrooted output topology, exact final triplet score, and the number of
scored splits in **every** backbone/final weight call. An S4-strict control on
replicate 01 also matches all three checks. Scope tests check sister-locus
eligibility, exclusion of ancestral siblings and absent species, nested
duplications, word boundaries, input immutability and unchanged S4 scope.
The replay harness is [test/check_s6_promotion_parity.py](../test/check_s6_promotion_parity.py).
Saved evidence and the broader suite's unrelated missing-`dendropy` limitation
are in [the validation record](../temporary-builds/s6-promotion-validation-20261119/README.md).

## What the existing measurements establish

The following phase numbers come from the already completed, local 2000-taxon
R4 logs: S4 from 2026-11-19 and the promoted S6 variant from 2026-11-18.
These are reference measurements, not a new controlled speed experiment.
RAM is the process's live RSS sampled during the phase; it includes retained
objects from preceding phases. The RAM column is therefore **not** memory
allocated exclusively by that phase and must not be added across rows.

| Phase | S4 seconds | S4 peak RAM, MiB | S6 subset-local duplication seconds | S6 peak RAM, MiB |
|---|---:|---:|---:|---:|
| Parse tagged genes | 1.027 | 1,162 | 1.506 | 1,300 |
| Candidate/guide preparation, Phase 1b | 15.756 | 4,726 | 188.936 | 5,762 |
| Prefix hashing, Phase 2 | 0.049 | 4,726 | 0.029 | 4,395 |
| Candidate hashes/clusters, Phase 3 | 1.816 | 5,640 | 0.008 | 4,013 |
| Original child partitions, Phase 4 | 1.224 | 6,072 | 3.258 | 3,715 |
| Final local transitions, Phase 5 | 2.519 | 6,148 | 0.001 | 3,715 |
| Final weights, Phase 6 | 2,450.423 | 6,700 | 41.358 | 3,808 |
| Final inference, Phase 7 | 0.321 | 6,312 | 0.007 | 3,412 |

Common external resolution/rooting/tagging and between-phase work also belong
to total runtime; the table is not a complete decomposition of total elapsed
time. S6 Phase 1b includes its backbone weight calls and NNI. Its final weight
row does not include those earlier backbone calls.

S4 scored 601,946 final splits in this run. S6 scored 2,074 final splits, plus
four small backbone DP weight calls. This difference explains why their best
optimization targets differ: S4 is dominated by final GPU weights, whereas
the S6 construction phase accounts for much more of its total runtime.

## Phase-by-phase implementation

### 0–1. Preprocessing and original gene parsing

Both modes resolve/uniquify input through the Python preprocessor, then invoke
ASTRAL-Pro3 `-T` for rooting/tagging. `GeneTreeRooterTagger.run` launches one
external process and does not forward the configured Java `--threads` value.
Setting 28 Java threads consequently does not imply that this external stage
uses 28 workers. Afterward `TreeParser` parses individual tagged trees on the
shared fixed CPU pool, with independent output slots. CUDA is not used here.

Each original `Tree` retains its node graph, species ordering, a representative
`positionMap`, and a CSR index of **all** positions of each species. The latter
is necessary for duplicate-invariant range membership. Original trees remain
available for weight scoring. Temporary resolved/tagged files are removed in
`Main`'s `finally` block.

Opportunities: time the external steps separately; parallelize independent
families in bounded chunks if rooting/tagging semantics and deterministic
output order can be preserved. Avoid parsing the whole input repeatedly for
S6 backbones by restricting the original in-memory gene forest while preserving
event tags exactly. Keep a parity test for restrictions before replacing the
current parser-based restriction.

Code: [Main](../src/stelarx/Main.java),
[GeneTreeRooterTagger](../src/stelarx/pro/GeneTreeRooterTagger.java),
[TreeParser](../src/stelarx/tree/TreeParser.java),
[Tree](../src/stelarx/tree/Tree.java),
[TaxonPositionIndex](../src/stelarx/tree/TaxonPositionIndex.java).

### 1b-A. DISCO decomposition and scope storage

`DiscoDecomposer.decomposeAll` currently visits gene families sequentially.
Within one copied family, a bottom-up walk merges smaller species sets into
larger sets and removes consumed child entries from the live map. Cut decisions
and scope masks are saved before cuts mutate the topology. Unary duplication
nodes are suppressed, and components below four observed leaves are filtered
as before. The input scoring tree is never cut.

Working set merges are transient. Only component masks survive; mask storage
is at most O(D ceil(n/64)) words for D retained components, with empty/trailing
words removed. However, the retained components are still full `Tree` objects:
each has an O(n) representative map **and** O(n) CSR offsets, even if it has
very few leaves. S4 retains this pool through similarity construction and
then releases batches; S6 retains the uncompleted pool across its four rounds.
S6 must retain or compact it because each round restricts it differently.

Opportunities: parallelize independent families with a bounded number of
workers, save each family's output in its original slot, then assign final
indices in input order. Unbounded parallel decomposition would multiply live
copies and transient sets. The largest memory opportunity is a lighter raw
DISCO component representation: topology, observed taxon list and missing
mask, with dense indexes created only when a component is restricted/consumed.
At n=2,001, 65,942 components' two dense n-sized integer indexes alone imply
about 1.0 GiB of payload; this arithmetic excludes nodes and leaf arrays.

Code: [DiscoDecomposer](../src/stelarx/pro/DiscoDecomposer.java),
[S6S2BBGuideData.prepare](../src/stelarx/pro/S6S2BBGuideData.java).

### 1b-B. Triplet similarity and completion anchors

S4 uses one matrix over all species. For a rooted single-copy component t with
k leaves, pair (a,b) contributes numerator `k - size(LCA_t(a,b))` and
denominator `k - 2` when both species are present. Summing these counts before
normalization gives the triplet sister-frequency similarity. No triplet tensor
is constructed. S6 applies the same calculation to each restricted backbone
or partition group, using only observed component topology as matrix evidence.

CPU tour/RMQ construction is parallel over components within bounded batches.
CPU accumulation is parallel over pair-owning rows; each worker writes distinct
cells. CUDA assigns taxon pairs to threads, loops over batched tree data, and
uses tiled matrix outputs. Compact versus wide blocked RMQ prevents narrow
index overflow. Java flattening and tour buffers are bounded by host batches;
native device buffers and pinned tile output are freed when their JNI call
returns. Contributions are normalized after all batches.

`SimilarityMatrix` stores symmetric pairs in segmented double triangles:
about 8n² bytes while numerator and denominator coexist, and 4n² afterward.
The denominator reference is dropped at normalization. `SortedRowsBuilder`
sorts rows on CPU workers using reusable per-worker primitive scratch.
Anchor ranks occupy 2n² bytes up to 65,536 taxa and 4n² above that.
The matrix is consumed in place by UPGMA when it has no later readers.

In S6, matrices below 64 taxa use the CPU to avoid a separate small CUDA launch.
Partition matrices always use the CPU today. S4 uses bounded CUDA similarity
when GPU mode is available; GPU mode is not evidence that every small operation
executes on the GPU.

Opportunities: balance triangular CPU pair work by pair count rather than equal
row counts; pipeline host tour construction with CUDA consumption using bounded
buffers; batch many small partition matrices into one CPU task/CUDA launch.
Benchmark the small-matrix CPU/GPU threshold. The global S4 matrix is modest
at 2,001 species (about 15.3 MiB after normalization), so optimizing component
objects and indexes is likely more useful there than reducing matrix precision.

Code: [SimilarityMatrixBuilder](../src/stelarx/completion/SimilarityMatrixBuilder.java),
[TripletTourBuilder](../src/stelarx/completion/TripletTourBuilder.java),
[SimilarityMatrix](../src/stelarx/completion/SimilarityMatrix.java),
[SortedRowsBuilder](../src/stelarx/completion/SortedRowsBuilder.java),
[CUDA similarity](../src/native/stelarx_similarity.cu).

### 1b-C/5. S4 scoped completion and candidate streaming

Completion mutates independent copied components in parallel **within** a batch.
Taxa are inserted in deterministic species-ID order, beside the best present
anchor from the precomputed row. This is the existing nearest-anchor procedure.
The normal inference path batches 64–256 components. It then serially hashes
completed nodes, deduplicates clades and local transitions into X/DP, and
retains a compact exemplar only if that tree contributed a new clade.

The exemplar has no node graph: it shares the ordering and position indexes
needed to recover candidate membership. After a batch, source references are
set to null and completed node graphs become collectible. UPGMA is emitted
after all DISCO batches. Completion rows are then dropped. Java reclamation
occurs through GC; nulling a reference does not guarantee an immediate RSS drop.

UPGMA's n−1 merge decisions are sequentially dependent. Within each merge,
the global clusterer parallelizes best-pair reduction, row updates and stale
row rescans. It consumes the matrix triangle instead of copying it when safe.

Opportunities: build per-worker emission buffers, followed by deterministic
deduplication, to overlap completion and extraction. Memory quotas must bound
buffers and simultaneously live batches. Compact or sparsify exemplar indexes;
nulling completed roots alone does not remove their retained O(n) position
maps. Exact duplicate components may be completed once if topology, labels,
scope and anchor ordering all agree; their multiplicities must still count
in similarity evidence. Deduplicating only by observed topology is insufficient.

Code: [ScopedTreeCompleter](../src/stelarx/pro/ScopedTreeCompleter.java),
[S2CandidateStreamer](../src/stelarx/pro/S2CandidateStreamer.java),
[UPGMAClusterer](../src/stelarx/completion/UPGMAClusterer.java).

### 1b-D. S6 backbones, edge placement, partition trees and NNI

The production four-round loop is currently sequential. At n≥100, each round
selects approximately `sqrt(n) log2(n) / 4` backbone species, preserving the
existing size rule and seed. It restricts original genes for weight evidence
and raw DISCO components for candidate generation. Completion intersects the
duplication scope with the selected species. Backbone completion is parallel
over components. A small local triplet matrix/UPGMA supplies candidate guides;
STELAR's triplet DP infers the backbone. Its temporary hash indexes are released
after cluster/partition/transition construction, and backbone scoring buffers
are local to that call.

Triplet scores then assign each remaining species to one edge of the fixed
backbone. These assignments are independent except for deterministic RNG
tie-breaking, but the current assignment loop is serial. The focused CUDA
query returns two m×m count slices for m backbone species. It scans indexed
memberships containing the inserted species, rather than all possible triplets.
CPU code evaluates edge scores from those counts.

For n≤256, immutable focus lists reside on the GPU. For larger n, they reside
in native host vectors; a CUDA query compacts/uploads only relevant partitions
and selected species. This deliberately trades some host work/transfers for
bounded VRAM. Per-query device scratch is freed after copying focused results.
There is no resident n³ tensor. A separate small guide-focus cache, when used,
has a 32 MiB admission cap.

After edge assignment, a single component traversal distributes observed
leaves to disjoint groups using small-to-large maps. No partition completion
is performed. Currently restricted trees for **all** groups are collected,
then each group's CPU similarity matrix is built serially; that group's source
list is nulled afterward. The returned matrices also coexist until their
UPGMA trees are constructed, when matrix references are nulled one by one.
`upgmaGroup` uses its own serial priority queue and a dense `(2m−1)²` double
scratch table; this is distinct from the parallel global `UPGMAClusterer`.

NNI uses cached local decisions and up to three borrowed CPU colour states.
Each native state holds three colour counts per speciation child and updates
only partitions touched by species whose colour changed. JNI does not create
an independent 28-thread pool for each query. Java workers evaluate bounded
windows and commit changes in a deterministic order. More workers need more
mutable state and can disrupt colour locality, so increasing them blindly can
increase memory and reduce speed.

Only completed full greedy guides survive Phase 1b. The original uncompleted
DISCO pool is released with `guideData`; the native CPU/GPU support is explicitly
closed before the final weight phase. The CUDA index destructor calls cudaFree;
the native CPU destructor deletes vectors and may call malloc_trim for large
indexes. Java arrays become collectible, with an explicit GC hint at the phase
boundary. Native allocated pages, Java live heap and OS RSS are distinct metrics.

Opportunities, in likely useful order:

1. Replace raw components' dense indexes as described above; reuse them across
   rounds without regenerating DISCO.
2. Process a bounded group window through restriction → matrix → UPGMA, returning
   only partition trees. This shortens both restricted-tree and matrix lifetimes.
   Preserve the current single distribution pass or use compact per-group tokens
   to avoid repeatedly scanning every component.
3. Batch fixed-backbone taxon assignments into CUDA tiles and reuse selected maps
   and device scratch. Precompute RNG draws/tie choices in input order so parallel
   scheduling does not change reproducibility. Bound memory by batch size × m².
4. For CPU mode, parallelize assignments with a small measured number of borrowed
   states. Avoid nested tasks on the same fixed thread pool.
5. Build independent group matrices/trees concurrently with a memory-based worker
   cap and a serial inner path. Current pair builders submit tasks to the shared
   executor, so simply submitting the outer groups there risks pool starvation.
6. Consider concurrent rounds only after profiling. It multiplies backbone trees,
   matrices and scorer states; bounded assignment/group parallelism is a safer
   first improvement.

For n<100, S6 uses direct triplet insertion/NNI. The current orchestration still
prepares its DISCO pool before that decision; lazily skipping this unused pool
is another small-input opportunity.

Code: [S6S2BBGuideData](../src/stelarx/pro/S6S2BBGuideData.java),
[S6GreedySearch](../src/stelarx/pro/S6GreedySearch.java),
[S6DynamicTripletSupport](../src/stelarx/pro/S6DynamicTripletSupport.java),
[native colour scorer](../src/native/stelar_pro_s6_triplets.cpp),
[focused CUDA scorer](../src/native/stelar_pro_s6_triplets.cu).

### 2–5. Hashing, canonical clades, original partitions and DP transitions

These construction loops are mostly serial today. Per-tree occurrence prefix
arrays coexist with the duplicate-invariant subtree hashes. The latter merge
species sets small-to-large, retain hashes only for required speciation-rooted
bipartitions, and release donor bags. Candidate and gene partitions are
deduplicated by canonical hashes; repeated observed partitions retain frequency.
DP transitions are deduplicated before the weight phase. S4 uses the original
tree hashes for observations; S6 constructs separate gene and guide hashes.

After their consumers finish, `UniqueTaxonSubtreeHashes.release()` drops its
node lookup maps. Both prefix-array references are nulled before final weights,
and Main requests GC. Candidate exemplars and original membership indexes remain
needed by weights; DP and cluster tables remain needed by reconstruction.

Opportunities: parallel per-tree hashing with worker-owned occurrence scratch,
then stable table merges; primitive cluster/split IDs and CSR transitions instead
of many HashMap/HashSet/hash-array objects; fuse traversals where the same nodes
are scanned repeatedly. Test equal-score tie order as well as maximum score:
changing hash-table insertion order may change the selected topology even when
the candidate set and optimum score are unchanged.

Code: [PrefixHashArrays](../src/stelarx/hash/PrefixHashArrays.java),
[UniqueTaxonSubtreeHashes](../src/stelarx/pro/UniqueTaxonSubtreeHashes.java),
[ClusterTable](../src/stelarx/cluster/ClusterTable.java),
[PartitionTable](../src/stelarx/partition/PartitionTable.java),
[DPTable](../src/stelarx/dp/DPTable.java).

### 6. Final triplet weights: shared S4/S6 implementation

Reachability pruning selects unique splits whose parents can be reached from
the root. CPU mode scores independent splits on the fixed thread pool. CUDA's
active smaller-side kernel assigns **one thread per split**, and that thread
loops over every unique original rooted partition. Four child/cluster
intersections determine a binary partition's contribution, multiplied by its
frequency. It scans the smaller leaf range, suppresses repeated species copies,
and binary-searches sorted target position vectors. This is species-set
intersection, not copy-count intersection. Only original gene partitions supply
observations.

The GPU retains flattened orderings, taxon-position CSR, partitions and optional
range descriptors for one weight call. Split input/output buffers are batched;
static data is reused across batches and freed when JNI returns. Java receives
one score array, indexed by canonical split IDs. There is no full split ×
partition score matrix. The Java packing of the active smaller-side index is
currently serial. JNI may temporarily copy Java arrays, increasing peak RSS.

Let Q be reachable splits, P unique observed partitions, and I the average cost
of an intersection. Work is approximately O(Q P I); adding splits increases work
even when the final DP is tiny. In this kernel most parallelism is across Q;
the P loop inside each thread is serial. Therefore a small final Q can
underutilize a large GPU, while variable range sizes also produce uneven
thread work.

Opportunities:

1. Precompute each occurrence's previous-copy position. Testing `prev[pos] < lo`
   then identifies the first copy in a source range without another binary
   search. This adds one integer per occurrence but preserves exact semantics.
2. Use a sparse/bitset hybrid for species intersections. Bitsets work at higher
   taxon counts too: they need ceil(n/64) words, so the cost/memory advantage
   depends on n and range density. Retain multicopy distinct-species semantics
   and test polytomies/complements before enabling the dormant alternate path.
3. Tile the split × partition work and reduce exact integer partial scores.
   Warp/block cooperation can expose P parallelism for small Q and shorten long
   serial loops. Benchmark against extra reduction/traffic/register cost.
4. Cache reusable **cluster × child-set intersections in bounded tiles**, since
   many splits reuse a cluster. A full global intersection table is potentially
   enormous and should not be allocated.
5. Prefix-sum tree lengths, then pack tree/index segments in parallel; reuse
   native/device buffers between backbone calls where ownership is explicit.
6. Add immediate launch/copy error checks and optional batch-duration diagnostics
   around long native launches. The earlier R4 S2 stalls are not explained by
   this source audit; reduced batching alone was not a proven fix. Diagnostic
   coverage should precede throughput changes in that path.

A better GPU can help if this stage is compute or memory-throughput limited,
but launch/synchronization, host packing and data-dependent searches also matter.
Neither maximum Java threads nor more VRAM guarantees a proportional speedup.
Exact score parity is required for every proposed kernel replacement.

Code: [WeightTable](../src/stelarx/weight/WeightTable.java),
[CUDA smaller-side kernel and JNI](../src/native/stelarx_weight.cu).

### 7. DP and reconstruction

The inference DP is a serial memoized recursion over reachable clusters.
Scores and best splits are stored once per visited cluster; reconstruction uses
candidate exemplars. No CUDA DP is active in the default local modes.
Cross-tree candidate construction has a separate optional CPU/GPU path, which
does not run by default. Weight scores/DP tables live until reconstruction ends;
the JVM releases them after the analysis ends.

Opportunities: if profiling makes this phase significant, use dense cluster IDs,
size-ordered dependency levels and parallel per-level evaluation. Preserve the
current tie rule. On R4 this phase took milliseconds to 0.3 seconds, so it is a
low-priority target compared with S4 weights or S6 guide construction.

Code: [Inference](../src/stelarx/dp/Inference.java),
[Threading](../src/stelarx/util/Threading.java),
[PhaseLogger](../src/stelarx/PhaseLogger.java).

## Suggested next experiments

For memory, first replace dense indexes on uncompleted DISCO components and
shorten S6 group lifetimes. For S4 speed, first compare previous-copy indexing
and exact sparse/bitset intersections, then cooperative partition tiling.
For S6 speed, first batch fixed-backbone placements and profile its separate
restriction, backbone-weight, group-matrix and NNI work. Keep input, seed,
candidate hashes, optimum scores and topology in parity checks; measure RSS
and process VRAM separately, running one inference at a time.
