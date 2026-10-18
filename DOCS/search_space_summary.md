# STELAR-Pro search spaces S1–S6

This describes the **current implementation**, not six different statistical
objectives. S1–S6 change which rooted species-tree clusters and compatible
splits are available to the optimizer. The final optimizer in every path is
STELAR-Pro's rooted, speciation-driven **triplet** dynamic program (DP), scored
against the original rooted/tagged multicopy gene families. Completed DISCO
trees, consensus trees, UPGMA trees, and greedy proposal trees enlarge the
search space; none is treated as an additional gene-family observation.

## 1. Shared input, notation, and DP

### Preprocessing and the meaning of a candidate

For inference, `Main` first uniquifies copy labels as needed, arbitrarily
resolves input gene-tree polytomies, and invokes the bundled ASTRAL-Pro3 `-T`
backend to root and tag the resulting trees. An internal `D` is a duplication;
an unmarked internal node is a speciation. Artificial nodes created while
resolving a polytomy receive **no special exemption**: their subsequent tag
determines whether they can contribute. `-T` tag-only mode roots/tags and
exits without the inference-time polytomy resolution. See
[`Main.java`](../src/stelarx/Main.java),
[`GeneTreePolytomyResolver.java`](../src/stelarx/pro/GeneTreePolytomyResolver.java),
and [`GeneTreeRooterTagger.java`](../src/stelarx/pro/GeneTreeRooterTagger.java).

Let `G` be the *original* rooted/tagged multicopy gene trees, and `X` the
candidate descendant-cluster set. For a gene speciation node `u`, let its
rooted child species sets be `M_0(u)` and `M_1(u)`. Each set contains a species
at most once, even if that species has several gene copies below the child.
Duplication nodes contribute **no scoring partition** and emit no S1 candidate
split. A duplication-rooted child can still supply its descendant species set
as one side of its *parent's* speciation split. Candidate trees made by S2–S6
can emit their own synthetic clusters and splits, but cannot enter `G`.

`ClusterTable` hashes and deduplicates species sets, including singleton DP
base cases and the all-species root. `UniqueTaxonSubtreeHashes` merges child
sets small-to-large, so repeated copies add a species only once to a subtree's
sum/XOR signature; temporary child sets are released after merging. Hashes
are retained for candidate-producing speciation subtrees, not for standalone
duplication-rooted candidate splits. `PartitionTable` extracts rooted child
partitions from **G only**. `DPTable` assembles allowed rooted transitions.
Its local transitions come from candidate-source trees; in `FULL` mode it also
combines compatible clusters across trees. It does not admit overlapping
children as a species-tree split. See
[`ClusterTable.java`](../src/stelarx/cluster/ClusterTable.java),
[`UniqueTaxonSubtreeHashes.java`](../src/stelarx/pro/UniqueTaxonSubtreeHashes.java),
[`PartitionTable.java`](../src/stelarx/partition/PartitionTable.java), and
[`DPTable.java`](../src/stelarx/dp/DPTable.java).

### The rooted-triplet score actually optimized

Write a candidate rooted species split as `A | B`, with `A ∩ B = ∅`. For one
binary speciation node `u` in an original gene tree, define

```text
a_i = |A ∩ M_i(u)|,    b_i = |B ∩ M_i(u)|,    i ∈ {0,1}.
```

The number of matching, speciation-driven rooted triplets at this pair of
species/gene nodes is

```text
w_u(A|B) = C(a_0,2)b_1 + a_0 C(b_1,2)
         + C(a_1,2)b_0 + a_1 C(b_0,2)

2w_u(A|B) = a_0 b_1(a_0+b_1-2) + a_1 b_0(a_1+b_0-2).
```

`C(x,2)=x(x−1)/2`. For example, the first term chooses a pair of species in
`A ∩ M_0` and an outgroup species in `B ∩ M_1`; the second chooses the pair in
`B ∩ M_1` and outgroup in `A ∩ M_0`. The other two terms exchange gene children.
Thus the *gene* triplet is rooted at speciation node `u`, and the *species*
triplet is displayed by split `A|B`. Sets, rather than copy counts, implement
the rule that identical species-triplet observations with the same gene-tree
root count once. Observations at different roots/families retain their
multiplicity. The split weight is `W(A|B)=Σ_{u∈Spec(G)} w_u(A|B)`; the code may
deduplicate identical rooted child partitions and multiply by their frequency.

For an unresolved speciation node with child sets `M_i`, put
`a_i=|A∩M_i|`, `b_i=|B∩M_i|`, `a=Σ_i a_i`, and `b=Σ_i b_i`. The corresponding
linear-in-child-count formula is

```text
2w_u(A|B) = Σ_i [a_i(a_i−1)(b−b_i) + b_i(b_i−1)(a−a_i)].
```

Only the node's actual children appear: species *outside* its subtree are not
an extra child. In normal inference, input polytomies have already been
resolved; this formula matters for score-only mode and retained polytomies.
See [`MIGRATION_DOCS/01-objective-and-data-flow.md`](../MIGRATION_DOCS/01-objective-and-data-flow.md)
and [`WeightTable.java`](../src/stelarx/weight/WeightTable.java).

The rooted DP is, schematically,

```text
F({x}) = 0
F(C)   = max over allowed (A|B) with A∪B=C, A∩B=∅
         [ F(A) + F(B) + W(A|B) ].
```

`X` and the transition table constrain the maximum. A true clade being in `X`
does **not** force selection: it must be reachable through compatible child
transitions and have a favorable total triplet score. This is why adding
oracle true-tree clades is a diagnostic of search-space coverage, not a
guarantee of improved RF accuracy. The DP yields a rooted species tree.

The weight phase computes `|A∩M_i|` and `|B∩M_i|` through each gene tree's
postorder leaf ranges. A per-tree/per-species sorted **vector of all copy
positions** replaces the old one-position map. It walks the smaller side,
skips later copies of a species within that range, and binary-searches the
other range's position vector for membership. The GPU representation is a
packed CSR of the same vectors. See
[`TaxonPositionIndex.java`](../src/stelarx/tree/TaxonPositionIndex.java) and
[`WeightTable.java`](../src/stelarx/weight/WeightTable.java).

## 2. Guide scores are not the final objective

Three score types occur in this code. Keeping them separate is essential:

1. **STELAR triplet split weight `W`** above: final DP evidence from original
   multicopy tagged gene trees, shared by S1–S6.
2. **Taxon-pair similarity**: a guide matrix built from *single-copy DISCO*
   trees in S2–S4. It guides completion, UPGMA, and S3/S4 consensus resolution;
   it is never added to `W`.
3. **Greedy-tree proposal objective**: ASTRAL-Pro's speciation-quartet score
   in S5, or STELAR's rooted-triplet tensor in S6. It selects proposal trees;
   final DP scoring is still `W`.

### S2/S3 quartet pair similarity

For a single-copy tree `t` with `k_t` leaves, and a pair `a,b` both present,
consider the `C(k_t−2,2)` choices of two other leaves `x,y`. Let `N^Q_t(a,b)`
count those induced quartets displaying `ab|xy` (the pair `a,b` on one side).
Then

```text
sim_Q(a,b) = [Σ_{t containing a,b} N^Q_t(a,b)]
             / [Σ_{t containing a,b} C(k_t−2,2)],
dist_Q(a,b) = 1 − sim_Q(a,b).
```

If the denominator is zero, the matrix's missing-pair handling applies rather
than inventing support. The CPU implementation accumulates this same-side
count by scattering contributions from internal-node components; the GPU
implementation uses an Euler-tour/LCA distance identity and pair tiles.
`dist_Q` is a *derived view of the same matrix*, not a separately estimated
topological distance matrix. An optional legacy distance-completion/verification
route exists, but it is not what the S2/S3 presets select. See
[`SimilarityMatrixBuilder.java`](../src/stelarx/completion/SimilarityMatrixBuilder.java)
and [`DOCS/similarity-matrix-design.md`](similarity-matrix-design.md).

### S4 rooted-triplet pair similarity

For a rooted single-copy tree `t` with `k_t` leaves, the third leaf `c`
displays `ab|c` precisely when it is outside the subtree rooted at
`LCA_t(a,b)`. Hence

```text
N^T_t(a,b) = k_t − |descendants(LCA_t(a,b))|,
sim_T(a,b) = [Σ_{t containing a,b} N^T_t(a,b)]
             / [Σ_{t containing a,b} (k_t−2)],
dist_T(a,b) = 1 − sim_T(a,b).
```

This is a *pair affinity derived from rooted triplets*, not the full
three-index STELAR objective. S4 constructs **one** such matrix and shares it
among tree completion, UPGMA, and consensus Step A/Step B. It never builds
the S2/S3 quartet matrix on the normal S4 route. However, the existing tree
completion still uses a **four-point algebraic decision rule** on that triplet
matrix. For three representatives `a,b,c` and missing taxon `x`, it compares,
for example,

```text
s_a = sim[x,a] + sim[b,c] − sim[x,b] − sim[a,c]
s_b = sim[x,b] + sim[a,c] − sim[x,a] − sim[b,c]
s_c = sim[x,c] + sim[a,b] − sim[x,b] − sim[a,c].
```

Therefore “S4 has no quartet **similarity-counting** kernel” is accurate;
“S4 contains no four-point procedure” is not. See
[`TreeCompleter.java`](../src/stelarx/completion/TreeCompleter.java).

### ASTRAL-Pro quartet score used for S5 proposals

An ordinary gene-quartet count would count every gene-copy quartet separately.
ASTRAL-Pro instead optimizes shared induced-quartet **equivalence classes**
rooted at speciation events, avoiding duplication-driven overcounting. At the
conceptual level, for a proposed species tree `T`, its greedy objective is

```text
Q_AP(T) = Σ_{gene families g} Σ_{[q] ∈ E_g}
          1{ T displays the species quartet topology of [q] },
```

where `E_g` denotes the ASTRAL-Pro speciation-driven quartet equivalence
classes in family `g`. This is **not** `sim_Q` from S2/S3 and **not** `W` from
the final rooted-triplet DP. The native A-Pro implementation evaluates local
colorings (the three incident species-tree sides, colored `x,y,z`) using
incrementally updated gene-tree counts. Its exact node recurrence is useful
for understanding what is computed, without enumerating quartets:

```text
T_x(w) = T_x(u)+T_x(v)
       + u_y v_z(u_y+v_z−2) + u_z v_y(u_z+v_y−2)

Q(w) = Σ_{c∈{x,y,z}} [u_c T_c(v) + v_c T_c(u)]
     + u_x(u_x−1)v_y v_z
     + u_y(u_y−1)v_x v_z
     + u_z(u_z−1)v_x v_y.
```

Here `u,v` are the gene node's children; `u_c,v_c` are their color counts;
`T_y,T_z` are cyclic renamings of `T_x`. The implementation accumulates the
change in `Q` at **speciation** nodes and returns half the accumulated score
for the colored tripartition; at a duplication node, the `T` counters are
propagated but no local `Q` is added. Recoloring a taxon updates only affected
ancestor paths. This is the `Tripartition` counter in
[`ASTER-Linux/src/multitree.hpp`](../ASTER-Linux/src/multitree.hpp), reused by
the [`S5 native bridge`](../src/native/stelar_pro_s5_greedy.cpp). Quartet scores
choose insertions/NNI moves and proposal trees only. The STELAR backbone and
final DP never replace `W` with `Q_AP`.

### S6 rooted-triplet proposal score

S6 precomputes an outgroup-major tensor from original tagged genes:

```text
H(a,b;c) = Σ_{u∈Spec(G)} 1{a,b occur in one child species set of u,
                         c occurs in another child species set of u},
                         a≠b, a≠c, b≠c.
```

The pair `a,b` is unordered; a species is represented once per child set. For
a candidate rooted species split, S6 can evaluate

```text
W_H(A|B) = Σ_{a<b∈A, c∈B} H(a,b;c)
         + Σ_{a<b∈B, c∈A} H(a,b;c).
```

This equals the binary speciation-partition weight `W(A|B)` above, but offers
constant-time individual triplet lookups during greedy insertion and rooted
NNI. It requires `O(n³)` tensor cells (eight bytes each in the current
`long[]` representation), so S6 has a materially different memory profile
from S1–S5. It is not a second final objective. See
[`S6TripletSupport.java`](../src/stelarx/pro/S6TripletSupport.java).

## 3. What each search-space path builds

| Path | Candidate sources beyond original speciation-rooted S1 clusters | Guide/proposal score | DP mode |
| --- | --- | --- | --- |
| **S1** | None | None | Local |
| **S2** | Completed single-copy DISCO trees; one UPGMA tree | Quartet pair similarity | Full |
| **S3** | S2 sources plus greedy-consensus clades/refinements | Quartet pair similarity | Full |
| **S4** | S3-style sources/refinements | Rooted-triplet pair similarity | Full |
| **S5** | ASTRAL-Pro quartet-greedy proposal trees; accumulated adaptive proposals | A-Pro quartet-equivalence objective | Full |
| **S6** | Rooted-triplet greedy proposal trees; accumulated adaptive proposals | STELAR triplet tensor `H` | Full |

The presets are configured in [`CliPresets.java`](../src/stelarx/CliPresets.java).
S2→S3→S4 are related DISCO/consensus variants. S5 and S6 are **independent
paths**, not additions on top of S2–S4. Every path retains S1's input-derived
candidates. These branches share common DP/weight code; that does not mean
choosing S5 executes S2–S4 or incurs their matrix-building cost.

### S1: native speciation-rooted candidates

Traverse each original rooted/tagged gene tree. Register distinct species
sets below leaves and speciation-rooted internal nodes, and their valid local
rooted child transitions. Skip standalone duplication-rooted candidates, but
continue traversing through duplications: their descendant taxa are needed
for ancestor speciation nodes. Deduplicate equivalent sets across occurrences.
Construct local DP transitions, compute `W` against the original gene nodes,
then run the rooted DP. S1 is the default and has no DISCO, guide matrix,
UPGMA, consensus, or greedy-proposal stage. It is the narrowest search space;
an otherwise well-scored true clade/transition may simply be unavailable.

### S2: DISCO, completion, UPGMA, and cross-tree recombination

Start with **all S1 candidates**. Decompose each rooted/tagged multicopy gene
family with DISCO by cutting at duplications, retaining single-copy components
with at least four leaves. Build `sim_Q` from these components. Complete
incomplete components by inserting missing taxa under the matrix-guided
four-point procedure; completion preserves their supplied root and existing
relationships. Add the resulting complete single-copy trees as candidate
sources. Build one UPGMA guide from the *same* pair matrix (`dist=1−sim`) and
add its clusters too. In `FULL` mode, build both local and compatible
cross-tree DP transitions. Crucially, `PartitionTable` and `WeightTable`
still use the original rooted/tagged multicopy trees only. See
[`DiscoDecomposer.java`](../src/stelarx/pro/DiscoDecomposer.java),
[`TreeCompleter.java`](../src/stelarx/completion/TreeCompleter.java), and
[`UPGMAClusterer.java`](../src/stelarx/completion/UPGMAClusterer.java).

Normal S2 inference streams completion in bounded CPU-parallel batches. Each
batch's clades and local splits are merged into the global candidate/DP tables;
completed node graphs, raw DISCO components, and temporary hashes are then
released. Only leaf-order and position indexes for trees that introduced a new
candidate survive as weight/reconstruction exemplars. S3/S4 still materialize
completed trees because their consensus Step B reads source-tree topology.

UPGMA repeatedly merges the active clusters of greatest pair similarity;
the intercluster similarity is the size-weighted average of underlying pair
similarities. Its branch/merge order proposes a hierarchy; that hierarchy is
*not* interpreted as additional gene-tree evidence.

### S3: S2 plus greedy consensus

S3 performs S2's decomposition, quartet-matrix completion, UPGMA, and full
cross-tree recombination. It additionally counts clades/bipartitions in the
**completed DISCO trees**, sorts them by frequency, and greedily inserts
compatible clades into a laminar consensus forest. A single frequency-sorted
pass snapshots seven thresholded consensus trees at `0, .01, .02, .05, .10,
.20, 1/3`. Unresolved consensus polytomies produce extra candidates: Step A
forms a group similarity by averaging the existing taxon-pair matrix between
two polytomy arms, then runs small UPGMA; Step B samples/restricts source trees
and resolves additional relationships. Direct consensus clades and these
refinements are bridged into `X` (potentially as multi-range exemplars).
They receive cross-tree transitions in `FULL` mode, but consensus trees are
not appended to the gene evidence or mined as ordinary local-transition
trees. See [`GreedyConsensus.java`](../src/stelarx/greedy/GreedyConsensus.java),
[`PolytomyResolver.java`](../src/stelarx/greedy/PolytomyResolver.java), and
[`EmissionBridge.java`](../src/stelarx/greedy/EmissionBridge.java).

### S4: S3 construction with one triplet guide matrix

S4 keeps the same DISCO components, completion procedure, UPGMA guide,
consensus thresholds/refinements, and `FULL` DP mechanism as S3. Its deliberate
change is to replace `sim_Q` with `sim_T` **once** at matrix construction; the
triplet matrix is reused by completion, UPGMA, and consensus. Original tagged
genes still supply the final `W`. S4 therefore removes quartet *counting* from
this guide path, but inherits the four-point insertion logic as explained
above. It is not automatically faster than S3: both paths still do tree
completion, UPGMA, consensus, candidate handling, and the same final DP;
matrix construction is only one portion of runtime.

### S5: ASTRAL-Pro greedy proposals, STELAR triplet DP

S5 bypasses DISCO, pair-similarity completion, UPGMA, and consensus. It invokes
the bundled ASTRAL-Pro native placement/NNI search *in process* to generate
full-taxon greedy proposal trees using `Q_AP`, with randomized order/tie
behavior. For `n≥100`, each two-step round selects roughly
`floor(sqrt(n)·log₂(n)/4)` backbone taxa (27 when `n=200`). It first generates
quartet-greedy **subset guides**. A Java bridge then restricts original tagged
genes to that subset and infers the backbone by **STELAR triplet DP** with the
guides as candidates; it does not use A-Pro quartet DP for that backbone.
The A-Pro placement code partitions the remaining taxa by backbone branch,
places compatible groups, handles deferred taxa, and improves by NNI. The
resulting unrooted proposal edges contribute **both orientations** to `X`,
since their two sides are different possible rooted descendant clusters.
This does not mean choosing an arbitrary single root and discarding one side.
Final transition recombination and scoring are rooted STELAR DP.

Initially S5 requests four greedy proposals by default. After the first
triplet-DP inference it may make up to **five additional DP passes**. Each pass
adds another four proposals, seeding their subset backbones from the previous
STELAR-DP tree restricted to newly selected taxa, retains all earlier proposals
and S1 candidates, and reruns STELAR DP. It stops early when the *rooted
topology* is unchanged; it need not use all five passes. `--s5-rounds` changes
the proposals per pass. See
[`AProGreedySearch.java`](../src/stelarx/pro/AProGreedySearch.java),
[`S5OverallSearch.java`](../src/stelarx/pro/S5OverallSearch.java),
[`stelar_pro_s5_greedy.cpp`](../src/native/stelar_pro_s5_greedy.cpp), and
[`ASTER-Linux/src/algorithms.hpp`](../ASTER-Linux/src/algorithms.hpp).

### S6: rooted-triplet greedy proposals, STELAR triplet DP

S6 also bypasses S2–S4. Build `H(a,b;c)` once from original speciation-rooted
gene partitions. For each randomized proposal round, grow a **rooted** tree:
start with a two-taxon rooted tree, score insertion edges for each later
taxon using triplet-score changes, allow an insertion *above the current
root* when it wins, then apply rooted NNI refinement. Thus the root is not
fixed from the initial pair, but the algorithm does **not** exhaustively
reroot and rescore every completed tree at every insertion. Proposal trees
contribute rooted descendant clusters and transitions to `X`.

For `n≥100`, S6 uses the same `floor(sqrt(n)·log₂(n)/4)` subset size. Initial
rounds generate triplet-greedy subset guides and infer a rooted backbone by
STELAR DP on restricted original genes; remaining taxa are placed in
backbone-related groups and refined. Later adaptive passes seed their subset
backbone from the previous rooted STELAR-DP tree. Like S5, S6 initially uses
four proposals by default, may run five additional accumulate-and-rerun DP
passes, and stops when the rooted topology does not change. Its proposal
objective and root handling are different from S5; `H` is the only greedy
score. See [`S6GreedySearch.java`](../src/stelarx/pro/S6GreedySearch.java) and
[`S6OverallSearch.java`](../src/stelarx/pro/S6OverallSearch.java).

## 4. CPU and CUDA parallelization: precise scope

`--threads N` sets the shared CPU worker count. The compute mode is `AUTO`
(probe CUDA and use CPU if unavailable), `CPU`, or `GPU`; strict GPU mode
fails instead of silently falling back. A GPU run is **hybrid**: rooting,
tagging, DISCO, search orchestration, and rooted DP reconstruction are not
all CUDA kernels. Phase boundaries and actual backend choices are logged.

| Stage | CPU work | CUDA work when selected/available |
| --- | --- | --- |
| Preprocess; duplicate-invariant subtree hashes; partition extraction | Host-side tree traversals and small-to-large merges | No corresponding end-to-end CUDA tree traversal |
| S2/S3 quartet pair matrix | Reference scatter over source trees (the current `buildCPU` tree loop is sequential) | Euler/RMQ preparation on CPU workers; `(a,b)` pair-tiled similarity kernel, tree-data batching |
| S4 triplet pair matrix | Per-tree tour/RMQ builds in parallel, then disjoint pair-row workers | Parallel host tour/flattening; pair-tiled triplet-LCA kernels, compact/wide RMQ and batching |
| DISCO completion | Different source trees completed in parallel; S2 bounds the active batch and discards completed node graphs after candidate extraction; insertion within one tree remains ordered | No separate completion CUDA kernel |
| UPGMA | Parallel best-pair scans/reductions and row updates within a merge; merge sequence serial | No separate UPGMA CUDA kernel |
| S3/S4 consensus | Frequency-sorted laminar pass; distinct polytomy-resolution tasks scheduled across workers | No separate consensus CUDA kernel; it may consume a GPU-built matrix |
| Cross-tree DP-space construction | Parallel per-parent-cluster matching | CUDA complement/match builder; Java owns resulting transition table |
| Common STELAR `W` weight calculation | Parallel per-candidate-split workers; duplicate-aware position-vector membership | Batched one-thread-per-split smaller-side kernel, CSR copy positions, binary search and reduction; multi-range consensus clusters supported |
| Final rooted DP/reconstruction | Java memoized optimization and traceback | No CUDA kernel for the DP recurrence itself |
| S5 greedy proposals | A-Pro native worker partitions of gene families for quartet-color updates; search orchestration/placement/NNI | GPU color counters own gene families and reduce scores; host still controls placement/NNI and STELAR backbones |
| S6 tensor and greedy proposals | Native tensor builder partitions output by outgroup across threads (Java fallback exists); Java insertion/NNI | One CUDA block per outgroup; threads stream speciation partitions and atomically accumulate `H`; insertion/NNI stay Java |

The common GPU split-weight path uses the fixed
`SMALLER_SIDE_TRAVERSAL` intersection method in current `Config`, with
free-VRAM-aware candidate batching. Both CPU/GPU count duplicate taxa once
through sorted copy-position vectors. An alternative prefix-sum kernel remains
in the codebase, but it is **not the current default**. S5 GPU means native
quartet-counter acceleration, not that the entire greedy search runs on the
GPU. S6 GPU similarly accelerates tensor construction and common DP-space/
weight stages, not each sequential insertion/NNI decision. Kernel and host
entry points are in
[`stelarx_similarity.cu`](../src/native/stelarx_similarity.cu),
[`stelarx_dp.cu`](../src/native/stelarx_dp.cu),
[`stelarx_weight.cu`](../src/native/stelarx_weight.cu),
[`stelar_pro_s5_score.cu`](../src/native/stelar_pro_s5_score.cu),
[`stelar_pro_s6_triplets.cpp`](../src/native/stelar_pro_s6_triplets.cpp), and
[`stelar_pro_s6_triplets.cu`](../src/native/stelar_pro_s6_triplets.cu).

## 5. Interpretation and limits

- A larger candidate set may improve *best achievable* triplet score without
  improving RF against the true species tree. Search-space completeness and
  objective/model accuracy are separate questions.
- S2–S4 guide similarity is computed on **DISCO single-copy components**;
  original multicopy genes alone supply final weights. S5's A-Pro quartet
  score also stops at proposal generation. S6 is fully triplet-scored in its
  proposal search, but its tensor costs `O(n³)` memory.
- S5's bidirectional greedy edges are intentionally different from S6's
  rooted-cluster proposal policy. Neither changes the rooted triplet DP rule.
- GPU timing includes CPU preprocessing, host/device transfers, and Java
  search stages. “GPU mode” should not be interpreted as every phase being
  parallel on the device.

The implementation entry point and preset parser are
[`Main.java`](../src/stelarx/Main.java) and
[`CliPresets.java`](../src/stelarx/CliPresets.java). The mathematical rooted
objective is also summarized in the [migration documentation](../MIGRATION_DOCS/01-objective-and-data-flow.md);
the S5 greedy design follows the ASTRAL-Pro2 supplementary search-space
procedure with the explicit STELAR backbone/final-DP substitutions above.
