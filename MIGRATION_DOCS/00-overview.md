# ASTRAL-X to STELAR-X migration

This codebase was converted in place. The migration keeps the optimized data
plane and replaces only the semantics that necessarily differ between an
unrooted quartet method and a rooted triplet method.

## Retained unchanged in principle

- constant-size, multi-seed cluster hashes and prefix hash arrays;
- range and multi-range cluster exemplars;
- hash-binned DP tables and hash-subtraction cross-tree recombination;
- reachability pruning and memoized top-down inference;
- S1/S2/S3 presets, similarity matrices, UPGMA guidance, and consensus enrichment;
- incomplete-tree support, now with root-preserving completion;
- I1 smaller-side, I2 prefix-sum, I3 tree-walk, and I4 bitset scoring;
- exact `long`/`Int128` and optional large-score `double` accumulation;
- parallel parsing, CPU scoring, batching, diagnostics, and portable tooling.

The implementation identity is now STELAR-X end to end: Java sources live in
`src/stelarx`, the entry point is `stelarx.Main`, JNI exports use the
`Java_stelarx_*` namespace, and every native library is named `libstelarx_*`.
Launchers, test utilities, experiment output directories, environment variables,
portable JARs, and crash logs use the same name.

## Necessary semantic changes

1. Input roots are authoritative. No parser or completion path may reroot a
   supplied gene tree.
2. The cluster set contains descendant clades only. Complement orientations of
   unrooted edges are not inserted.
3. Local DP transitions are only `sub(u) -> sub(left(u)) | sub(right(u))`.
   Complement-side rotations and outgroup anchoring are removed.
4. Every rooted gene-tree internal node contributes a child bipartition,
   including the supplied root.
5. Split weights use rooted-triplet arithmetic, documented in
   `01-objective-and-data-flow.md`.
6. User-visible quartet labels are replaced by triplet labels.

For binary rooted trees these changes reduce cluster and transition counts:
one descendant orientation replaces two edge orientations, and complement
rotations disappear. The retained arrays and hash tables therefore have equal
or lower asymptotic memory requirements.
