# C++17 port validation

Validated locally on 2026-11-27. Only the new `cpp/` tree was changed; production
Java and original native/JNI sources were not edited. Reference classes were
compiled freshly from the existing Java sources, not taken from an old jar.

Environment: Linux x86-64, Intel Core i7-14700HX, GCC 15.2.0, CUDA 13.1.115,
RTX 4060 Laptop GPU (compute 8.9). Make and CMake CPU/CUDA builds both passed.
CMake Release test targets explicitly retain assertions. The CUDA build uses
a narrow compatibility preinclude for CUDA 13.1/glibc 2.43 math declarations;
it does not disable GNU/POSIX features needed by C++ threading.

## Correctness evidence

| Area | Checks | Outcome |
|---|---|---|
| CPU primitives | Exact set interning, trailing-zero normalization, unordered split deduplication, binary/polytomy brute-force triplet enumeration, repeated multicopy-root observations and >64-bit scores | Pass |
| Parser and workers | Quoted labels, nested comments, event tags, malformed Newick rejection, nested jobs, 1/2/4-thread reconfiguration, worker exceptions and pool recovery | Pass |
| DISCO and completion | Global/strict-lineage/duplication scopes, nested cuts, original-input immutability, exact triplet matrix values, anchor ranks and equal-distance ties | Pass |
| Raw Java S4 comparison | Nested fixture plus randomized 8/17/70-species families at 1 and 4 workers; scope masks, component trees, similarities, anchor orders, completed trees and UPGMA | Pass |
| Main CLI | 23 accepted inference cases on CPU **and** CUDA, plus 3 score-only, 2 taxa-extraction and 5 rejected configurations | 33 cases, 0 failures |
| Independent flags | Cross-tree on/off, exhaustive on/off, both together, `--search-mode` alias, explicit off override, pruning, hash seed count and compatibility selectors | Pass |
| Utility/diagnostic CLI | Native S1 matrix/UPGMA output versus Java, tag-only/mapping output, inference mapping rejection, log/weight dumps, score-only CUDA, input/hard-link/derived-output collision protection, legacy consensus reference bridge | 11 groups, pass |
| s3 above cutoff | 101/137/200 species × 1/4 rounds × 1/4 workers against Java S6; exact score, rooted topology and every backbone/final scored-split count | 12 CPU cases, pass |
| Cutoff/scale | 100, 101 and 2,000 species with one guide round/four workers against Java; 100 dispatches to S4, 101 uses S6 | Pass |
| s3 CUDA | 101/200/2,000 species; four-round cross-tree run at 200 species; non-default signed seed `-17` at 137 species | Pass |
| Placement/NNI oracles | Insertion deltas versus whole-tree rescoring, incremental colour support versus direct triplets, narrow/wide colour counts, one-pass group matrices versus explicit restrictions | Pass |
| CUDA weights | CPU parity across disjoint colourings, polytomies, multicopy observations, frequency aggregation, automatic/manual/disabled batching and forced 128-bit | Pass |
| CUDA overflow | Exact score `4294924345327050000000000`, exceeding unsigned 64 bits | Pass |
| CUDA similarity/focus | Exact matrix values, word boundaries 0/63/64/127/128/129/200, unsorted selected species, absent species, cache reuse/replacement and memory-limited nonresident staging | Pass |
| CUDA cross-tree | Exact subset/residual verification, bounded output tiles and deterministic transition order versus CPU | Pass |
| Sanitizers | AddressSanitizer + UndefinedBehaviorSanitizer CPU core suite; S6 algorithm suite also checked during development | No reported errors; leak detection disabled in this environment |

The main end-to-end checks compare **rooted topology**, not merely unrooted RF;
also exact final/objective scores, candidate clade sets and reachable scored
split counts. s3 compares every backbone weight call separately from its final
weight call. A topology match alone would not detect a changed candidate pool.

Compact interned sets intentionally replace randomized Java correctness keys.
Java-compatible insertion order, UPGMA ties, nearest-anchor ties, RNG and event
tags were retained where they affect output. In particular, synthetic completion
joins are admitted by S4 candidate streaming but retain their untagged status
when the S6 backbone path applies Java's event-sensitive filtering.

## Full real-input replay

The complete, species-mapped estimated input for `sim200_dup1_loss1`, replicate
01, contains 1,000 families and 201 observed labels (including the simulator's
outgroup). All three native modes match fresh Java S1/S4/S6 in:

- rooted output topology and exact final triplet score;
- complete candidate clade set;
- scored splits in every backbone/final weight call.

| C++ / Java mode | Exact final score | Scored splits |
|---|---:|---|
| s1 / S1 | 198,994,798 | 10,155 |
| s2 / S4 | 433,588,088 | 57,513 |
| s3 / S6 | 433,606,852 | Backbones: 5,708 / 5,803 / 5,793 / 5,860; final: 211 |

| Mode | Java time, s | C++ time, s | Java peak RSS, MiB | C++ peak RSS, MiB |
|---|---:|---:|---:|---:|
| s1 / S1 | 9.75 | 4.13 | 643.9 | 136.6 |
| s2 / S4 | 16.63 | 5.23 | 690.4 | 238.5 |
| s3 / S6 | 10.86 | 5.93 | 747.1 | 168.1 |

Resource measurements are stored in [the real-input results](validation/real-r1/results.json)
with commands, logs, trees, cluster dumps and GNU time reports alongside them.
They include preprocessing/rooting/tagging and inference, but exclude compilation
and input preparation. The harness launches exactly one inference at a time.
Four CPU workers and strict GPU execution were used for both implementations;
the Java reference used `-Xmx8g`. This is a single-replicate replay, not a
multi-replicate performance claim or evidence that JNI removal alone caused
the improvement.

Do not use the original copy-ID-labeled `est_g_trees.trees` directly for this
comparison. As in the existing benchmark, the input must use species labels.
Otherwise both implementations interpret thousands of copy IDs as distinct
species and can correctly agree on a biologically unintended zero-score tree.

## Reproduction and saved results

### Paired benchmark, 2026-11-27

The follow-up benchmark completed C++ s2/s3 and fresh Java S4/S6 on estimated
`sim200_dup1_loss1` replicates 01–10, sequentially on this host with 28 workers
and strict GPU execution (Java `-Xmx8g`). All 20 paired checks match rooted
topology, exact final score, every backbone/final scored-split count, candidate
count and RF. Candidate identities and individual weights were not dumped in
these timed runs; the earlier correctness tests above check those separately.

| Mode | Mean RF | Mean wall time, s | Mean peak RSS, MiB | Mean sampled GPU peak, MiB |
|---|---:|---:|---:|---:|
| C++ s2 | 7.1212% | 5.10 | 230.0 | 110.0 |
| Java S4 | 7.1212% | 14.80 | 638.9 | 115.8 |
| C++ s3 | 7.1212% | 5.87 | 167.3 | 101.4 |
| Java S6 | 7.1212% | 11.13 | 751.4 | 122.4 |

C++ s2 and s3 both completed the full 2,001-label/1,000-family R1 input.
s2 obtained RF 300/3,996 (7.5075%) in 688.86 s wall with 1,371.9 MiB peak RSS,
486 MiB sampled GPU peak and 576,617 final scored splits. s3 obtained RF
294/3,996 (7.3574%) in 120.64 s wall with 1,191.0 MiB peak RSS, 210 MiB
sampled GPU peak and 2,031 final scored splits. s3's rooted tree, exact score
(770415543956), candidate count and all four backbone/final weight-call split
counts match the archived local duplication-scope subset-local Java S6 run.
s2's rooted tree and exact score (770400631020) match the archived Java S4 run.
The user stopped the attempted fresh large Java S4 rerun; the saved summary
keeps archived Java measurements separate from the fresh C++ results.

Timing includes preprocessing/rooting/tagging and inference, but excludes
compilation, input preparation and comparison checks. GPU peaks are
process-attributed samples, hence lower bounds, not exact allocation maxima.
See [the complete measurements](validation/paired-s2-s3-20261127/SUMMARY.md),
[per-run CSV](validation/paired-s2-s3-20261127/results.csv) and
[parity checks](validation/paired-s2-s3-20261127/parity.json). Dataset-local
copies are also saved under each dataset's `stelar-cpp-java-s2-s3-20261127/`.

### Reproduction commands

All commands below run from the repository root:

```bash
make -C cpp -j4 test
make -C cpp -j4 test-gpu CUDA_ARCH=sm_89

javac -sourcepath src -d cpp/validation/java-reference cpp/test/JavaS4Reference.java
python3 cpp/test/check_s4_reference.py cpp/build/test_s4 cpp/validation/java-reference

python3 cpp/test/parity.py --native cpp/bin/stelar-pro-cpp-cpu \
  --gpu cpp/bin/stelar-pro-cpp --output cpp/validation/cli
python3 cpp/test/check_cli_utilities.py --native cpp/bin/stelar-pro-cpp \
  --java-classes cpp/validation/cli/java-classes --work cpp/validation/utilities --gpu

python3 cpp/test/check_s6_reference.py --native cpp/bin/stelar-pro-cpp-cpu \
  --java-classes cpp/validation/cli/java-classes --work cpp/validation/s3-cpu
python3 cpp/test/check_s6_reference.py --native cpp/bin/stelar-pro-cpp \
  --java-classes cpp/validation/cli/java-classes --work cpp/validation/s3-boundaries \
  --gpu --taxa 100 101 2000 --rounds 1 --threads 4
python3 cpp/test/check_s6_reference.py --native cpp/bin/stelar-pro-cpp \
  --java-classes cpp/validation/cli/java-classes --work cpp/validation/s3-cross \
  --gpu --cross-tree --taxa 200 --rounds 4 --threads 4
python3 cpp/test/check_s6_reference.py --native cpp/bin/stelar-pro-cpp \
  --java-classes cpp/validation/cli/java-classes --work cpp/validation/s3-seeded \
  --gpu --taxa 137 --rounds 4 --threads 4 --seed -17

python3 cpp/test/check_dataset_reference.py \
  --input ../data/sim200_500/sim200_dup1_loss1/stelar-astral-first10-20261020/01/input_species_labels.tre \
  --native cpp/bin/stelar-pro-cpp --java-classes cpp/validation/cli/java-classes \
  --output cpp/validation/real-r1 --gpu
```

Main saved CLI results: [parity JSON](validation/cli/results.json),
[utility/diagnostic JSON](validation/utilities/results.json).
Additional saved s3 evidence: [CPU rounds/workers](validation/s3-cpu/parity.json),
[cutoff and 2,000 species](validation/s3-boundaries/parity.json),
[CUDA cross-tree](validation/s3-cross/parity.json),
[non-default seed](validation/s3-seeded/parity.json).
Build products and raw test evidence are ignored by Git; the source fixtures,
harnesses and this validation record are retained.

## Limits and explicit compatibility decisions

This validation is broad, not a proof over arbitrary inputs, hardware or all
possible combinations of flags. The ten-replicate paired performance study
above is complete; both native modes finished the full 2,000-species/1,000-family
R1 input, while a fresh large Java comparison was stopped at the user's request.
The separate 2,000-species synthetic test is an integration/scale check with
12 families.
Large candidate spaces still require substantial bitset storage and weight work.

See [README.md](README.md#options-and-compatibility) for the native exact-score
policy, diagnostic TSV formats, dormant legacy selectors and the explicitly
Java-backed consensus-only verifier. Normal native inference needs Python and
the existing ASTRAL-Pro rooting executable, but no Java runtime or JNI library.
