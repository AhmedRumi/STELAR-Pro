# S6 promotion validation, 2026-11-19

Production S6 now defaults to subset-local duplication scope. Saved checks:

- `replay/parity.json`: all ten estimated-gene `sim200_dup1_loss1` replicates
  match the saved S6-subset-duplication runs in topology (RF 0), exact final
  triplet score and each backbone/final weight call's scored-split count.
  The eleventh row is S4 strict replicate 01, which also matches.
- `replay/<rep>/<reference>/`: commands, logs and output trees from this replay.
  Commands record their original `/tmp` output paths; results were subsequently
  copied here for retention. This was a correctness replay, not a timing trial.
- `dynamic-parity.log`: existing exhaustive CPU/CUDA focused-triplet, colour
  query and proposal parity test passed.
- `old-s6.log` / `old-s6.tre`: frozen S6 without subset-local successfully ran
  on the complete 201-species replicate 01 input with one guide round. Its log
  confirms a shared global strict-lineage guide pool and one global UPGMA
  candidate, independently of the new production build.
- Scope/local-restriction unit tests passed: sister-locus additions,
  ancestral-sibling and absent-species exclusions, nested duplications,
  64-bit word boundaries, unchanged strict S4 scope and input immutability.

The focused regression suite passed the inference/default/compatibility checks
before stopping in an unrelated ASTRAL-Pro3 runner notification test.
`astral-runner-trace.log` locates the failed assertion: it expects RF 0.0000,
but the legacy wrapper's `rf.py` returns NA because this host lacks `dendropy`.
Running `python3 -B rf.py ...` independently reports that missing dependency.
`rf.py` and the ASTRAL runner were not changed in this task. The replay's RF
comparison uses the dependency-free `scripts/calculate_rf_rate.py` instead.

Replay command, from the repository root:

```sh
python3 -B test/check_s6_promotion_parity.py \
  --archive ../data/sim200_500/sim200_dup1_loss1/stelar-s4-s6-duplication-scope-20261118 \
  --work /tmp/stelar-s6-promotion-replay
```

Implementation/phase audit: [DOCS/s4-s6-parallelism-memory-audit.md](../../DOCS/s4-s6-parallelism-memory-audit.md).
