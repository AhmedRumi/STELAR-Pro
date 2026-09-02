# R1 STELAR-Pro and ASTRAL-Pro3 benchmark

Input: `../all_gt.tre` (1,000 gene trees, 200 species)

Truth: `../s_tree.trees`

Environment: 28 CPU threads, CUDA unavailable. GNU `time` wall time includes
each program's complete inference command. STELAR-Pro was built before timing;
the build is excluded. Its required resolve/root/tag preprocessing is included.

| Method | Raw RF | RF rate | Topological agreement (`1-RF`) | Wall time | Max RSS |
|---|---:|---:|---:|---:|---:|
| STELAR-Pro S1/I1, local DP | 28/394 | 0.0710659898 | 92.8934% | 128.93 s | 1,131,124 KiB |
| ASTRAL-Pro3 v1.25.3.8 | 24/394 | 0.0609137056 | 93.9086% | 53.56 s | 494,196 KiB |

Lower RF is better. On this replicate ASTRAL-Pro3 had an RF rate lower by
0.0101522842 and was 2.41 times faster by wall time.

## Commands

```bash
REPL=pro_data/simphy-200taxa-1000gt-10rep-dupmean-0.25-lb-0.8e-10-noloss/R1

/usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kb=%M\nexit_status=%x' \
  ./stelar-pro --no-build --cpu --threads 28 \
  --search-space S1 --intersection-method I1 --search-mode local -q \
  -i "$REPL/all_gt.tre" -o "$REPL/benchmark-results/out-stelar-pro.tre"

/usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kb=%M\nexit_status=%x' \
  ASTER-Linux/bin/astral-pro3 --thread 28 --seed 42 --verbose 1 \
  --output "$REPL/benchmark-results/out-astral-pro3.tre" "$REPL/all_gt.tre"

scripts/calculate_rf_rate.py "$REPL/benchmark-results/out-stelar-pro.tre" "$REPL/s_tree.trees"
scripts/calculate_rf_rate.py "$REPL/benchmark-results/out-astral-pro3.tre" "$REPL/s_tree.trees"
```

The RF values were also cross-checked using the independent edge-removal split
oracle in `test/test_rf_rate.py`.
