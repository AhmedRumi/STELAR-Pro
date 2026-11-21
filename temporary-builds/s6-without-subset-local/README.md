# S6 without subset-local (temporary build)

Frozen on 2026-11-19 immediately before promoting subset-local,
duplication-scope S6. This jar contains the previous compiled production
classes. The changed source files are copied under `src/` for review;
the production build compiles only the repository's top-level `src/`.

Original jar SHA-256: `3de19a7b3027b34144e1a5e7b44cc33ab24c8265bb38e7d32fde0421f9d63fba`.
The launcher now uses `stelar-pro-s6-small-s4.jar`, made by compiling this
archive's patched `Main.java` against the original jar with
`rebuild-small-s4.sh`. The original jar remains unchanged. At 100 or fewer
distinct species, the patched entry point runs the ordinary S4 pipeline;
above 100, it retains this archive's original global S6 behavior.
Patched jar SHA-256: `e53b9f81b095630d6457bbc80755c72c734fa6537ffc714409c50ac8be7c3fc4`.

Run from any directory:

```sh
bash temporary-builds/s6-without-subset-local/run.sh \
  --input genes.trees --output old-s6.tre --gpu-strict --threads 28 --xmx 9g
```

This variant uses strict lineage scope, a shared global triplet similarity
matrix, globally prepared lineage-completed DISCO guides, and a global UPGMA
candidate in the final DP. It preserves four rounds and the previous seed
defaults. It uses the repository's existing native libraries and ASTRAL-Pro3
rooting/tagging executable. `--subsetlocal` is rejected by this wrapper so
the archived name has one meaning. The normal `./stelar-pro --search-space S6`
uses the newly promoted subset-local duplication-scope implementation.

The jar freezes the Java implementation, not the operating system, JVM, or
native libraries. Those shared dependencies can affect timing across machines.
