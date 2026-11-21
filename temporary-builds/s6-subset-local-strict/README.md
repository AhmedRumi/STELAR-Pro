# S6 subset-local, strict lineage scope (temporary build)

This launcher selects the pre-promotion S6 subset-local path with strict
lineage-scoped DISCO components. It shares the patched archive jar under
[`../s6-without-subset-local/`](../s6-without-subset-local/README.md).
For 100 or fewer distinct species, the shared entry point runs the ordinary
S4 pipeline. Above 100, it runs the previous subset-local strict variant;
the current production S6 instead uses duplication scope.

```sh
bash temporary-builds/s6-subset-local-strict/run.sh \
  --input genes.trees --output strict-s6.tre --gpu-strict --threads 28 --xmx 9g
```

The original frozen archive jar remains available unchanged. This variant
uses the repository's current ASTRAL-Pro3 and native libraries.
