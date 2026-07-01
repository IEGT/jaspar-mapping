# JASPAR 2026 Motif Resource

The JASPAR 2026 CORE non-redundant motif matrix file is intentionally not
tracked in Git. For the updated run, retrieve it with:

```sh
make jaspar
```

Resource details verified on 2026-07-01:

- Source URL: https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt
- Local filename: `JASPAR2026_CORE_non-redundant_pfms_jaspar.txt`
- Motif records: `2633`
- SHA-256: `0dc9b7f9e159a8376c2e52edf373863ae21518bb760ceeea494ad6b746092e53`

The publication-era source snapshot used JASPAR 2022. Use the `v0.0.0` tag
for that exact code state, or set `JASPAR_VERSION=2022` when preparing the old
matrix with the updated Makefile.
