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

## Matrix taxonomy and source species

The official JASPAR 2026 CORE metadata table is tracked as
[`resources/jaspar/JASPAR2026_CORE_metadata.tsv`](../resources/jaspar/JASPAR2026_CORE_metadata.tsv)
because it is small and is required to interpret an all-taxon scan. It is
upstream data, retained byte-for-byte rather than hand-edited.

- Source URL: https://mencius.uio.no/JASPAR/JASPAR_metadata/2026/ultimate_metadata_table_CORE.tsv
- Matrix-version rows: `4572`
- SHA-256: `1b5e1d3818ac3f58796a7936d9d7614c0c0961ea55a2e025fb512a0b42ef7b6e`

The upstream `tax_id` and `species` fields can contain parallel `::`-separated
lists, and JASPAR can leave source species unspecified. Build the normalized
catalog with:

```sh
scripts/build_jaspar_metadata_catalog.py \
  --metadata resources/jaspar/JASPAR2026_CORE_metadata.tsv \
  --pssm JASPAR2026_CORE_non-redundant_pfms_jaspar.txt \
  --output "$HOME/resources/jaspar/2026/core_metadata" \
  --expected-metadata-sha256 \
    1b5e1d3818ac3f58796a7936d9d7614c0c0961ea55a2e025fb512a0b42ef7b6e \
  --expected-motif-count 2633
```

The catalog contains one `jaspar_matrix` row per matrix version, one
`jaspar_matrix_species` row per source species, and an exact
`jaspar_motif_set_matrix` bridge for the scanned non-redundant PSSM file.
`tax_group = 'vertebrates'` is the primary filter for interpreting putative
human cofactors. `includes_homo_sapiens` is a stricter provenance sensitivity,
not a claim that non-human vertebrate matrices are biologically inapplicable.

The publication-era source snapshot used JASPAR 2022. Use the `v0.0.0` tag
for that exact code state, or set `JASPAR_VERSION=2022` when preparing the old
matrix with the updated Makefile.
