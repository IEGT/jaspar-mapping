# Human-source TP73 isoform cofactor highlights, autosomes

## Scope and provenance

This report summarizes the human-source sensitivity subset of the completed
JASPAR 2026 vertebrate-motif analysis around TP73 anchors on GRCh38 autosomes
1-22. It compares TAp73 and DNp73 CUT&RUN support in six mutually exclusive
motif interval-distance bands. `skmel29_1` is excluded; the valid series are
Saos-2 and `skmel29_2`.

The production manifest records:

- run ID `jaspar2026_grch38_tp73_distance_species_vertebrates_v1`;
- schema version 3;
- scientific finalizer commit `ee331ee718592233525bbb6954fe44e22a80f0b1`;
- completion at `2026-08-23T11:40:26Z`;
- 856 JASPAR matrices whose source-species metadata include *Homo sapiens*;
- 5,136 motif-band rows, of which 5,134 are estimable;
- Mantel-Haenszel odds ratios stratified by series and TP73 score band;
- paired delete-one-5-Mb-block jackknife uncertainty for TA-versus-DN
  contrasts.

The two non-estimable rows are ZNF354C (`MA0130.1`) at 51-100 bp and
101-150 bp because the strict-negative class is below 1%.

The transferred source files are kept outside Git under
`dry_runs/tp73_distance_species_isoform_highlights_autosomes_20260823/`:

- `cofactor_distance_isoform_comparison_human_source.tsv`: all 5,136 rows;
- `cofactor_distance_highlight_20_human_source.tsv`: 773 labelled presentation
  rows;
- `manifest.json`: run and source provenance.

Their SHA-256 sums were checked byte-for-byte against the files on Haumea:

```text
eeb699b568b55125089b66d36ee5dc28223410b285ef439fbcf7d2f42c7f81fc  cofactor_distance_highlight_20_human_source.tsv
97bfe213361e34e10769964a766c8ceb15d984121980a51ed35354ad8278433d  cofactor_distance_isoform_comparison_human_source.tsv
8e3dd90b209c42dc90bfbbcc61e55b4aeed97bdb5e0365a1af16ec24ab943ad9  manifest.json
```

Here, "human-source" means that human material contributed to construction of
the JASPAR matrix. It is a provenance sensitivity subset, not a claim that the
sequence preference is human-specific.

## How to read the tables

For each TP73 anchor, motif, and distance band, the analysis compares a
motif-positive context at its empirical operating threshold with a strict
negative context having no retained hit at score `>= -1`. Intermediate contexts
are excluded from that contrast. An odds ratio above one is enrichment among
anti-p73-supported anchors relative to the matched negative-control contrast;
an odds ratio below one is depletion.

The highlight table contains four independent data-driven selections per
distance band: strongest significant TA enrichment/depletion, strongest
significant DN enrichment/depletion, strongest direct TA-versus-DN differences
in either direction, and strict opposite-direction effects. It then appends the
prespecified POU2F2, SP1, PATZ1, REST, and E2F1 panel without changing any rank.
Association lists are ranked by absolute log odds ratio, not by p-value alone.

## Numerical overview

Significance is common because the analysis has many anchors. Effect size,
cross-band persistence, series consistency, and motif-family redundancy are
therefore more informative than q-values alone.

| Distance band | TA enriched/depleted at q<0.05 | DN enriched/depleted at q<0.05 | Shared TA/DN top-20 enriched/depleted | Strict reversals |
|---|---:|---:|---:|---:|
| overlap | 306 / 477 | 338 / 474 | 17 / 14 | 9 |
| adjacent 0-5 bp | 343 / 416 | 381 / 419 | 17 / 16 | 7 |
| gap 6-20 bp | 367 / 440 | 385 / 432 | 13 / 19 | 3 |
| gap 21-50 bp | 372 / 443 | 380 / 439 | 10 / 19 | 1 |
| gap 51-100 bp | 364 / 450 | 381 / 447 | 8 / 20 | 2 |
| gap 101-150 bp | 362 / 450 | 374 / 447 | 9 / 18 | 1 |

The depleted top-20 sets remain highly concordant between isoforms, especially
beyond 6 bp. Enriched top-20 concordance falls from 17 shared motifs at overlap
and adjacency to 8-10 at 21-150 bp, suggesting more isoform-specific distal
sequence contexts.

## Stable highlights

The following are the most stable entries when motifs are ranked first by the
number of distance bands in which they enter the respective top 20, then by
their mean within-band rank:

| Selection | Stable leading motifs (number of six bands in top 20) |
|---|---|
| TA enriched | CGGBP1 (6), ZNF142 (6), SLC2A4RG (6), ZBED4 (6), TFAP2E (6) |
| DN enriched | CGGBP1 (6), ZBED4 (6), ZNF213 (6), PLAG1 (6), ZNF395 (6) |
| TA depleted | SRY (6), FOXD1 (6), PBX1 (6), ZNF35 (6), GSX1 (6) |
| DN depleted | SRY (6), ZNF35 (6), FOXD1 (6), PBX1 (6), GSX1 (5) |
| Larger DN odds ratio | HELT (6), FEZF2 (6), ZNF213 (6), TCF12 (6), INSM1 (5) |
| Larger TA odds ratio | POU4F3 (6), POU3F1 (5), ZBTB40 (5), POU3F3 (4), FOXB1 (4) |

CGGBP1 is the clearest stable enriched signature: it is TA rank 1 in every
band, with TA odds ratios 2.858-4.488 and DN odds ratios 3.320-5.400. The
TA-larger group above must not be misread as TA-specific enrichment: these
leading POU/ZBTB/FOX motifs are depleted in both isoforms but more strongly
depleted for DNp73, so their numerical odds ratio is larger for TAp73.

Strict sign reversals are rare: 23 of 5,134 estimable motif-band rows, covering
20 motifs. All are DN-enriched/TA-depleted; none is TA-enriched/DN-depleted.
Their effects are modest (TA OR 0.910-0.986; DN OR 1.014-1.157). Only EOMES,
TBX2, and ZBTB5 recur in two bands, so these should be treated as focused
follow-up candidates rather than a dominant genome-wide pattern.

## Prespecified panel

All 30 panel motif-band rows have significant individual TA and DN associations,
significant direct isoform contrasts, and consistent directions in Saos-2 and
`skmel29_2`. None is a strict sign reversal.

| Motif | Direction in both isoforms | TA OR range | DN OR range | TA/DN OR ratio | Relevant TA rank range | Relevant DN rank range |
|---|---|---:|---:|---:|---:|---:|
| PATZ1 (`MA1961.2`) | enriched | 1.439-1.729 | 1.826-2.084 | 0.783-0.829 | 11-66 | 12-96 |
| SP1 (`MA0079.5`) | enriched | 1.384-1.599 | 1.674-1.860 | 0.821-0.861 | 24-78 | 37-130 |
| E2F1 (`MA0024.3`) | enriched | 1.248-1.319 | 1.326-1.553 | 0.846-0.941 | 91-162 | 110-204 |
| REST (`MA0138.3`) | enriched | 1.184-1.234 | 1.346-1.399 | 0.866-0.905 | 156-218 | 159-213 |
| POU2F2 (`MA0507.3`) | depleted | 0.720-0.831 | 0.501-0.681 | 1.221-1.438 | 70-269 | 58-240 |

The directly useful biological summary is therefore not merely that PATZ1,
SP1, E2F1, and REST are enriched and POU2F2 is depleted. Every one of these
associations is stronger in magnitude for DNp73: the four enriched motifs have
`OR_TA / OR_DN < 1`, while POU2F2 has `OR_TA / OR_DN > 1` because its depletion
is deeper for DNp73. PATZ1 is strongest immediately adjacent to TP73 (0-5 bp:
TA OR 1.729, DN OR 2.084); POU2F2 depletion is strongest at overlap (TA OR
0.720, DN OR 0.501).

POU4F1 is also depleted for both isoforms and much more strongly for DNp73
(TA OR 0.556-0.747, DN OR 0.388-0.567). It enters the top-20 direct isoform
difference list in four bands and reaches rank 6 at overlap. TCF7 behaves in
the same qualitative direction (TA OR 0.421-0.805, DN OR 0.312-0.650) but does
not reach an association top 20 because many other depleted motifs have still
larger effect sizes.

## Interpretation limits

These are associations of JASPAR motif signatures with the technical TP73
CUT&RUN support contrast, not proof of physical protein recruitment or causal
cofactor action. Closely related matrices can represent correlated sequence
preferences, and the two valid cell-line series are not biological replicates.
The complete side-by-side table should therefore remain the source for family
consolidation, cell-line sensitivity work, gene-relation stratification, and
subsequent multivariable models; the highlight table is a presentation layer,
not a replacement for it.
