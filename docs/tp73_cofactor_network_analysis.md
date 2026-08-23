# TP73 cofactor co-occurrence and outcome network plan

## Scientific question and terminology

The network analysis asks which motif contexts occur together around a TP73
anchor and which combinations predict TP73 CUT&RUN support or signal intensity
beyond their individual effects. It does not infer that the corresponding
proteins physically interact, nor that one transcription factor regulates
another.

CUT&RUN coverage is an occupancy-related assay signal, not an expression
level. Its continuous role is analogous to a quantitative response, but output
columns and reports should call it `cutandrun_signal` or `occupancy_signal`, not
expression. Gene expression becomes a separate downstream outcome after
anchors have been related to promoters and transcripts.

## Recommendation

Use each of the six exclusive motif interval-distance classes as a separate
primary analysis. Overlap, adjacency, and distal proximity encode different
physical hypotheses, and pooling them can cancel or blur effects. Add a joint
distance-interaction model as a secondary analysis so that differences between
bands are tested directly rather than inferred from separate significance
labels.

The primary network should have two explicitly different edge layers:

1. a **co-occurrence layer**, describing motif pairs found together more or
   less often than expected from their marginal prevalences; and
2. an **outcome-interaction layer**, describing motif pairs whose joint
   presence changes TP73 CUT&RUN support or intensity beyond additive motif
   effects.

This distinction prevents genomic sequence similarity or shared GC/repeat
preference from being interpreted as functional cooperation.

## Published methods and their fit

### TF-COMB: primary co-occurrence reference

[TF-COMB](https://pmc.ncbi.nlm.nih.gov/articles/PMC9358416/) was designed for
transcription-factor co-occurrence, distance, orientation, repeated sites, and
TF-pair networks. Its adapted market-basket framework reports support,
confidence, lift, and related association measures. This closely matches the
anchor-as-transaction representation already present in this project.

Use TF-COMB's definitions as an external methodological reference and parity
check. The project still needs its own outcome extension because TF-COMB does
not model a matched anti-p73-versus-negative-control CUT&RUN response or direct
TA-versus-DN differences.

### WGCNA: secondary module sensitivity only

[WGCNA](https://pmc.ncbi.nlm.nih.gov/articles/PMC2631488/) clusters variables
from a quantitative sample-by-feature matrix and summarizes modules with
eigengenes. A binary anchor-by-motif matrix can technically be supplied, but
its Pearson correlation becomes a prevalence-sensitive phi coefficient;
millions of spatially dependent anchors are not independent expression
samples, and the scale-free soft-threshold assumption has no privileged
biological meaning here.

WGCNA is therefore useful only as an exploratory comparison after aggregating
standardized motif prevalence or maximum-score summaries by 5 Mb genomic block,
series, isoform, and distance band. Such modules describe broad sequence-context
covariation, not local cofactor partnerships. They should be compared with, not
substituted for, the local co-occurrence network.

### Mixed graphical models: conditional-dependence sensitivity

The [`mgm` mixed graphical model](https://www.jstatsoft.org/article/view/v093i08/0)
supports categorical, count, and continuous nodes. It is a reasonable
conditional-dependence sensitivity analysis for binary motif presence plus
continuous CUT&RUN signal. It is not the first production method because more
than 1,000 partially redundant motif nodes make stability selection and motif
family consolidation essential.

### Gene-regulatory-network tools: later expression layer

PANDA starts with motif-derived TF-to-gene priors and integrates expression and
optional protein interactions by message passing
([original method](https://pmc.ncbi.nlm.nih.gov/articles/PMC3669401/),
[pandaR](https://pmc.ncbi.nlm.nih.gov/articles/PMC5870629/)). It becomes a
particularly natural option once the promoter/TSS bridge and gene-expression
matrix are joined.

[GENIE3](https://pmc.ncbi.nlm.nih.gov/articles/PMC2946910/),
[Inferelator](https://pmc.ncbi.nlm.nih.gov/articles/PMC9048651/), and
[SCENIC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5937676/) infer regulator-to-
target relationships from expression variation, optionally constrained by
motif or chromatin priors. They should not be applied directly to the present
motif-presence/CUT&RUN table and labelled as a gene-regulatory network. SCENIC+
would become relevant if matched single-cell chromatin accessibility and RNA
data are added.

## Analysis grain

The base observation is one physical TP73 anchor in one biological series and
isoform. Motif features are generated independently in each distance band:

- `strict_negative`: no retained locus at score `>= -1`;
- `intermediate`: a retained locus exists but is below the positive threshold;
- `positive`: the band maximum reaches the positive threshold;
- `maximum_score`, occurrence count, best-hit distance, orientation, and pair
  state remain available as continuous or descriptive features.

Run two threshold encodings:

1. motif-specific calibrated operating thresholds, the primary model; and
2. a common score-0 sensitivity analysis, which makes POU-family and other
   cross-motif comparisons easier to interpret.

Intermediate observations must not be silently recoded as strict negatives.
For binary association rules they are omitted; score-aware models retain their
continuous maximum score.

## Phase 1: co-occurrence edges

For every distance band, treat each TP73 anchor as a transaction containing its
positive motifs. For every sufficiently supported motif pair, calculate:

- pair support and expected support under marginal independence;
- confidence in both directions;
- lift, leverage, conviction, Jaccard similarity, and the pair odds ratio;
- counts stratified by chromosome, 5 Mb block, series, isoform, TP73 score band,
  and gene-relation class;
- orientation and same-motif pair summaries where applicable.

Use chromosome/block-preserving permutations or a 5 Mb block bootstrap for
uncertainty. Filter only on a prespecified minimum pair count and stability,
not on the same outcome later used for evaluation. Consolidate highly similar
JASPAR matrices into motif families while retaining the original matrix-level
edges as provenance.

## Phase 2: CUT&RUN outcome edges

Retain the existing binary strict-immersion/support response and add a
continuous matched-control outcome such as

```text
log2((anti_p73_window_signal + alpha) /
     (negative_control_window_signal + alpha))
```

with normalization fixed before model fitting. A two-part model is preferable:

1. logistic or conditional-logistic occupancy component for presence of
   technically valid TP73 support;
2. robust continuous/count component for signal intensity among informative
   anchors.

For a candidate pair `A`, `B`, fit within each distance band and isoform:

```text
outcome ~ TP73_score + A + B + A:B + series + genomic_covariates
```

The interaction coefficient `A:B`, not raw co-occurrence, defines the
outcome-interaction edge. Include gene relation, GC, mappability, repeat/ALU
annotation, chromosome class, and local background coverage as extensible
covariates. Do not adjust the primary total-effect model for variables that are
downstream mediators of motif context.

Estimate TA and DN separately, then test their interaction-coefficient
difference directly. A joint hierarchical model with `distance_band` and
`isoform` interactions is the secondary comparison and can borrow information
without erasing the six primary strata.

## Candidate-pair control and validation

Testing all roughly 500,000 matrix pairs in every band and outcome is possible
but scientifically wasteful because many JASPAR matrices are redundant. Use
this order:

1. collapse or tag motif-similarity families;
2. require minimum positive support in both CUT&RUN outcome classes;
3. retain every prespecified POU, PATZ1, SP1, E2F1, REST, TFAP2C, and TCF7 pair;
4. screen remaining pairs only by outcome-independent co-occurrence support;
5. fit additive and interaction models with stability selection;
6. control FDR separately by distance band and declared edge family.

Validation must use held-out chromosomes, with 5 Mb block resampling inside
training chromosomes. Fit Saos-2 and `skmel29_2` separately before pooling and
require directional consistency. Report effect sizes and interval estimates;
millions of anchors must not turn negligible effects into headline findings.

## Modules and graph summaries

Build communities from stable signed edges using a standard graph community
method such as Leiden. Keep co-occurrence and outcome-interaction networks
separate, then compare their modules. A WGCNA block-prevalence analysis and an
`mgm` conditional graph are sensitivity layers. Agreement across these methods
is stronger evidence than any single module assignment.

## Proposed query tables

Keep long tables and materialize only compact sufficient statistics:

- `cofactor_pair_cooccurrence`: threshold policy, distance band, motif pair,
  support/confidence/lift/odds ratio, block uncertainty, and q-value;
- `cofactor_pair_outcome_interaction`: outcome definition, isoform, series,
  motif pair, additive and interaction coefficients, uncertainty, q-value, and
  validation status;
- `cofactor_network_edge`: normalized edge view with an explicit
  `edge_type` of `cooccurrence` or `outcome_interaction`;
- `cofactor_module_membership`: method, network version, module, motif/family,
  and membership strength;
- `cofactor_module_outcome`: module score association with binary and continuous
  TP73 outcomes;
- `cofactor_network_run_config`: source manifests, score policy, response
  definition, covariates, split manifest, software versions, and random seeds.

Anchor-level feature Parquet remains the source of truth. The network tables
are rebuildable summaries and must never replace the underlying score,
distance, orientation, pair-state, and annotation fields.

## Execution order

1. Complete the all-POU calibrated-versus-score-0 sensitivity analysis.
2. Produce per-band motif prevalence and pair-support inventories without a
   CUT&RUN outcome.
3. Validate TF-COMB metric parity on a small chromosome fixture and a selected
   real chromosome.
4. Fit the prespecified-pair outcome models, including every POU-family motif.
5. Expand via outcome-independent support screening and stability selection.
6. Add WGCNA and `mgm` sensitivity analyses.
7. Join promoter/gene/expression data and evaluate PANDA or Inferelator as a
   genuinely gene-regulatory second network layer.
