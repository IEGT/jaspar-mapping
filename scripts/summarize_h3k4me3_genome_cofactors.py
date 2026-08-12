#!/usr/bin/env python3

"""Join genome-wide TP73 enrichment with H3K4me3 cofactor effects."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


class SummaryError(RuntimeError):
    pass


OUTPUT_TABLES = (
    "h3_autosome_baseline",
    "motif_coverage",
    "joint_primary_motif",
    "joint_reference_zero_motif",
    "context_primary_effect",
    "gene_relation_primary_effect",
    "score_gradient_primary_effect",
    "tp73_interaction_primary_effect",
    "summary_metrics",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def run_sql(duckdb: Path, sql: str) -> None:
    process = subprocess.run(
        [str(duckdb), "-light-mode", "-bail", ":memory:"],
        input=sql, text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise SummaryError(process.stderr.strip() or "DuckDB extraction failed")


def query_scalar(duckdb: Path, parquet: Path) -> int:
    process = subprocess.run(
        [str(duckdb), "-light-mode", "-csv", "-noheader", ":memory:", "-c",
         f"SELECT count(*) FROM read_parquet({sql_string(parquet)});"],
        text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise SummaryError(process.stderr.strip() or "DuckDB validation failed")
    try:
        return int(process.stdout.strip())
    except ValueError as error:
        raise SummaryError("DuckDB returned an invalid row count") from error


def git_identity(source: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False,
    )
    dirty = subprocess.run(
        ["git", "-C", str(source), "status", "--short", "--untracked-files=no"],
        text=True, capture_output=True, check=False,
    )
    return (
        commit.stdout.strip() if commit.returncode == 0 else "unavailable",
        dirty.returncode != 0 or bool(dirty.stdout.strip()),
    )


def output_sql(inputs: dict[str, Path], staging: Path) -> str:
    enrichment = sql_string(inputs["enrichment"])
    intensity = sql_string(inputs["intensity"])
    context = sql_string(inputs["context"])
    gene_relation = sql_string(inputs["gene_relation"])
    gradient = sql_string(inputs["score_gradient"])
    interaction = sql_string(inputs["interaction"])
    h3_summary = sql_string(inputs["h3_summary"])

    copies = []
    for table in OUTPUT_TABLES:
        parquet = sql_string(staging / f"{table}.parquet")
        tsv = sql_string(staging / f"{table}.tsv")
        copies.append(
            f"COPY {table} TO {parquet} "
            "(FORMAT PARQUET, COMPRESSION ZSTD);\n"
            f"COPY {table} TO {tsv} "
            "(FORMAT CSV, HEADER true, DELIMITER '\\t', NULL 'NA');"
        )

    return f"""
SET preserve_insertion_order=false;
CREATE VIEW enrichment AS SELECT * FROM read_parquet({enrichment});
CREATE VIEW intensity AS SELECT * FROM read_parquet({intensity});
CREATE VIEW context_effect AS SELECT * FROM read_parquet({context});
CREATE VIEW gene_relation_effect AS SELECT * FROM read_parquet({gene_relation});
CREATE VIEW score_gradient AS SELECT * FROM read_parquet({gradient});
CREATE VIEW tp73_interaction AS SELECT * FROM read_parquet({interaction});
CREATE VIEW h3_summary AS SELECT * FROM read_parquet({h3_summary});

CREATE TABLE h3_autosome_baseline AS
SELECT series_id, min(cell_line) AS cell_line,
       sum(anchors)::BIGINT AS anchors,
       sum(anchors_with_h3k4me3_signal)::BIGINT AS anchors_with_h3k4me3_signal,
       sum(uninformative_all_h3k4me3_zero)::BIGINT
         AS uninformative_all_h3k4me3_zero,
       (sum(anchors*mean_delta_ta_vs_gfp)/sum(anchors))::DOUBLE
         AS mean_delta_ta_vs_gfp,
       (sum(anchors*mean_delta_dn_vs_gfp)/sum(anchors))::DOUBLE
         AS mean_delta_dn_vs_gfp
FROM h3_summary
WHERE try_cast(chrom AS INTEGER) BETWEEN 1 AND 22
GROUP BY series_id
ORDER BY series_id;

CREATE TABLE motif_coverage AS
WITH h AS (
  SELECT motif_id,
    count(*) FILTER (WHERE negative_reference_threshold = -1
                     AND evaluation_status = 'ok') AS strict_ok_rows,
    count(*) FILTER (WHERE negative_reference_threshold = 0
                     AND evaluation_status = 'ok') AS zero_ok_rows,
    min(source_score_floor) AS source_score_floor
  FROM intensity GROUP BY motif_id
)
SELECT e.motif_id, e.factor_name, e.positive_threshold,
       e.source_score_floor, e.evaluation_status AS enrichment_status,
       e.negative_reference_observable AS enrichment_negative_observable,
       e.adjusted_odds_ratio, e.confidence_interval_95_lower,
       e.confidence_interval_95_upper, e.q_value_bh_all_jaspar,
       e.association_direction,
       coalesce(h.strict_ok_rows, 0) AS h3_strict_ok_rows,
       coalesce(h.zero_ok_rows, 0) AS h3_zero_ok_rows
FROM enrichment e LEFT JOIN h USING (motif_id)
ORDER BY e.motif_id;

CREATE TEMP TABLE binary_primary AS
SELECT motif_id, min(factor_name) AS factor_name,
  count(*) AS strict_ok_rows,
  count(DISTINCT (series_id,isoform)) AS strict_designs,
  max(estimate) FILTER (WHERE series_id='saos2' AND isoform='TA') AS ta_saos,
  max(confidence_interval_95_lower) FILTER
    (WHERE series_id='saos2' AND isoform='TA') AS ta_saos_ci_lower,
  max(confidence_interval_95_upper) FILTER
    (WHERE series_id='saos2' AND isoform='TA') AS ta_saos_ci_upper,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='saos2' AND isoform='TA') AS ta_saos_q,
  max(anchors_positive) FILTER
    (WHERE series_id='saos2' AND isoform='TA') AS ta_saos_positive_anchors,
  max(anchors_negative) FILTER
    (WHERE series_id='saos2' AND isoform='TA') AS ta_saos_negative_anchors,
  max(estimate) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_skmel,
  max(confidence_interval_95_lower) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_skmel_ci_lower,
  max(confidence_interval_95_upper) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_skmel_ci_upper,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_skmel_q,
  max(anchors_positive) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_skmel_positive_anchors,
  max(anchors_negative) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_skmel_negative_anchors,
  max(estimate) FILTER (WHERE series_id='saos2' AND isoform='DN') AS dn_saos,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='saos2' AND isoform='DN') AS dn_saos_q,
  max(estimate) FILTER
    (WHERE series_id='skmel29_2' AND isoform='DN') AS dn_skmel,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='skmel29_2' AND isoform='DN') AS dn_skmel_q
FROM intensity
WHERE negative_reference_threshold=-1 AND evaluation_status='ok'
GROUP BY motif_id;

CREATE TEMP TABLE gradient_primary AS
SELECT motif_id,
  max(estimate) FILTER (WHERE series_id='saos2' AND isoform='TA')
    AS ta_score_gradient_saos,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='saos2' AND isoform='TA') AS ta_score_gradient_saos_q,
  max(estimate) FILTER (WHERE series_id='skmel29_2' AND isoform='TA')
    AS ta_score_gradient_skmel,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_score_gradient_skmel_q
FROM score_gradient
WHERE score_clamp_reference=-1 AND evaluation_status='ok'
GROUP BY motif_id;

CREATE TEMP TABLE intergenic_primary AS
SELECT motif_id,
  max(estimate) FILTER (WHERE series_id='saos2' AND isoform='TA')
    AS ta_intergenic_saos,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='saos2' AND isoform='TA') AS ta_intergenic_saos_q,
  max(estimate) FILTER (WHERE series_id='skmel29_2' AND isoform='TA')
    AS ta_intergenic_skmel,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_intergenic_skmel_q
FROM gene_relation_effect
WHERE negative_reference_threshold=-1 AND evaluation_status='ok'
  AND gene_relation_class='intergenic'
GROUP BY motif_id;

CREATE TEMP TABLE interaction_primary AS
SELECT motif_id,
  max(estimate) FILTER (WHERE series_id='saos2' AND isoform='TA')
    AS ta_confirmation_interaction_saos,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='saos2' AND isoform='TA')
    AS ta_confirmation_interaction_saos_q,
  max(estimate) FILTER (WHERE series_id='skmel29_2' AND isoform='TA')
    AS ta_confirmation_interaction_skmel,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA')
    AS ta_confirmation_interaction_skmel_q
FROM tp73_interaction
WHERE negative_reference_threshold=-1 AND evaluation_status='ok'
  AND contrast='cofactor_by_tp73_confirmation_interaction'
GROUP BY motif_id;

CREATE TABLE joint_primary_motif AS
SELECT e.motif_id, e.factor_name, e.positive_threshold,
       e.source_score_floor, e.anchors_positive, e.anchors_negative_reference,
       e.adjusted_odds_ratio AS tp73_adjusted_odds_ratio,
       e.confidence_interval_95_lower AS tp73_ci_lower,
       e.confidence_interval_95_upper AS tp73_ci_upper,
       e.q_value_bh_all_jaspar AS tp73_q,
       e.association_direction AS tp73_association_direction,
       b.* EXCLUDE (motif_id, factor_name, strict_ok_rows, strict_designs),
       g.* EXCLUDE (motif_id), i.* EXCLUDE (motif_id),
       x.* EXCLUDE (motif_id),
       (b.ta_saos+b.ta_skmel)/2 AS ta_mean_effect,
       (b.dn_saos+b.dn_skmel)/2 AS dn_mean_effect,
       CASE WHEN b.ta_saos>0 AND b.ta_skmel>0 THEN 'both_positive'
            WHEN b.ta_saos<0 AND b.ta_skmel<0 THEN 'both_negative'
            ELSE 'opposite_direction' END AS ta_direction,
       CASE WHEN b.dn_saos>0 AND b.dn_skmel>0 THEN 'both_positive'
            WHEN b.dn_saos<0 AND b.dn_skmel<0 THEN 'both_negative'
            ELSE 'opposite_direction' END AS dn_direction,
       b.ta_saos_q<=0.05 AND b.ta_skmel_q<=0.05 AS ta_both_q05,
       b.dn_saos_q<=0.05 AND b.dn_skmel_q<=0.05 AS dn_both_q05
FROM enrichment e JOIN binary_primary b USING (motif_id)
LEFT JOIN gradient_primary g USING (motif_id)
LEFT JOIN intergenic_primary i USING (motif_id)
LEFT JOIN interaction_primary x USING (motif_id)
WHERE e.evaluation_status='ok' AND b.strict_ok_rows=4 AND b.strict_designs=4
ORDER BY e.motif_id;

CREATE TEMP TABLE binary_zero AS
SELECT motif_id, count(*) AS zero_ok_rows,
  count(DISTINCT (series_id,isoform)) AS zero_designs,
  max(estimate) FILTER (WHERE series_id='saos2' AND isoform='TA') AS ta_saos,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='saos2' AND isoform='TA') AS ta_saos_q,
  max(estimate) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_skmel,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='skmel29_2' AND isoform='TA') AS ta_skmel_q,
  max(estimate) FILTER (WHERE series_id='saos2' AND isoform='DN') AS dn_saos,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='saos2' AND isoform='DN') AS dn_saos_q,
  max(estimate) FILTER
    (WHERE series_id='skmel29_2' AND isoform='DN') AS dn_skmel,
  max(q_value_bh_all_motifs) FILTER
    (WHERE series_id='skmel29_2' AND isoform='DN') AS dn_skmel_q
FROM intensity
WHERE negative_reference_threshold=0 AND evaluation_status='ok'
GROUP BY motif_id;

CREATE TABLE joint_reference_zero_motif AS
SELECT e.motif_id, e.factor_name, e.positive_threshold,
       e.source_score_floor, e.evaluation_status AS enrichment_status,
       e.adjusted_odds_ratio AS tp73_adjusted_odds_ratio,
       e.q_value_bh_all_jaspar AS tp73_q,
       e.association_direction AS tp73_association_direction,
       z.* EXCLUDE (motif_id, zero_ok_rows, zero_designs)
FROM enrichment e JOIN binary_zero z USING (motif_id)
WHERE z.zero_ok_rows=4 AND z.zero_designs=4
ORDER BY e.motif_id;

CREATE TABLE context_primary_effect AS
SELECT e.adjusted_odds_ratio AS tp73_adjusted_odds_ratio,
       e.q_value_bh_all_jaspar AS tp73_q,
       e.association_direction AS tp73_association_direction,
       c.*
FROM enrichment e JOIN context_effect c USING (motif_id)
WHERE e.evaluation_status='ok' AND c.negative_reference_threshold=-1
ORDER BY c.motif_id,c.isoform,c.series_id,c.genomic_context_class;

CREATE TABLE gene_relation_primary_effect AS
SELECT e.adjusted_odds_ratio AS tp73_adjusted_odds_ratio,
       e.q_value_bh_all_jaspar AS tp73_q,
       e.association_direction AS tp73_association_direction,
       g.*
FROM enrichment e JOIN gene_relation_effect g USING (motif_id)
WHERE e.evaluation_status='ok' AND g.negative_reference_threshold=-1
ORDER BY g.motif_id,g.isoform,g.series_id,g.gene_relation_class;

CREATE TABLE score_gradient_primary_effect AS
SELECT e.adjusted_odds_ratio AS tp73_adjusted_odds_ratio,
       e.q_value_bh_all_jaspar AS tp73_q,
       e.association_direction AS tp73_association_direction,
       s.*
FROM enrichment e JOIN score_gradient s USING (motif_id)
WHERE e.evaluation_status='ok' AND s.score_clamp_reference=-1
ORDER BY s.motif_id,s.isoform,s.series_id;

CREATE TABLE tp73_interaction_primary_effect AS
SELECT e.adjusted_odds_ratio AS tp73_adjusted_odds_ratio,
       e.q_value_bh_all_jaspar AS tp73_q,
       e.association_direction AS tp73_association_direction,
       t.*
FROM enrichment e JOIN tp73_interaction t USING (motif_id)
WHERE e.evaluation_status='ok' AND t.negative_reference_threshold=-1
ORDER BY t.motif_id,t.isoform,t.series_id,t.contrast;

CREATE TABLE summary_metrics AS
SELECT
  (SELECT count(*) FROM motif_coverage) AS planned_motifs,
  (SELECT count(*) FROM motif_coverage
    WHERE enrichment_status='ok') AS enrichment_estimable,
  (SELECT count(*) FROM motif_coverage
    WHERE enrichment_status='ok' AND q_value_bh_all_jaspar<=0.05
      AND association_direction='anti_p73_enriched') AS tp73_enriched_q05,
  (SELECT count(*) FROM motif_coverage
    WHERE enrichment_status='ok' AND q_value_bh_all_jaspar<=0.05
      AND association_direction='anti_p73_depleted') AS tp73_depleted_q05,
  (SELECT count(*) FROM joint_primary_motif) AS strict_primary_motifs,
  (SELECT count(*) FROM joint_reference_zero_motif) AS reference_zero_motifs,
  (SELECT count(*) FROM joint_primary_motif
    WHERE ta_both_q05) AS ta_both_series_q05,
  (SELECT count(*) FROM joint_primary_motif
    WHERE ta_both_q05 AND ta_direction<>'opposite_direction')
    AS ta_both_series_q05_same_direction,
  (SELECT count(*) FROM joint_primary_motif
    WHERE dn_both_q05) AS dn_both_series_q05,
  (SELECT corr(ln(tp73_adjusted_odds_ratio),ta_saos)
    FROM joint_primary_motif WHERE tp73_q<=0.05) AS tp73_logor_ta_saos_correlation,
  (SELECT corr(ln(tp73_adjusted_odds_ratio),ta_skmel)
    FROM joint_primary_motif WHERE tp73_q<=0.05) AS tp73_logor_ta_skmel_correlation,
  (SELECT corr(ln(tp73_adjusted_odds_ratio),dn_saos)
    FROM joint_primary_motif WHERE tp73_q<=0.05) AS tp73_logor_dn_saos_correlation,
  (SELECT corr(ln(tp73_adjusted_odds_ratio),dn_skmel)
    FROM joint_primary_motif WHERE tp73_q<=0.05) AS tp73_logor_dn_skmel_correlation;

{os.linesep.join(copies)}
"""


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--enrichment", type=Path, required=True,
                        help="final cofactor_primary_occupancy Parquet")
    result.add_argument("--intensity", type=Path, required=True,
                        help="final H3K4me3 intensity_effect Parquet")
    result.add_argument("--context", type=Path, required=True,
                        help="final context-stratified effect Parquet")
    result.add_argument("--gene-relation", type=Path, required=True,
                        help="final four-way gene-relation effect Parquet")
    result.add_argument("--score-gradient", type=Path, required=True,
                        help="final score_gradient Parquet")
    result.add_argument("--interaction", type=Path, required=True,
                        help="final TP73-interaction Parquet")
    result.add_argument("--h3-summary", type=Path, required=True,
                        help="final H3K4me3 chromosome-summary Parquet")
    result.add_argument("--output-dir", type=Path, required=True)
    result.add_argument("--duckdb", type=Path, default=Path("duckdb"))
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        duckdb = arguments.duckdb
        if not duckdb.is_absolute():
            located = shutil.which(str(duckdb))
            duckdb = Path(located) if located else duckdb
        if not duckdb.is_file() or not os.access(duckdb, os.X_OK):
            raise SummaryError(f"DuckDB executable not found: {arguments.duckdb}")

        inputs = {
            "enrichment": arguments.enrichment.expanduser().resolve(),
            "intensity": arguments.intensity.expanduser().resolve(),
            "context": arguments.context.expanduser().resolve(),
            "gene_relation": arguments.gene_relation.expanduser().resolve(),
            "score_gradient": arguments.score_gradient.expanduser().resolve(),
            "interaction": arguments.interaction.expanduser().resolve(),
            "h3_summary": arguments.h3_summary.expanduser().resolve(),
        }
        for label, path in inputs.items():
            if not path.is_file():
                raise SummaryError(f"{label} input is missing: {path}")

        output = arguments.output_dir.expanduser().resolve()
        if output.exists():
            raise SummaryError(f"output directory already exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = output.with_name(f".{output.name}.staging-{os.getpid()}")
        if staging.exists():
            raise SummaryError(f"staging directory already exists: {staging}")
        staging.mkdir()

        run_sql(duckdb, output_sql(inputs, staging))
        outputs = {}
        for table in OUTPUT_TABLES:
            parquet = staging / f"{table}.parquet"
            tsv = staging / f"{table}.tsv"
            if not parquet.is_file() or not tsv.is_file():
                raise SummaryError(f"DuckDB omitted output table: {table}")
            outputs[parquet.name] = {
                "bytes": parquet.stat().st_size,
                "rows": query_scalar(duckdb, parquet),
                "sha256": sha256(parquet),
            }
            outputs[tsv.name] = {
                "bytes": tsv.stat().st_size,
                "sha256": sha256(tsv),
            }

        source = Path(__file__).resolve().parent.parent
        commit, dirty = git_identity(source)
        manifest = {
            "schema_version": 2,
            "analysis": "joint_tp73_enrichment_h3k4me3_cofactor_summary",
            "analysis_partition": "autosome",
            "included_chromosomes": [str(value) for value in range(1, 23)],
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "primary_negative_reference": -1,
            "historical_negative_reference": 0,
            "source_commit": commit,
            "source_dirty": dirty,
            "inputs": {
                label: {
                    "path": str(path), "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for label, path in inputs.items()
            },
            "outputs": outputs,
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, output)
        print(f"I: wrote joint genome-wide cofactor summary: {output}", file=sys.stderr)
        return 0
    except (SummaryError, OSError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
