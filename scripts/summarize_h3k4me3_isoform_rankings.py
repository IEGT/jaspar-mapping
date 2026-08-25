#!/usr/bin/env python3

"""Build score-zero human TF rankings for GFP-referenced H3K4me3 change."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path


class RankingError(RuntimeError):
    pass


OUTPUTS = (
    "human_score0_effect_matrix",
    "human_score0_positive_top",
    "human_score0_negative_top",
    "human_score0_isoform_difference_top",
)


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_duckdb(duckdb: Path, database: Path, sql: str) -> None:
    process = subprocess.run(
        [str(duckdb), "-light-mode", "-readonly", "-bail", str(database)],
        cwd=database.parent, input=sql, text=True, capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise RankingError(process.stderr.strip() or "DuckDB ranking query failed")


def query_count(duckdb: Path, parquet: Path) -> int:
    process = subprocess.run(
        [
            str(duckdb), "-light-mode", "-csv", "-noheader",
            ":memory:", "-c",
            f"SELECT count(*) FROM read_parquet({sql_string(parquet)});",
        ],
        text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise RankingError(process.stderr.strip() or "cannot count ranking output")
    return int(process.stdout.strip())


def ranking_sql(metadata: Path, staging: Path, top: int) -> str:
    copies = []
    for table in OUTPUTS:
        copies.extend((
            f"COPY {table} TO {sql_string(staging / (table + '.parquet'))} "
            "(FORMAT PARQUET, COMPRESSION ZSTD);",
            f"COPY {table} TO {sql_string(staging / (table + '.tsv'))} "
            "(FORMAT CSV, HEADER true, DELIMITER '\\t', NULL 'NA');",
        ))
    return f"""
SET preserve_insertion_order=false;
CREATE TEMP VIEW jaspar_metadata_raw AS
SELECT * FROM read_csv_auto(
  {sql_string(metadata)}, delim='\t', header=true, all_varchar=true,
  sample_size=-1, nullstr=''
);
CREATE TEMP TABLE jaspar_metadata AS
SELECT matrix_id::VARCHAR AS motif_id,
       min(name)::VARCHAR AS jaspar_name,
       min(tax_group)::VARCHAR AS tax_group,
       bool_or(regexp_matches(tax_id, '(^|;)9606(;|$)'))
         AS includes_homo_sapiens,
       string_agg(DISTINCT species, '; ' ORDER BY species)
         FILTER (WHERE species IS NOT NULL AND species <> '') AS source_species
FROM jaspar_metadata_raw
GROUP BY matrix_id;

CREATE TEMP TABLE intensity_wide AS
SELECT i.motif_id, min(i.factor_name)::VARCHAR AS factor_name,
       i.distance_band, min(i.positive_threshold)::DOUBLE AS positive_threshold,
       min(i.source_score_floor)::DOUBLE AS source_score_floor,
       max(i.estimate) FILTER (WHERE i.series_id='saos2' AND i.isoform='TA')
         AS ta_saos_effect,
       max(i.q_value_bh_all_motifs)
         FILTER (WHERE i.series_id='saos2' AND i.isoform='TA') AS ta_saos_q,
       max(i.adjusted_mean_change_negative)
         FILTER (WHERE i.series_id='saos2' AND i.isoform='TA')
         AS ta_saos_adjusted_negative,
       max(i.adjusted_mean_change_positive)
         FILTER (WHERE i.series_id='saos2' AND i.isoform='TA')
         AS ta_saos_adjusted_positive,
       max(i.estimate) FILTER (WHERE i.series_id='skmel29_2' AND i.isoform='TA')
         AS ta_skmel_effect,
       max(i.q_value_bh_all_motifs)
         FILTER (WHERE i.series_id='skmel29_2' AND i.isoform='TA') AS ta_skmel_q,
       max(i.adjusted_mean_change_negative)
         FILTER (WHERE i.series_id='skmel29_2' AND i.isoform='TA')
         AS ta_skmel_adjusted_negative,
       max(i.adjusted_mean_change_positive)
         FILTER (WHERE i.series_id='skmel29_2' AND i.isoform='TA')
         AS ta_skmel_adjusted_positive,
       max(i.estimate) FILTER (WHERE i.series_id='saos2' AND i.isoform='DN')
         AS dn_saos_effect,
       max(i.q_value_bh_all_motifs)
         FILTER (WHERE i.series_id='saos2' AND i.isoform='DN') AS dn_saos_q,
       max(i.adjusted_mean_change_negative)
         FILTER (WHERE i.series_id='saos2' AND i.isoform='DN')
         AS dn_saos_adjusted_negative,
       max(i.adjusted_mean_change_positive)
         FILTER (WHERE i.series_id='saos2' AND i.isoform='DN')
         AS dn_saos_adjusted_positive,
       max(i.estimate) FILTER (WHERE i.series_id='skmel29_2' AND i.isoform='DN')
         AS dn_skmel_effect,
       max(i.q_value_bh_all_motifs)
         FILTER (WHERE i.series_id='skmel29_2' AND i.isoform='DN') AS dn_skmel_q,
       max(i.adjusted_mean_change_negative)
         FILTER (WHERE i.series_id='skmel29_2' AND i.isoform='DN')
         AS dn_skmel_adjusted_negative,
       max(i.adjusted_mean_change_positive)
         FILTER (WHERE i.series_id='skmel29_2' AND i.isoform='DN')
         AS dn_skmel_adjusted_positive,
       min(i.anchors_positive)::BIGINT AS minimum_positive_anchors,
       min(i.anchors_negative)::BIGINT AS minimum_negative_anchors,
       count(*)::BIGINT AS estimable_designs
FROM intensity_effect i
WHERE i.negative_reference_threshold=-1
  AND i.positive_threshold=0
  AND i.evaluation_status='ok'
GROUP BY i.motif_id, i.distance_band
HAVING count(*)=4;

CREATE TEMP TABLE contrast_wide AS
SELECT c.motif_id, c.distance_band,
       max(c.estimate) FILTER (WHERE c.series_id='saos2') AS isoform_saos_effect,
       max(c.q_value_bh_all_motifs) FILTER (WHERE c.series_id='saos2')
         AS isoform_saos_q,
       max(c.estimate) FILTER (WHERE c.series_id='skmel29_2')
         AS isoform_skmel_effect,
       max(c.q_value_bh_all_motifs) FILTER (WHERE c.series_id='skmel29_2')
         AS isoform_skmel_q,
       count(*)::BIGINT AS estimable_contrasts
FROM isoform_contrast c
WHERE c.negative_reference_threshold=-1
  AND c.positive_threshold=0
  AND c.contrast='TA_minus_DN'
  AND c.evaluation_status='ok'
GROUP BY c.motif_id, c.distance_band
HAVING count(*)=2;

CREATE TEMP TABLE human_score0_effect_matrix AS
SELECT i.motif_id, coalesce(m.jaspar_name, i.factor_name) AS factor_name,
       i.distance_band, i.positive_threshold, i.source_score_floor,
       m.tax_group, m.includes_homo_sapiens, m.source_species,
       i.minimum_positive_anchors, i.minimum_negative_anchors,
       i.ta_saos_effect, i.ta_saos_q, i.ta_skmel_effect, i.ta_skmel_q,
       (i.ta_saos_effect+i.ta_skmel_effect)/2 AS ta_macro_effect,
       (i.ta_saos_adjusted_negative+i.ta_skmel_adjusted_negative)/2
         AS ta_macro_adjusted_negative,
       (i.ta_saos_adjusted_positive+i.ta_skmel_adjusted_positive)/2
         AS ta_macro_adjusted_positive,
       i.dn_saos_effect, i.dn_saos_q, i.dn_skmel_effect, i.dn_skmel_q,
       (i.dn_saos_effect+i.dn_skmel_effect)/2 AS dn_macro_effect,
       (i.dn_saos_adjusted_negative+i.dn_skmel_adjusted_negative)/2
         AS dn_macro_adjusted_negative,
       (i.dn_saos_adjusted_positive+i.dn_skmel_adjusted_positive)/2
         AS dn_macro_adjusted_positive,
       c.isoform_saos_effect, c.isoform_saos_q,
       c.isoform_skmel_effect, c.isoform_skmel_q,
       (c.isoform_saos_effect+c.isoform_skmel_effect)/2
         AS isoform_macro_effect,
       CASE WHEN c.isoform_saos_effect>0 AND c.isoform_skmel_effect>0
              THEN 'TA_stronger_both_systems'
            WHEN c.isoform_saos_effect<0 AND c.isoform_skmel_effect<0
              THEN 'DN_stronger_both_systems'
            ELSE 'cell_system_heterogeneous' END AS isoform_direction,
       CASE WHEN (i.ta_saos_effect+i.ta_skmel_effect)/2 > 0
                  AND (i.dn_saos_effect+i.dn_skmel_effect)/2 < 0
              THEN 'TA_positive_DN_negative'
            WHEN (i.ta_saos_effect+i.ta_skmel_effect)/2 < 0
                  AND (i.dn_saos_effect+i.dn_skmel_effect)/2 > 0
              THEN 'TA_negative_DN_positive'
            ELSE 'not_opposite_macro_sign' END AS opposite_isoform_sign,
       greatest(
         (i.ta_saos_effect+i.ta_skmel_effect)/2,
         (i.dn_saos_effect+i.dn_skmel_effect)/2
       ) AS positive_rank_effect,
       least(
         (i.ta_saos_effect+i.ta_skmel_effect)/2,
         (i.dn_saos_effect+i.dn_skmel_effect)/2
       ) AS negative_rank_effect,
       abs((c.isoform_saos_effect+c.isoform_skmel_effect)/2)
         AS isoform_difference_rank_effect
FROM intensity_wide i
JOIN contrast_wide c USING (motif_id, distance_band)
JOIN jaspar_metadata m USING (motif_id)
WHERE m.tax_group='vertebrates' AND m.includes_homo_sapiens
ORDER BY i.distance_band, factor_name, i.motif_id;

CREATE TEMP TABLE ranked_positive AS
SELECT *, row_number() OVER (
  PARTITION BY distance_band, factor_name
  ORDER BY positive_rank_effect DESC, motif_id
) AS factor_motif_rank
FROM human_score0_effect_matrix
WHERE positive_rank_effect>0;
CREATE TEMP TABLE human_score0_positive_top AS
SELECT * EXCLUDE (factor_motif_rank, overall_rank), overall_rank AS rank
FROM (
  SELECT *, row_number() OVER (
    PARTITION BY distance_band
    ORDER BY positive_rank_effect DESC, factor_name, motif_id
  ) AS overall_rank
  FROM ranked_positive WHERE factor_motif_rank=1
)
WHERE overall_rank<={top}
ORDER BY distance_band, rank;

CREATE TEMP TABLE ranked_negative AS
SELECT *, row_number() OVER (
  PARTITION BY distance_band, factor_name
  ORDER BY negative_rank_effect ASC, motif_id
) AS factor_motif_rank
FROM human_score0_effect_matrix
WHERE negative_rank_effect<0;
CREATE TEMP TABLE human_score0_negative_top AS
SELECT * EXCLUDE (factor_motif_rank, overall_rank), overall_rank AS rank
FROM (
  SELECT *, row_number() OVER (
    PARTITION BY distance_band
    ORDER BY negative_rank_effect ASC, factor_name, motif_id
  ) AS overall_rank
  FROM ranked_negative WHERE factor_motif_rank=1
)
WHERE overall_rank<={top}
ORDER BY distance_band, rank;

CREATE TEMP TABLE ranked_isoform AS
SELECT *, row_number() OVER (
  PARTITION BY distance_band, factor_name
  ORDER BY isoform_difference_rank_effect DESC, motif_id
) AS factor_motif_rank
FROM human_score0_effect_matrix
WHERE isoform_difference_rank_effect>0;
CREATE TEMP TABLE human_score0_isoform_difference_top AS
SELECT * EXCLUDE (factor_motif_rank, overall_rank), overall_rank AS rank
FROM (
  SELECT *, row_number() OVER (
    PARTITION BY distance_band
    ORDER BY isoform_difference_rank_effect DESC, factor_name, motif_id
  ) AS overall_rank
  FROM ranked_isoform WHERE factor_motif_rank=1
)
WHERE overall_rank<={top}
ORDER BY distance_band, rank;

{chr(10).join(copies)}
"""


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--analysis-database", type=Path, required=True)
    result.add_argument("--jaspar-metadata", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    result.add_argument("--duckdb", type=Path, default=Path("duckdb"))
    result.add_argument("--top", type=int, default=25)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        database = arguments.analysis_database.expanduser().resolve()
        metadata = arguments.jaspar_metadata.expanduser().resolve()
        output = arguments.output_dir.expanduser().resolve()
        duckdb = arguments.duckdb.expanduser()
        if not duckdb.is_absolute():
            found = shutil.which(str(duckdb))
            if found is None:
                raise RankingError(f"DuckDB is unavailable: {duckdb}")
            duckdb = Path(found)
        duckdb = duckdb.resolve()
        if not database.is_file() or not metadata.is_file():
            raise RankingError("analysis database or JASPAR metadata is missing")
        if arguments.top < 1:
            raise RankingError("--top must be positive")
        if output.exists():
            raise RankingError(f"output directory already exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
        try:
            run_duckdb(
                duckdb, database,
                ranking_sql(metadata, staging, arguments.top),
            )
            row_counts = {
                table: query_count(duckdb, staging / f"{table}.parquet")
                for table in OUTPUTS
            }
            if row_counts["human_score0_effect_matrix"] == 0:
                raise RankingError("no estimable exact-human score-zero motifs")
            manifest = {
                "schema_version": 1,
                "state": "complete",
                "analysis": "human_score0_h3k4me3_isoform_rankings",
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "analysis_database": str(database),
                "analysis_database_sha256": sha256(database),
                "jaspar_metadata": str(metadata),
                "jaspar_metadata_sha256": sha256(metadata),
                "positive_threshold": 0,
                "negative_reference_rule": "context_score < -1 or absent",
                "intermediate_exclusion": "-1 <= context_score < 0",
                "source_scope": "tax_group=vertebrates AND includes_homo_sapiens",
                "tf_consolidation": (
                    "one strongest representative motif per factor and criterion"
                ),
                "top_per_distance_band": arguments.top,
                "row_counts": row_counts,
            }
            (staging / "manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            os.replace(staging, output)
        finally:
            if staging.exists():
                shutil.rmtree(staging)
        print(f"I: wrote human score-zero H3K4me3 rankings: {output}",
              file=sys.stderr)
        return 0
    except (RankingError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
