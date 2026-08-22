#!/usr/bin/env bash

set -euo pipefail

if ! command -v duckdb >/dev/null 2>&1; then
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
fi

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/tp73-distance-manager.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

scan="$temporary/scan"
catalog="$temporary/catalog"
run_root="$temporary/run"
scratch="$temporary/scratch"
mkdir -p "$scan/tables/jaspar2026/scan_file_inventory" \
    "$scan/tables/jaspar2026/motif_metadata" "$scan/hits" "$catalog" "$scratch"

duckdb -light-mode -batch :memory: >/dev/null <<SQL
COPY (
  SELECT * FROM (VALUES
    ('MA0001.1', 'HUMAN_TEST', 9::BIGINT),
    ('MA0002.1', 'MOUSE_TEST', 9::BIGINT),
    ('MA0003.1', 'PLANT_TEST', 9::BIGINT)
  ) AS v(motif_id, motif_name, motif_length)
) TO '$scan/tables/jaspar2026/motif_metadata/part-000000.parquet'
  (FORMAT PARQUET);

COPY (
  WITH anchor AS (
    SELECT i, (100 + 200 * i)::BIGINT AS anchor_start,
           (116 + 200 * i)::BIGINT AS anchor_end,
           CASE WHEN i < 4 THEN 'positive'
                WHEN i = 4 THEN 'intermediate' ELSE 'negative' END AS class
    FROM range(8) r(i)
  )
  SELECT '1'::VARCHAR AS chrom, anchor_start, anchor_end,
         3.0::FLOAT AS anchor_score,
         (class = 'positive' AND i <> 3 OR class = 'negative' AND i = 5)
           AS supported_tp73_saos2_TA,
         (class = 'positive' AND i = 3 OR class = 'negative' AND i <> 5)
           AS supported_negative_control_saos2_TA,
         (class = 'positive' AND i = 3 OR class = 'negative' AND i <> 7)
           AS supported_tp73_saos2_DN,
         (class = 'positive' AND i <> 3 OR class = 'negative' AND i = 7)
           AS supported_negative_control_saos2_DN,
         (class = 'positive' AND i <> 3 OR class = 'negative' AND i = 5)
           AS supported_tp73_skmel29_2_TA,
         (class = 'positive' AND i = 3 OR class = 'negative' AND i <> 5)
           AS supported_negative_control_skmel29_2_TA,
         (class = 'positive' AND i = 3 OR class = 'negative' AND i <> 7)
           AS supported_tp73_skmel29_2_DN,
         (class = 'positive' AND i <> 3 OR class = 'negative' AND i = 7)
           AS supported_negative_control_skmel29_2_DN
  FROM anchor
) TO '$temporary/anchors.parquet' (FORMAT PARQUET);

COPY (
  SELECT (116 + 200 * i)::BIGINT AS start,
         (125 + 200 * i)::BIGINT AS "end",
         CASE WHEN i < 4 THEN 3.0 ELSE 0.0 END::FLOAT AS score
  FROM range(5) r(i)
) TO '$temporary/hits-plus.parquet' (FORMAT PARQUET);
COPY (
  SELECT 116::BIGINT AS start, 125::BIGINT AS "end", 2.5::FLOAT AS score
) TO '$temporary/hits-minus.parquet' (FORMAT PARQUET);

COPY (
  SELECT * FROM (VALUES
    ('synthetic_thresholds', 'MA0001.1', 2.0::DOUBLE, -1.0::DOUBLE, 'synthetic'),
    ('synthetic_thresholds', 'MA0002.1', 2.0::DOUBLE, -1.0::DOUBLE, 'synthetic'),
    ('synthetic_thresholds', 'MA0003.1', 2.0::DOUBLE, -1.0::DOUBLE, 'synthetic')
  ) AS v(threshold_set_id, motif_id, recommended_threshold,
         source_minimum_score, calibration_status)
) TO '$temporary/thresholds.parquet' (FORMAT PARQUET);
SQL

for motif in MA0001.1 MA0002.1 MA0003.1; do
    mkdir -p "$scan/hits/$motif"
    cp "$temporary/hits-plus.parquet" "$scan/hits/$motif/plus.parquet"
    cp "$temporary/hits-minus.parquet" "$scan/hits/$motif/minus.parquet"
done

python3 - "$scan" > "$temporary/inventory.tsv" <<'PY'
import hashlib
from pathlib import Path
import sys

root = Path(sys.argv[1])
print("motif_id\tchrom\tstrand\toutput_relative_path\tbytes\tsha256\temitted_hits\tminimum_score")
for motif in ("MA0001.1", "MA0002.1", "MA0003.1"):
    for strand, name in (("+", "plus.parquet"), ("-", "minus.parquet")):
        relative = Path("hits") / motif / name
        path = root / relative
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        print(motif, "1", strand, relative, path.stat().st_size, digest, 5, -1,
              sep="\t")
PY

duckdb -light-mode -batch :memory: >/dev/null <<SQL
COPY (
  SELECT motif_id::VARCHAR AS motif_id, chrom::VARCHAR AS chrom,
         strand::VARCHAR AS strand,
         output_relative_path::VARCHAR AS output_relative_path,
         bytes::BIGINT AS bytes, sha256::VARCHAR AS sha256,
         emitted_hits::BIGINT AS emitted_hits,
         minimum_score::DOUBLE AS minimum_score
  FROM read_csv_auto('$temporary/inventory.tsv', delim='\t', header=true)
) TO '$scan/tables/jaspar2026/scan_file_inventory/part-000000.parquet'
  (FORMAT PARQUET);
SQL
printf '{}\n' > "$scan/manifest.json"

duckdb -light-mode -batch "$catalog/jaspar_metadata.duckdb" >/dev/null <<SQL
CREATE TABLE jaspar_matrix AS
SELECT * FROM (VALUES
  ('CORE', 'vertebrates', 'MA0001.1', 'MA0001', 1, 'HUMAN_TEST',
   'Class A', 'Family A', true, false, false, 1::BIGINT),
  ('CORE', 'vertebrates', 'MA0002.1', 'MA0002', 1, 'MOUSE_TEST',
   'Class B', 'Family B', false, true, false, 1::BIGINT),
  ('CORE', 'plants', 'MA0003.1', 'MA0003', 1, 'PLANT_TEST',
   'Class C', 'Family C', false, false, false, 1::BIGINT)
) AS v(collection, tax_group, matrix_id, base_id, version, name,
       class, family, includes_homo_sapiens, includes_mus_musculus,
       includes_rattus_norvegicus, species_count);
CREATE TABLE jaspar_matrix_species AS
SELECT * FROM (VALUES
  ('MA0001.1', 1::BIGINT, 9606::BIGINT, 'Homo sapiens'),
  ('MA0002.1', 1::BIGINT, 10090::BIGINT, 'Mus musculus'),
  ('MA0003.1', 1::BIGINT, 3702::BIGINT, 'Arabidopsis thaliana')
) AS v(matrix_id, species_ordinal, tax_id, species_scientific_name);
CREATE TABLE jaspar_motif_set_matrix AS
SELECT 'synthetic_set'::VARCHAR AS motif_set_id, matrix_id,
       row_number() OVER (ORDER BY matrix_id)::BIGINT AS motif_order
FROM jaspar_matrix;
CHECKPOINT;
SQL
printf '{}\n' > "$catalog/manifest.json"

if "$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" prepare \
    --run-root "$temporary/failed-prepare" --scan-package "$scan" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --thresholds "$temporary/thresholds.parquet" \
    --threshold-set-id synthetic_thresholds --jaspar-catalog "$catalog" \
    --source "$repository_root" \
    --source-commit aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa \
    --run-id synthetic_failed_prepare --tax-group vertebrates \
    --chromosomes 1 --duckdb false >/dev/null 2>&1; then
    echo "E: deliberately failed prepare unexpectedly succeeded" >&2
    exit 1
fi
[[ ! -e $temporary/failed-prepare ]] || {
    echo "E: failed prepare published a partial run root" >&2
    exit 1
}

task_count=$(
    "$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" prepare \
        --run-root "$run_root" --scan-package "$scan" \
        --anchor-evidence "$temporary/anchors.parquet" \
        --thresholds "$temporary/thresholds.parquet" \
        --threshold-set-id synthetic_thresholds --jaspar-catalog "$catalog" \
        --source "$repository_root" \
        --source-commit aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa \
        --run-id synthetic_distance_species_v1 --tax-group vertebrates \
        --chromosomes 1
)
[[ $task_count == 2 ]] || {
    echo "E: expected two vertebrate motif tasks, found $task_count" >&2
    exit 1
}
grep -q $'^0\tMA0001.1\tHUMAN_TEST' "$run_root/plan/tasks.tsv"
grep -q $'^1\tMA0002.1\tMOUSE_TEST' "$run_root/plan/tasks.tsv"
if grep -q 'MA0003.1' "$run_root/plan/tasks.tsv"; then
    echo "E: plant motif entered the vertebrate task family" >&2
    exit 1
fi

"$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" \
    prepare-anchors --run-root "$run_root" --threads 1 --memory-limit 1GB
"$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" \
    prepare-anchors --run-root "$run_root" --threads 1 --memory-limit 1GB
for task in 0 1; do
    "$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" run-task \
        --run-root "$run_root" --task-index "$task" --scratch "$scratch" \
        --threads 1 --memory-limit 1GB --max-temp-size 1GB
done
[[ -f $run_root/checkpoints/task-000000-MA0001.1/chrom-1/complete.json ]] || {
    echo "E: chromosome checkpoint was not published" >&2
    exit 1
}
rm -rf "$run_root/tasks/task-000000-MA0001.1"
resume_log=$(
    "$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" run-task \
        --run-root "$run_root" --task-index 0 --scratch "$scratch" \
        --threads 1 --memory-limit 1GB --max-temp-size 1GB 2>&1
)
grep -q 'Reusing chromosome checkpoint' <<<"$resume_log"

status=$(
    "$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" \
        status --run-root "$run_root"
)
grep -q $'planned\t2' <<<"$status"
grep -q $'complete\t2' <<<"$status"

"$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" finalize \
    --run-root "$run_root"
"$repository_root/scripts/manage_tp73_distance_cofactor_enrichment.py" finalize \
    --run-root "$run_root"

database="$run_root/final/distance_enrichment/tp73_distance_cofactor_enrichment.duckdb"
duckdb -light-mode -batch "$database" >/dev/null <<SQL
SELECT CASE WHEN (SELECT count(*) FROM cofactor_distance_enrichment) <> 24
  THEN error('result cardinality is wrong') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM cofactor_distance_enrichment WHERE tax_group <> 'vertebrates'
) THEN error('non-vertebrate motif entered the primary result') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM cofactor_distance_enrichment
  WHERE motif_id='MA0001.1' AND includes_homo_sapiens
    AND source_species='Homo sapiens'
) THEN error('human source-species annotation is absent') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM cofactor_distance_enrichment
  WHERE motif_id='MA0002.1' AND NOT includes_homo_sapiens
    AND includes_mus_musculus AND NOT includes_rattus_norvegicus
    AND source_species='Mus musculus'
) THEN error('mouse source-species annotation is absent') END;
SELECT CASE WHEN (SELECT count(*) FROM jaspar_matrix_species) <> 3
  THEN error('normalized species dimension was not retained') END;
SELECT CASE WHEN (SELECT count(*) FROM jaspar_motif_set_matrix) <> 3
  THEN error('motif-set membership was not retained') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM information_schema.tables
  WHERE table_name='cofactor_distance_top_bottom_20_human_source'
) THEN error('exact-human sensitivity ranking table is absent') END;
SQL

python3 - "$run_root" <<'PY'
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
config = json.loads((root / "plan/run_config.json").read_text())
manifest = json.loads((root / "final/distance_enrichment/manifest.json").read_text())
assert config["chromosomes"] == ["1"]
assert config["tax_group"] == "vertebrates"
assert config["human_source_task_count"] == 1
assert config["source_species_unspecified_task_count"] == 0
assert "not a binding-species restriction" in config["source_species_semantics"]
assert manifest["chromosomes"] == ["1"]
assert manifest["tax_group"] == "vertebrates"
PY

echo "TP73 distance cofactor enrichment manager tests passed."
