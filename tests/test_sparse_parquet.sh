#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
scanner="$repository_root/pssm_scan_parquet"
fixture_dir="$repository_root/test_files/synthetic_dense"
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-sparse-parquet.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required." >&2
    exit 1
}
[[ -x $scanner ]] || {
    echo "E: Build the Arrow scanner with 'make pssm_scan_parquet'." >&2
    exit 1
}

fasta_index="$temporary/genome.fna.fai"
python3 "$repository_root/scripts/build_fasta_index.py" \
    "$fixture_dir/genome.fna" --output "$fasta_index" >/dev/null
printf 'MA9999.1\n' > "$temporary/motifs.txt"

common=(
    --sparse-parquet
    --genome "$fixture_dir/genome.fna"
    --fasta-index "$fasta_index"
    --pssm "$fixture_dir/JASPAR2026_synthetic.jaspar"
    --motif-list "$temporary/motifs.txt"
    --motif-set-id synthetic_jaspar2026
    --genome-id synthetic_acgtn_v1
    --chr 1
    --strand both
    --coordinate-mode bed
    --score-mode log2_relative_risk
    --pseudocount 0
    --skip-N
)

included="$temporary/nested/output/included"
omitted="$temporary/omitted"
included_task="$included/task_data/task_id=synthetic-threshold-0"
omitted_task="$omitted/task_data/task_id=synthetic-threshold-m1"
"$scanner" "${common[@]}" --threshold 0 --show-sequence \
    --scan-file-stats "$temporary/included-stats.jsonl" \
    --outdir "$included_task" >/dev/null
"$scanner" "${common[@]}" --threshold -1 \
    --outdir "$omitted_task" >/dev/null

included_glob="$included/task_data/*/tables/jaspar2026/motif_hit/**/*.parquet"
omitted_glob="$omitted/task_data/*/tables/jaspar2026/motif_hit/**/*.parquet"
included_plus=$(find "$included" -path '*/strand=plus/*.parquet' -print -quit)
[[ -n $included_plus ]] || {
    echo "E: Sparse writer did not create the plus-strand partition." >&2
    exit 1
}

metadata_dir="$included/tables/jaspar2026/motif_metadata"
mkdir -p "$metadata_dir"
duckdb :memory: -bail -c "COPY (
    SELECT 'MA9999.1'::VARCHAR AS motif_id,
           'synthetic_jaspar2026'::VARCHAR AS motif_set_id,
           'A_RICH'::VARCHAR AS motif_name,
           1::INTEGER AS motif_length,
           2026::INTEGER AS jaspar_version,
           'synthetic-fixture'::VARCHAR AS source_sha256
) TO '$metadata_dir/part-000000.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);" >/dev/null

duckdb :memory: -bail -c "
    CREATE VIEW included AS
    SELECT * FROM read_parquet('$included_glob', hive_partitioning=1);
    CREATE VIEW omitted AS
    SELECT * FROM read_parquet('$omitted_glob', hive_partitioning=1);

    SELECT CASE WHEN (SELECT count(*) FROM included) <> 2
        THEN error('threshold-zero sparse row count is incorrect') END;
    SELECT CASE WHEN EXISTS (
        SELECT 1 FROM included
        WHERE \"end\" <> start + 1
           OR abs(score - 1.0) > 1e-6
           OR abs(pwm_relative_score - 1.0) > 1e-6
           OR matched_seq <> 'A'
           OR motif_id <> 'MA9999.1'
           OR motif_set_id <> 'synthetic_jaspar2026'
           OR genome_id <> 'synthetic_acgtn_v1'
           OR score_mode <> 'log2_relative_risk'
           OR pseudocount <> 0
           OR background_model_id <> 'uniform_acgt_v1'
           OR pseudocount_scheme <> 'additive_per_base'
           OR minimum_score <> 0
           OR chrom <> 1
           OR n_policy <> 'skip'
           OR matched_sequence <> 'included'
    ) THEN error('threshold-zero sparse values or partitions are incorrect') END;
    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM included WHERE start = 0 AND strand = 'plus'
    ) OR NOT EXISTS (
        SELECT 1 FROM included WHERE start = 3 AND strand = 'minus'
    ) THEN error('sparse strand coordinates are incorrect') END;

    SELECT CASE WHEN (SELECT count(*) FROM omitted) <> 8
        THEN error('threshold-minus-one sparse row count is incorrect') END;
    SELECT CASE WHEN EXISTS (
        SELECT 1 FROM omitted
        WHERE matched_seq IS NOT NULL
           OR score < -1
           OR minimum_score <> -1
           OR matched_sequence <> 'omitted'
    ) THEN error('sequence-omitted sparse values are incorrect') END;

    SELECT CASE WHEN EXISTS (
        SELECT 1 FROM parquet_schema('$included_plus')
        WHERE name IN ('motif_set_id', 'genome_id', 'motif_id', 'score_mode',
                       'pseudocount', 'background_model_id',
                       'pseudocount_scheme', 'chrom', 'strand')
    ) THEN error('row-constant identity leaked into the physical schema') END;
    SELECT CASE WHEN (
        SELECT count(*) FROM parquet_schema('$included_plus')
        WHERE name IN ('start', 'end', 'score', 'pwm_relative_score', 'matched_seq')
    ) <> 5 THEN error('physical sparse schema is incomplete') END;
" >/dev/null

python3 - "$temporary/included-stats.jsonl" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as stream:
    records = [json.loads(line) for line in stream if line.strip()]
assert len(records) == 2, records
assert {record["strand"] for record in records} == {"+", "-"}
for record in records:
    assert record["motif_set_id"] == "synthetic_jaspar2026"
    assert record["genome_id"] == "synthetic_acgtn_v1"
    assert record["score_mode"] == "log2_relative_risk"
    assert record["pseudocount"] == 0
    assert record["background_model_id"] == "uniform_acgt_v1"
    assert record["pseudocount_scheme"] == "additive_per_base"
    assert record["minimum_score"] == 0
    assert record["minimum_pwm_relative_score"] is None
    assert record["maximum_pwm_relative_score"] is None
    assert record["coordinate_mode"] == "bed"
    assert record["n_policy"] == "skip"
    assert record["matched_sequence_policy"] == "included"
    assert record["expected_windows"] == 5
    assert record["skipped_n_windows"] == 1
    assert record["sentinel_score_windows"] == 0
    assert record["emitted_hits"] == 1
    assert record["bytes"] > 0
    assert record["elapsed_seconds"] >= 0
    assert (
        record["skipped_n_windows"]
        + record["sentinel_score_windows"]
        + record["below_minimum_score_windows"]
        + record["pwm_filtered_windows"]
        + record["emitted_hits"]
    ) == record["expected_windows"]
PY

contract_sql="$temporary/motif-hit-contract.sql"
awk '
    /^CREATE OR REPLACE VIEW motif_metadata AS/ { emit = 1 }
    /^-- Immutable run and file inventory/ { exit }
    emit { print }
' "$repository_root/sql/schema.sql" >"$contract_sql"
(
    cd "$included"
    duckdb contract.duckdb -bail -f "$contract_sql" >/dev/null
    duckdb contract.duckdb -bail -c "
        SELECT CASE WHEN (SELECT count(*) FROM motif_hit) <> 2
            THEN error('logical motif_hit view has the wrong row count') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM motif_hit
            WHERE motif_name <> 'A_RICH'
               OR task_id <> 'synthetic-threshold-0'
               OR strand NOT IN ('+', '-')
               OR \"end\" <> start + 1
        ) THEN error('logical motif_hit view is inconsistent') END;
    " >/dev/null
)

if "$scanner" "${common[@]}" --threshold 0 --show-sequence \
    --outdir "$included_task" >"$temporary/retry.out" 2>&1; then
    echo "E: Sparse writer replaced an existing result." >&2
    exit 1
fi
grep -Fq 'Refusing to replace existing sparse Parquet file' "$temporary/retry.out" || {
    echo "E: Sparse writer did not explain its overwrite refusal." >&2
    exit 1
}

if "$scanner" "${common[@]}" --outdir "$temporary/no-threshold" \
    >"$temporary/no-threshold.out" 2>&1; then
    echo "E: Sparse writer accepted a run without --threshold." >&2
    exit 1
fi
grep -Fq 'requires an explicit --threshold' "$temporary/no-threshold.out"

if "$scanner" --sparse-parquet --genome "$fixture_dir/genome.fna" \
    --fasta-index "$fasta_index" \
    --pssm "$fixture_dir/JASPAR2026_synthetic.jaspar" \
    --motif MA9999.1 --chr 1 --strand + --coordinate-mode bed --threshold 0 \
    --outdir "$temporary/no-identities" >"$temporary/no-identities.out" 2>&1; then
    echo "E: Sparse writer accepted missing explicit identities." >&2
    exit 1
fi
grep -Fq 'requires explicit --motif-set-id and --genome-id' \
    "$temporary/no-identities.out"

echo "Direct sparse-Parquet tests passed."
