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

common=(
    --sparse-parquet
    --genome "$fixture_dir/genome.fna"
    --fasta-index "$fasta_index"
    --pssm "$fixture_dir/JASPAR2026_synthetic.jaspar"
    --motif MA9999.1
    --chr 1
    --strand both
    --coordinate-mode bed
    --score-mode log2_relative_risk
    --pseudocount 0
    --skip-N
)

included="$temporary/nested/output/included"
omitted="$temporary/omitted"
"$scanner" "${common[@]}" --threshold 0 --show-sequence \
    --outdir "$included" >/dev/null
"$scanner" "${common[@]}" --threshold -1 \
    --outdir "$omitted" >/dev/null

included_glob="$included/tables/jaspar2026/motif_hit/**/*.parquet"
omitted_glob="$omitted/tables/jaspar2026/motif_hit/**/*.parquet"
included_plus=$(find "$included" -path '*/strand=plus/*.parquet' -print -quit)
[[ -n $included_plus ]] || {
    echo "E: Sparse writer did not create the plus-strand partition." >&2
    exit 1
}

metadata_dir="$included/tables/jaspar2026/motif_metadata"
mkdir -p "$metadata_dir"
duckdb :memory: -bail -c "COPY (
    SELECT 'MA9999.1'::VARCHAR AS motif_id,
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
           OR score_mode <> 'log2_relative_risk'
           OR pseudocount <> 0
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
        WHERE name IN ('motif_id', 'score_mode', 'pseudocount', 'chrom', 'strand')
    ) THEN error('row-constant identity leaked into the physical schema') END;
    SELECT CASE WHEN (
        SELECT count(*) FROM parquet_schema('$included_plus')
        WHERE name IN ('start', 'end', 'score', 'pwm_relative_score', 'matched_seq')
    ) <> 5 THEN error('physical sparse schema is incomplete') END;
" >/dev/null

contract_sql="$temporary/motif-hit-contract.sql"
awk '
    /^CREATE OR REPLACE VIEW motif_metadata AS/ { emit = 1 }
    /^-- Dense score blocks/ { exit }
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
               OR strand NOT IN ('+', '-')
               OR \"end\" <> start + 1
        ) THEN error('logical motif_hit view is inconsistent') END;
    " >/dev/null
)

if "$scanner" "${common[@]}" --threshold 0 --show-sequence \
    --outdir "$included" >"$temporary/retry.out" 2>&1; then
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

echo "Direct sparse-Parquet tests passed."
