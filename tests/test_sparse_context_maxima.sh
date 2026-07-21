#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-maxima.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping sparse context-maximum test." >&2
    exit 0
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * FROM (VALUES
        ('1', 100::BIGINT, 110::BIGINT),
        ('1', 500::BIGINT, 510::BIGINT)
    ) AS t(chrom, anchor_start, anchor_end)
) TO '$temporary/anchors.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        (260::BIGINT, 270::BIGINT, 5.0::FLOAT),
        (661::BIGINT, 671::BIGINT, 9.0::FLOAT)
    ) AS t(start, "end", score)
) TO '$temporary/m1-plus.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        (50::BIGINT, 60::BIGINT, 6.0::FLOAT),
        (340::BIGINT, 350::BIGINT, 4.0::FLOAT)
    ) AS t(start, "end", score)
) TO '$temporary/m1-minus.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        (260::BIGINT, 300::BIGINT, 7.0::FLOAT)
    ) AS t(start, "end", score)
) TO '$temporary/m2-plus.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        (650::BIGINT, 660::BIGINT, 8.0::FLOAT)
    ) AS t(start, "end", score)
) TO '$temporary/m2-minus.parquet' (FORMAT PARQUET);
SQL

"$repository_root/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$temporary/anchors.parquet" \
    --cofactor M1 "$temporary/m1-plus.parquet" "$temporary/m1-minus.parquet" \
    --cofactor M2 "$temporary/m2-plus.parquet" "$temporary/m2-minus.parquet" \
    --output "$temporary/maxima.parquet" \
    --flank 150 --source-score-floor -1 \
    --duckdb "$duckdb" --threads 1 --memory-limit 256MB \
    --temp-directory "$temporary"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW maxima AS SELECT * FROM read_parquet('$temporary/maxima.parquet');

SELECT CASE WHEN (SELECT COUNT(*) FROM maxima) <> 4
    THEN error('sparse context maxima are not rectangular') END;
SELECT CASE WHEN (SELECT COUNT(DISTINCT (chrom, anchor_start, motif_id))
                  FROM maxima) <> 4
    THEN error('sparse context maxima keys are not unique') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM maxima
    WHERE anchor_start = 100 AND motif_id = 'M1' AND context_score = 6
) THEN error('highest M1 score was not selected') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM maxima
    WHERE anchor_start = 500 AND motif_id = 'M1' AND context_score = 4
) THEN error('a distance-151 M1 hit was admitted') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM maxima
    WHERE anchor_start = 100 AND motif_id = 'M2' AND context_score = 7
) THEN error('wide interval-near M2 hit was lost by the center prefilter') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM maxima
    WHERE anchor_start = 500 AND motif_id = 'M2' AND context_score = 8
) THEN error('second M2 anchor maximum is incorrect') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM maxima
    WHERE capture_prefilter_center_bp <> 200
       OR observed_max_anchor_span_bp <> 10
       OR observed_max_context_span_bp <> 40
       OR context_distance_metric <> 'signed_interval_edge_distance'
) THEN error('derived interval-prefilter provenance is incorrect') END;
SQL

[[ -s $temporary/maxima.parquet.run_config.tsv ]] || {
    echo "E: Sparse context maxima lack their source-file run config." >&2
    exit 1
}
grep -Fq $'M1\t' "$temporary/maxima.parquet.run_config.tsv" || {
    echo "E: Sparse context run config lacks the motif accession." >&2
    exit 1
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * FROM (VALUES
        ('1', 100::BIGINT, 110::BIGINT),
        ('2', 500::BIGINT, 510::BIGINT)
    ) AS t(chrom, anchor_start, anchor_end)
) TO '$temporary/multichrom-anchors.parquet' (FORMAT PARQUET);
SQL

if "$repository_root/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$temporary/multichrom-anchors.parquet" \
    --cofactor M1 "$temporary/m1-plus.parquet" "$temporary/m1-minus.parquet" \
    --output "$temporary/multichrom.parquet" --duckdb "$duckdb" \
    >"$temporary/multichrom.out" 2>"$temporary/multichrom.err"
then
    echo "E: Multi-chromosome anchors were accepted for chromosome-free hits." >&2
    exit 1
fi
grep -Fq 'exactly one chromosome' "$temporary/multichrom.err" || {
    echo "E: Multi-chromosome failure was not explained." >&2
    exit 1
}

if "$repository_root/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$temporary/anchors.parquet" \
    --cofactor M1 "$temporary/m1-plus.parquet" "$temporary/m1-minus.parquet" \
    --output "$temporary/maxima.parquet" --duckdb "$duckdb" \
    >"$temporary/reuse.out" 2>"$temporary/reuse.err"
then
    echo "E: Existing output was overwritten." >&2
    exit 1
fi
grep -Fq 'output already exists' "$temporary/reuse.err" || {
    echo "E: Existing-output failure was not explained." >&2
    exit 1
}

echo "Sparse context-maximum tests passed."
