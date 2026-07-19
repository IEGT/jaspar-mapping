#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
genome="$repository_root/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"
fasta_index=""
pssm="$repository_root/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt"
scanner="$repository_root/pssm_scan_parquet"
output=""
range_start=3600000
range_end=3700000
full_chromosome=0
range_was_set=0
readonly motif_set_id=jaspar2026_core_nonredundant
readonly genome_id=homo_sapiens_grch38_ensembl113_primary

usage() {
    cat <<'EOF'
Usage: run_chr1_patz1_tp73_dry_run.sh [OPTIONS]

Options:
  --output DIR       Package directory (derived from the range by default)
  --from START       First chr1 alignment start (default: 3600000)
  --to END           Exclusive end of the scanned DNA interval (default: 3700000)
  --full-chr1        Scan all of chromosome 1; must be requested explicitly
  --genome FASTA     Uncompressed genome FASTA
  --fasta-index FAI  FASTA index (default: FASTA.fai)
  --pssm FILE        JASPAR-format motif file
  --scanner FILE     Arrow-enabled pssm_scan binary
  -h, --help         Show this help

The run includes PATZ1/MA1961.2 and TP73/MA0861.2, both strands,
log2_relative_risk and log_odds, pseudocount 1, and skip-N handling.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output) output=$2; shift 2 ;;
        --from) range_start=$2; range_was_set=1; shift 2 ;;
        --to) range_end=$2; range_was_set=1; shift 2 ;;
        --full-chr1) full_chromosome=1; shift ;;
        --genome) genome=$2; shift 2 ;;
        --fasta-index) fasta_index=$2; shift 2 ;;
        --pssm) pssm=$2; shift 2 ;;
        --scanner) scanner=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ $full_chromosome -eq 0 || $range_was_set -eq 0 ]] || {
    echo "E: --full-chr1 cannot be combined with --from or --to." >&2
    exit 2
}
[[ $range_start =~ ^[0-9]+$ && $range_end =~ ^[0-9]+$ ]] || {
    echo "E: --from and --to must be non-negative integers." >&2
    exit 2
}

command -v python3 >/dev/null 2>&1 || { echo "E: python3 is required." >&2; exit 1; }
command -v duckdb >/dev/null 2>&1 || { echo "E: duckdb CLI is required." >&2; exit 1; }

resolve_repo_path() {
    if [[ $1 = /* ]]; then
        printf '%s\n' "$1"
    else
        printf '%s/%s\n' "$repository_root" "$1"
    fi
}

genome=$(resolve_repo_path "$genome")
pssm=$(resolve_repo_path "$pssm")
scanner=$(resolve_repo_path "$scanner")
[[ -n $fasta_index ]] && fasta_index=$(resolve_repo_path "$fasta_index")
fasta_index=${fasta_index:-$genome.fai}

[[ -f $genome ]] || { echo "E: Genome FASTA not found: $genome" >&2; exit 1; }
[[ -f $pssm ]] || { echo "E: JASPAR file not found: $pssm" >&2; exit 1; }

if [[ ! -x $scanner ]]; then
    echo "I: Building Arrow-enabled scanner." >&2
    make -C "$repository_root" pssm_scan_parquet
fi

if [[ ! -f $fasta_index ]]; then
    python3 "$repository_root/scripts/build_fasta_index.py" \
        "$genome" --output "$fasta_index"
fi

chromosome_length=$(awk -F '\t' '$1 == "1" { print $2; exit }' "$fasta_index")
[[ $chromosome_length =~ ^[0-9]+$ ]] || {
    echo "E: Chromosome 1 is absent from FASTA index: $fasta_index" >&2
    exit 1
}

range_arguments=()
part_range="from=$range_start-to=$range_end"
if [[ $full_chromosome -eq 1 ]]; then
    range_start=0
    range_end=$chromosome_length
    part_range="from=0-to=end"
else
    (( range_end > range_start )) || { echo "E: --to must be greater than --from." >&2; exit 2; }
    (( range_end <= chromosome_length )) || { echo "E: --to exceeds chromosome 1 length." >&2; exit 2; }
    range_arguments=(--from "$range_start" --to "$range_end")
fi
(( range_end - range_start >= 16 )) || {
    echo "E: Requested interval is shorter than the TP73 motif." >&2
    exit 2
}

if [[ -z $output ]]; then
    if [[ $full_chromosome -eq 1 ]]; then
        output="$repository_root/dry_runs/chr1_patz1_tp73_full"
    else
        output="$repository_root/dry_runs/chr1_patz1_tp73_from-${range_start}-to-${range_end}"
    fi
else
    output=$(resolve_repo_path "$output")
fi
mkdir -p "$output"
output=$(cd "$output" && pwd -P)

if [[ $full_chromosome -eq 1 ]]; then
    available_kib=$(df -Pk "$output" | awk 'NR == 2 { print $4 }')
    required_kib=$((12 * 1024 * 1024))
    (( available_kib >= required_kib )) || {
        echo "E: Full chr1 mode requires at least 12 GiB free in the output filesystem." >&2
        exit 1
    }
fi

for motif in MA1961.2 MA0861.2; do
    grep -q "^>${motif}[[:space:]]" "$pssm" || {
        echo "E: Motif $motif is absent from $pssm" >&2
        exit 1
    }
done

part_file="part-${part_range}-n_policy=skip-000000.parquet"
for motif in MA1961.2 MA0861.2; do
    for mode in log2_relative_risk log_odds; do
        partition="$output/tables/jaspar2026/motif_score_dense/motif_set_id=$motif_set_id/genome_id=$genome_id/motif_id=$motif/score_mode=$mode/pseudocount=1/background_model_id=uniform_acgt_v1/pseudocount_scheme=additive_per_base/chrom=1"
        plus_file="$partition/strand=plus/$part_file"
        minus_file="$partition/strand=minus/$part_file"
        if [[ -s $plus_file && -s $minus_file ]]; then
            echo "I: Reusing complete $motif $mode Parquet pair." >&2
            continue
        fi
        "$scanner" \
            --dense-scores \
            --dense-block-size 65536 \
            --outdir "$output" \
            --genome "$genome" \
            --fasta-index "$fasta_index" \
            --pssm "$pssm" \
            --motif "$motif" \
            --motif-set-id "$motif_set_id" \
            --genome-id "$genome_id" \
            --chr 1 \
            --strand both \
            --coordinate-mode bed \
            --score-mode "$mode" \
            --pseudocount 1 \
            --skip-N \
            "${range_arguments[@]}"
    done
done

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{ print $1 }'
    else
        shasum -a 256 "$1" | awk '{ print $1 }'
    fi
}

jaspar_sha256=$(sha256_file "$pssm")
scanner_sha256=$(sha256_file "$scanner")
fasta_index_sha256=$(sha256_file "$fasta_index")

metadata_dir="$output/tables/jaspar2026/motif_metadata"
metadata_file="$metadata_dir/part-000000.parquet"
mkdir -p "$metadata_dir"
if [[ ! -f $metadata_file ]]; then
    (
        cd "$output"
        duckdb :memory: -bail -c "COPY (
            SELECT * FROM (VALUES
                ('$motif_set_id', 'MA1961.2', 'PATZ1', 11, 2026, '$jaspar_sha256'),
                ('$motif_set_id', 'MA0861.2', 'TP73', 16, 2026, '$jaspar_sha256')
            ) AS t(motif_set_id, motif_id, motif_name, motif_length,
                   jaspar_version, source_sha256)
        ) TO 'tables/jaspar2026/motif_metadata/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);" >/dev/null
    )
fi

database_name=jaspar2026_chr1_patz1_tp73.duckdb
schema="$repository_root/sql/chr1_dense_dry_run_schema.sql"
source_commit=$(git -C "$repository_root" rev-parse HEAD)
source_tree_dirty=false
if ! git -C "$repository_root" diff --quiet --ignore-submodules -- ||
   ! git -C "$repository_root" diff --cached --quiet --ignore-submodules --; then
    source_tree_dirty=true
fi
generated_at=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
duckdb_version=$(duckdb --version)
run_id="chr1_patz1_tp73_${range_start}_${range_end}"
full_literal=false
[[ $full_chromosome -eq 1 ]] && full_literal=true

sql_escape() {
    printf '%s' "$1" | sed "s/'/''/g"
}

manifest_sql="CREATE OR REPLACE TABLE run_manifest AS SELECT
    '$(sql_escape "$run_id")'::VARCHAR AS run_id,
    '$(sql_escape "$generated_at")'::VARCHAR AS generated_at_utc,
    '$(sql_escape "$source_commit")'::VARCHAR AS source_commit,
    $source_tree_dirty::BOOLEAN AS source_tree_dirty,
    '$(sql_escape "$(basename "$scanner")")'::VARCHAR AS scanner_file,
    '$(sql_escape "$scanner_sha256")'::VARCHAR AS scanner_sha256,
    '$(sql_escape "$duckdb_version")'::VARCHAR AS duckdb_version,
    2026::INTEGER AS jaspar_version,
    '$(sql_escape "$jaspar_sha256")'::VARCHAR AS jaspar_sha256,
    '$(sql_escape "$(basename "$pssm")")'::VARCHAR AS pssm_file,
    '$(sql_escape "$(basename "$genome")")'::VARCHAR AS genome_file,
    '$(sql_escape "$(basename "$fasta_index")")'::VARCHAR AS fasta_index_file,
    '$(sql_escape "$fasta_index_sha256")'::VARCHAR AS fasta_index_sha256,
    '1'::VARCHAR AS chrom,
    $chromosome_length::BIGINT AS chrom_length,
    $range_start::BIGINT AS range_start,
    $range_end::BIGINT AS range_end,
    $full_literal::BOOLEAN AS full_chromosome,
    1.0::DOUBLE AS pseudocount,
    '+,-'::VARCHAR AS strands,
    'log2_relative_risk,log_odds'::VARCHAR AS score_modes,
    'MA1961.2,MA0861.2'::VARCHAR AS motif_ids,
    'bed'::VARCHAR AS coordinate_mode,
    'skip'::VARCHAR AS n_policy;"

validation_sql="SELECT CASE
    WHEN (SELECT COUNT(*) FROM dense_run_inventory) <> 8
    THEN error('expected exactly eight motif/mode/strand configurations')
END;
SELECT CASE
    WHEN EXISTS (
        SELECT 1
        FROM dense_run_inventory i
        JOIN motif_metadata m USING (motif_set_id, motif_id)
        CROSS JOIN run_manifest r
        WHERE i.chrom <> '1'
           OR i.pseudocount <> 1.0
           OR i.part_files <> 1
           OR i.alignment_start_begin <> r.range_start
           OR i.alignment_start_end <> r.range_end - m.motif_length + 1
           OR i.n_windows <> r.range_end - r.range_start - m.motif_length + 1
    )
    THEN error('dense Parquet inventory does not match the declared run interval')
END;
SELECT CASE
    WHEN (SELECT COUNT(*) FROM motif_metadata) <> 2
      OR EXISTS (SELECT 1 FROM motif_metadata WHERE source_sha256 <> '$jaspar_sha256')
    THEN error('motif metadata does not match the selected JASPAR file')
END;"

(
    cd "$output"
    duckdb -bail "$database_name" -f "$schema" >/dev/null
    duckdb -bail "$database_name" -c "$manifest_sql" >/dev/null
    duckdb -bail "$database_name" -c "$validation_sql" >/dev/null
)

manifest_tmp=$(mktemp "$output/.manifest.json.XXXXXX")
trap 'rm -f "$manifest_tmp"' EXIT HUP INT TERM
(
    cd "$output"
    duckdb -readonly -json "$database_name" -c "SELECT * FROM run_manifest;" > "$manifest_tmp"
)
mv "$manifest_tmp" "$output/manifest.json"
trap - EXIT HUP INT TERM

echo "I: Dense chr1 dry-run package is ready: $output" >&2
bash "$repository_root/scripts/inspect_chr1_dense_dry_run.sh" \
    --package "$output" overview
