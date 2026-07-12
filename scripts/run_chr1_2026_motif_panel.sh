#!/usr/bin/env bash

set -euo pipefail

readonly JASPAR_URL="https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt"
readonly JASPAR_SHA256="0dc9b7f9e159a8376c2e52edf373863ae21518bb760ceeea494ad6b746092e53"
readonly GENOME_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
readonly EXPECTED_CHR1_LENGTH=248956422

readonly -a MOTIF_IDS=(MA0861.2 MA0024.3 MA0079.5 MA1961.2 MA0507.3)
readonly -a MOTIF_NAMES=(TP73 E2F1 SP1 PATZ1 POU2F2)
readonly -a MOTIF_LENGTHS=(16 12 9 11 13)
readonly -a SCORE_MODES=(log2_relative_risk log_odds)
readonly -a STRANDS=(+ -)

usage() {
    cat <<'EOF'
Usage: run_chr1_2026_motif_panel.sh --run-root DIR [OPTIONS]

Run a dense JASPAR 2026 chromosome 1 panel inside an existing Slurm allocation.
The panel contains TP73/MA0861.2, E2F1/MA0024.3, SP1/MA0079.5,
PATZ1/MA1961.2, and POU2F2/MA0507.3 in both score modes and orientations:
20 independent scans in total.

Options:
  --run-root DIR    Dedicated directory for every generated file (required)
  --source DIR      Public jaspar-mapping Git checkout (default: script parent)
  --micromamba FILE Micromamba executable (default: $MAMBA_EXE or micromamba)
  --print-config    Print the pinned panel and public inputs, then exit
  -h, --help        Show this help and exit

The job downloads version-pinned public JASPAR and Ensembl inputs, creates an
isolated Conda environment below DIR, builds the Arrow-enabled scanner from the
checked-out Git commit, and writes validated Parquet through per-task staging.
No existing file is replaced or removed. A rerun reuses only outputs that pass
structural and window-count validation.

Recommended Haumea submission:
  sbatch --account=cluster --partition=requeue --ntasks=20 --mem=20G \
    --time=1-00:00:00 --requeue --chdir=DIR \
    --output=DIR/logs/slurm-%j.out --error=DIR/logs/slurm-%j.err \
    scripts/run_chr1_2026_motif_panel.sh --run-root DIR
EOF
}

print_config() {
    printf 'Key\tValue\n'
    printf 'JASPAR_URL\t%s\n' "$JASPAR_URL"
    printf 'JASPAR_SHA256\t%s\n' "$JASPAR_SHA256"
    printf 'GENOME_URL\t%s\n' "$GENOME_URL"
    printf 'Chromosome\t1\n'
    printf 'ChromosomeLength\t%s\n' "$EXPECTED_CHR1_LENGTH"
    printf 'Pseudocount\t1\n'
    printf 'ScoreModes\t%s\n' "${SCORE_MODES[*]}"
    printf 'Strands\t%s\n' "${STRANDS[*]}"
    for index in "${!MOTIF_IDS[@]}"; do
        printf 'Motif\t%s:%s:%s\n' \
            "${MOTIF_NAMES[$index]}" "${MOTIF_IDS[$index]}" \
            "${MOTIF_LENGTHS[$index]}"
    done
}

run_root=""
source_dir=$(cd "$(dirname "$0")/.." && pwd)
micromamba_executable=${MAMBA_EXE:-micromamba}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root)
            [[ $# -ge 2 ]] || { echo "E: --run-root requires a directory." >&2; exit 2; }
            run_root=$2
            shift 2
            ;;
        --source)
            [[ $# -ge 2 ]] || { echo "E: --source requires a directory." >&2; exit 2; }
            source_dir=$2
            shift 2
            ;;
        --micromamba)
            [[ $# -ge 2 ]] || { echo "E: --micromamba requires a file." >&2; exit 2; }
            micromamba_executable=$2
            shift 2
            ;;
        --print-config)
            print_config
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "E: Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -n $run_root ]] || { echo "E: --run-root is required." >&2; exit 2; }
[[ -n ${SLURM_JOB_ID:-} ]] || {
    echo "E: Run this script inside a Slurm allocation or submit it with sbatch." >&2
    exit 1
}
[[ ${SLURM_NTASKS:-0} =~ ^[0-9]+$ && ${SLURM_NTASKS:-0} -ge 20 ]] || {
    echo "E: The panel requires at least 20 Slurm tasks; found ${SLURM_NTASKS:-0}." >&2
    exit 1
}
[[ -d $source_dir/.git ]] || { echo "E: Git checkout not found: $source_dir" >&2; exit 1; }

source_dir=$(cd "$source_dir" && pwd -P)
mkdir -p "$run_root"
run_root=$(cd "$run_root" && pwd -P)

readonly input_dir="$run_root/input/public"
readonly output_dir="$run_root/output"
readonly logs_dir="$run_root/logs"
readonly staging_root="$run_root/staging/job-${SLURM_JOB_ID}-restart-${SLURM_RESTART_COUNT:-0}"
readonly environment_dir="$run_root/env"
readonly mamba_root="$run_root/mamba-root"
readonly package_dir="$output_dir/jaspar2026_chr1_dense_5motifs_full"
mkdir -p "$input_dir" "$output_dir" "$logs_dir" "$staging_root" "$mamba_root"

source_commit=$(git -C "$source_dir" rev-parse HEAD)
git -C "$source_dir" diff --quiet --ignore-submodules -- || {
    echo "E: Tracked source changes are present in $source_dir." >&2
    exit 1
}
git -C "$source_dir" diff --cached --quiet --ignore-submodules -- || {
    echo "E: Staged source changes are present in $source_dir." >&2
    exit 1
}

readonly code_dir="$run_root/code-$source_commit"
if [[ ! -d $code_dir ]]; then
    source_stage="$staging_root/source-$source_commit"
    [[ ! -e $source_stage ]] || {
        echo "E: Source staging path already exists: $source_stage" >&2
        exit 1
    }
    mkdir -p "$source_stage"
    git -C "$source_dir" archive "$source_commit" | tar -x -C "$source_stage"
    printf '%s\n' "$source_commit" > "$source_stage/SOURCE_COMMIT"
    [[ ! -e $code_dir ]] || { echo "E: Code directory appeared concurrently." >&2; exit 1; }
    mv "$source_stage" "$code_dir"
fi
[[ $(<"$code_dir/SOURCE_COMMIT") == "$source_commit" ]] || {
    echo "E: Code snapshot does not match source commit $source_commit." >&2
    exit 1
}

if [[ $micromamba_executable != */* ]]; then
    micromamba_executable=$(command -v "$micromamba_executable" || true)
fi
[[ -x $micromamba_executable ]] || {
    echo "E: Micromamba executable not found: ${micromamba_executable:-unset}" >&2
    exit 1
}

if [[ ! -x $environment_dir/bin/pkg-config ]]; then
    "$micromamba_executable" --root-prefix "$mamba_root" create \
        --prefix "$environment_dir" --override-channels --channel conda-forge --yes \
        libarrow=25.0.0 libparquet=25.0.0 duckdb-cli=1.5.4 \
        pkg-config=0.29.2 zlib=1.3.2 bzip2=1.0.8
fi

export PATH="$environment_dir/bin:$PATH"
export PKG_CONFIG_PATH="$environment_dir/lib/pkgconfig"
export LD_LIBRARY_PATH="$environment_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

[[ $(pkg-config --modversion arrow) == 25.0.0 ]] || { echo "E: Arrow 25.0.0 is required." >&2; exit 1; }
[[ $(pkg-config --modversion parquet) == 25.0.0 ]] || { echo "E: Parquet 25.0.0 is required." >&2; exit 1; }
duckdb --version | grep -q 'v1\.5\.4' || { echo "E: DuckDB 1.5.4 is required." >&2; exit 1; }

record_url() {
    local path=$1
    local url=$2
    local url_file="$path.url"
    if [[ -e $url_file ]]; then
        [[ $(<"$url_file") == "$url" ]] || {
            echo "E: Recorded URL differs for $path." >&2
            return 1
        }
    else
        printf '%s\n' "$url" > "$url_file"
    fi
}

download_once() {
    local url=$1
    local destination=$2
    local partial="$destination.download-${SLURM_JOB_ID}-${SLURM_RESTART_COUNT:-0}"
    record_url "$destination" "$url"
    if [[ -s $destination ]]; then
        echo "I: Reusing public input $destination"
        return
    fi
    [[ ! -e $destination ]] || { echo "E: Empty input exists: $destination" >&2; return 1; }
    [[ ! -e $partial ]] || { echo "E: Download staging file exists: $partial" >&2; return 1; }
    curl --fail --location --retry 5 --output "$partial" "$url"
    [[ -s $partial ]] || { echo "E: Download is empty: $url" >&2; return 1; }
    [[ ! -e $destination ]] || { echo "E: Input appeared concurrently: $destination" >&2; return 1; }
    mv "$partial" "$destination"
}

readonly jaspar_file="$input_dir/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt"
readonly genome_gz="$input_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
readonly genome_file="$input_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"
readonly fasta_index="$genome_file.fai"

download_once "$JASPAR_URL" "$jaspar_file"
[[ $(sha256sum "$jaspar_file" | awk '{print $1}') == "$JASPAR_SHA256" ]] || {
    echo "E: JASPAR checksum mismatch." >&2
    exit 1
}

download_once "$GENOME_URL" "$genome_gz"
gzip -t "$genome_gz"

if [[ ! -s $genome_file ]]; then
    genome_partial="$genome_file.decompress-${SLURM_JOB_ID}-${SLURM_RESTART_COUNT:-0}"
    [[ ! -e $genome_file ]] || { echo "E: Empty genome FASTA exists: $genome_file" >&2; exit 1; }
    [[ ! -e $genome_partial ]] || { echo "E: Genome staging file exists: $genome_partial" >&2; exit 1; }
    gzip -dc "$genome_gz" > "$genome_partial"
    [[ -s $genome_partial ]] || { echo "E: Decompressed genome is empty." >&2; exit 1; }
    mv "$genome_partial" "$genome_file"
fi

if [[ ! -s $fasta_index ]]; then
    [[ ! -e $fasta_index ]] || { echo "E: Empty FASTA index exists: $fasta_index" >&2; exit 1; }
    python3 "$code_dir/scripts/build_fasta_index.py" "$genome_file" --output "$fasta_index"
fi
chromosome_length=$(awk -F '\t' '$1 == "1" { print $2; exit }' "$fasta_index")
[[ $chromosome_length == "$EXPECTED_CHR1_LENGTH" ]] || {
    echo "E: Expected chr1 length $EXPECTED_CHR1_LENGTH, found ${chromosome_length:-missing}." >&2
    exit 1
}

for motif in "${MOTIF_IDS[@]}"; do
    grep -q "^>${motif}[[:space:]]" "$jaspar_file" || {
        echo "E: Motif $motif is absent from the public JASPAR file." >&2
        exit 1
    }
done

make -C "$code_dir" pssm_scan_parquet
readonly scanner="$code_dir/pssm_scan_parquet"
[[ -x $scanner ]] || { echo "E: Scanner build failed." >&2; exit 1; }

sql_quote() {
    printf '%s' "$1" | sed "s/'/''/g"
}

validate_dense_file() {
    local file=$1
    local expected_windows=$2
    [[ -s $file ]] || return 1
    local quoted_file
    quoted_file=$(sql_quote "$file")
    local result
    result=$(duckdb :memory: -noheader -list -c "
        SELECT CASE
            WHEN COUNT(*) > 0
             AND MIN(block_start) = 0
             AND MAX(block_start + len(scores)) = $expected_windows
             AND SUM(len(scores)) = $expected_windows
            THEN 'ok' ELSE 'invalid' END
        FROM read_parquet('$quoted_file');" 2>/dev/null) || return 1
    [[ $result == ok ]]
}

strand_label() {
    [[ $1 == + ]] && printf 'plus\n' || printf 'minus\n'
}

run_one() {
    local motif_id=$1
    local motif_length=$2
    local score_mode=$3
    local strand=$4
    local strand_name
    strand_name=$(strand_label "$strand")
    local label="${motif_id}_${score_mode}_${strand_name}"
    local expected_windows=$((chromosome_length - motif_length + 1))
    local relative_file="tables/jaspar2026/motif_score_dense/motif_id=${motif_id}/score_mode=${score_mode}/pseudocount=1/chrom=1/strand=${strand_name}/part-from=0-to=end-n_policy=skip-000000.parquet"
    local final_file="$package_dir/$relative_file"

    if [[ -e $final_file ]]; then
        validate_dense_file "$final_file" "$expected_windows" || {
            echo "E: Existing final output is incomplete or invalid: $final_file" >&2
            return 1
        }
        echo "I: Reusing validated $label"
        return
    fi

    local task_stage="$staging_root/$label"
    [[ ! -e $task_stage ]] || { echo "E: Task staging path exists: $task_stage" >&2; return 1; }
    mkdir -p "$task_stage"
    echo "I: Starting $label"
    srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=none \
        "$scanner" -v --dense-scores --dense-block-size 65536 \
        --outdir "$task_stage" \
        --genome "$genome_file" --fasta-index "$fasta_index" \
        --pssm "$jaspar_file" --motif "$motif_id" --chr 1 \
        --strand "$strand" --coordinate-mode bed \
        --score-mode "$score_mode" --pseudocount 1 --skip-N \
        > "$logs_dir/${SLURM_JOB_ID}_${label}.out" \
        2> "$logs_dir/${SLURM_JOB_ID}_${label}.err"

    local staged_file="$task_stage/$relative_file"
    validate_dense_file "$staged_file" "$expected_windows" || {
        echo "E: Newly generated output failed validation: $staged_file" >&2
        return 1
    }
    mkdir -p "$(dirname "$final_file")"
    [[ ! -e $final_file ]] || { echo "E: Final file appeared concurrently: $final_file" >&2; return 1; }
    mv "$staged_file" "$final_file"
    echo "I: Completed and promoted $label"
}

mkdir -p "$package_dir"
pids=()
labels=()
for index in "${!MOTIF_IDS[@]}"; do
    for score_mode in "${SCORE_MODES[@]}"; do
        for strand in "${STRANDS[@]}"; do
            label="${MOTIF_IDS[$index]}_${score_mode}_$(strand_label "$strand")"
            run_one "${MOTIF_IDS[$index]}" "${MOTIF_LENGTHS[$index]}" \
                "$score_mode" "$strand" &
            pids+=("$!")
            labels+=("$label")
        done
    done
done

failed=0
for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
        echo "E: Scan failed: ${labels[$index]}" >&2
        failed=1
    fi
done
[[ $failed -eq 0 ]] || exit 1

metadata_dir="$package_dir/tables/jaspar2026/motif_metadata"
metadata_file="$metadata_dir/part-000000.parquet"
if [[ ! -e $metadata_file ]]; then
    metadata_stage="$staging_root/motif_metadata.parquet"
    [[ ! -e $metadata_stage ]] || { echo "E: Metadata staging file exists." >&2; exit 1; }
    duckdb :memory: -bail -c "COPY (
        SELECT * FROM (VALUES
            ('MA0861.2', 'TP73', 16, 2026, '$JASPAR_SHA256'),
            ('MA0024.3', 'E2F1', 12, 2026, '$JASPAR_SHA256'),
            ('MA0079.5', 'SP1', 9, 2026, '$JASPAR_SHA256'),
            ('MA1961.2', 'PATZ1', 11, 2026, '$JASPAR_SHA256'),
            ('MA0507.3', 'POU2F2', 13, 2026, '$JASPAR_SHA256')
        ) AS t(motif_id, motif_name, motif_length, jaspar_version, source_sha256)
    ) TO '$(sql_quote "$metadata_stage")' (FORMAT PARQUET, COMPRESSION ZSTD);"
    mkdir -p "$metadata_dir"
    [[ ! -e $metadata_file ]] || { echo "E: Metadata appeared concurrently." >&2; exit 1; }
    mv "$metadata_stage" "$metadata_file"
fi

scanner_sha256=$(sha256sum "$scanner" | awk '{print $1}')
genome_gz_sha256=$(sha256sum "$genome_gz" | awk '{print $1}')
fasta_index_sha256=$(sha256sum "$fasta_index" | awk '{print $1}')
generated_at=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
database_name=jaspar2026_chr1_dense_5motifs.duckdb
database_file="$package_dir/$database_name"

if [[ ! -e $database_file ]]; then
    database_stage="$staging_root/$database_name"
    (
        cd "$package_dir"
        duckdb "$database_stage" -bail -f "$code_dir/sql/chr1_dense_dry_run_schema.sql" >/dev/null
        duckdb "$database_stage" -bail -c "
            CREATE OR REPLACE TABLE run_manifest AS SELECT
                'jaspar2026_chr1_dense_5motifs_full'::VARCHAR AS run_id,
                '$generated_at'::VARCHAR AS generated_at_utc,
                '$source_commit'::VARCHAR AS source_commit,
                '$JASPAR_URL'::VARCHAR AS jaspar_url,
                '$JASPAR_SHA256'::VARCHAR AS jaspar_sha256,
                '$GENOME_URL'::VARCHAR AS genome_url,
                '$genome_gz_sha256'::VARCHAR AS genome_gz_sha256,
                '$fasta_index_sha256'::VARCHAR AS fasta_index_sha256,
                '$scanner_sha256'::VARCHAR AS scanner_sha256,
                '1'::VARCHAR AS chrom,
                $chromosome_length::BIGINT AS chrom_length,
                1.0::DOUBLE AS pseudocount,
                '+,-'::VARCHAR AS strands,
                'log2_relative_risk,log_odds'::VARCHAR AS score_modes,
                'MA0861.2,MA0024.3,MA0079.5,MA1961.2,MA0507.3'::VARCHAR AS motif_ids,
                'bed'::VARCHAR AS coordinate_mode,
                'skip'::VARCHAR AS n_policy,
                '${SLURM_JOB_ID}'::VARCHAR AS slurm_job_id;

            SELECT CASE WHEN (SELECT COUNT(*) FROM dense_run_inventory) <> 20
                THEN error('expected 20 dense configurations') END;
            SELECT CASE WHEN EXISTS (
                SELECT 1 FROM dense_run_inventory i
                JOIN motif_metadata m USING (motif_id)
                WHERE i.part_files <> 1 OR i.alignment_start_begin <> 0
                   OR i.alignment_start_end <> $chromosome_length - m.motif_length + 1
                   OR i.n_windows <> $chromosome_length - m.motif_length + 1
            ) THEN error('dense inventory failed final validation') END;" >/dev/null
    )
    [[ ! -e $database_file ]] || { echo "E: Database appeared concurrently." >&2; exit 1; }
    mv "$database_stage" "$database_file"
fi

manifest_file="$package_dir/manifest.json"
if [[ ! -e $manifest_file ]]; then
    manifest_stage="$staging_root/manifest.json"
    cat > "$manifest_stage" <<EOF
{
  "run_id": "jaspar2026_chr1_dense_5motifs_full",
  "generated_at_utc": "$generated_at",
  "source_commit": "$source_commit",
  "jaspar_url": "$JASPAR_URL",
  "jaspar_sha256": "$JASPAR_SHA256",
  "genome_url": "$GENOME_URL",
  "genome_gz_sha256": "$genome_gz_sha256",
  "fasta_index_sha256": "$fasta_index_sha256",
  "scanner_sha256": "$scanner_sha256",
  "chrom": "1",
  "chrom_length": $chromosome_length,
  "motif_ids": ["MA0861.2", "MA0024.3", "MA0079.5", "MA1961.2", "MA0507.3"],
  "score_modes": ["log2_relative_risk", "log_odds"],
  "strands": ["+", "-"],
  "pseudocount": 1,
  "coordinate_mode": "bed",
  "n_policy": "skip",
  "slurm_job_id": "${SLURM_JOB_ID}"
}
EOF
    mv "$manifest_stage" "$manifest_file"
fi

echo "I: Full chromosome 1 motif panel completed: $package_dir"
du -sh "$package_dir"
