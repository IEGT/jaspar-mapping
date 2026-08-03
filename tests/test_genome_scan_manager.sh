#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-genome-scan-manager.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required." >&2
    exit 1
}
scanner="$repository_root/pssm_scan_parquet"
[[ -x $scanner ]] || {
    echo "E: Build pssm_scan_parquet before this test." >&2
    exit 1
}

cat > "$temporary/genome.fna" <<'EOF'
>1 synthetic chromosome
ACGTNACGT
EOF
python3 "$repository_root/scripts/build_fasta_index.py" \
    "$temporary/genome.fna" --output "$temporary/genome.fna.fai" >/dev/null

cat > "$temporary/motifs.jaspar" <<'EOF'
>MA0861.2 TP73_TEST
A [ 4 1 ]
C [ 1 4 ]
G [ 1 1 ]
T [ 1 1 ]
>MA0001.1 OTHER_TEST
A [ 1 ]
C [ 1 ]
G [ 4 ]
T [ 1 ]
EOF

source_commit=0123456789abcdef0123456789abcdef01234567
run_root="$temporary/run"
manager="$repository_root/scripts/manage_genome_scan.py"

task_count=$(
    "$manager" prepare \
        --run-root "$run_root" --run-id synthetic_jaspar_scan_v1 \
        --source-commit "$source_commit" \
        --genome "$temporary/genome.fna" \
        --fasta-index "$temporary/genome.fna.fai" \
        --jaspar "$temporary/motifs.jaspar" \
        --genome-id synthetic_genome_v1 --motif-set-id synthetic_jaspar_v1 \
        --taxon-id 9606 --species 'Homo sapiens synthetic' \
        --assembly-name SYN1 --assembly-accession GCA_SYNTHETIC \
        --ensembl-release 113 \
        --fasta-url https://example.invalid/genome.fna \
        --jaspar-url https://example.invalid/motifs.jaspar \
        --chrom 1 --motif-batch-size 1 --minimum-free-gib 0 \
        | awk -F= '$1 == "TASK_COUNT" { print $2 }'
)
[[ $task_count == 2 ]] || {
    echo "E: Expected two non-overlapping motif tasks, got $task_count." >&2
    exit 1
}

"$manager" run-chromosome --run-root "$run_root" --chrom 1 \
    --scanner "$scanner" --duckdb duckdb --allow-local \
    --scratch-directory "$temporary/chromosome-scratch" \
    --minimum-scratch-free-gib 0 --batch-workers 2 \
    --scratch-task-output --allow-scanner-provenance-mismatch >/dev/null
[[ -s $temporary/chromosome-scratch/1.fa ]] || {
    echo "E: Chromosome worker did not stage its FASTA." >&2
    exit 1
}
[[ -s $temporary/chromosome-scratch/1.fa.fai ]] || {
    echo "E: Chromosome worker did not index its staged FASTA." >&2
    exit 1
}
[[ -s $temporary/chromosome-scratch/motifs.jaspar ]] || {
    echo "E: Chromosome worker did not stage its JASPAR source." >&2
    exit 1
}

# A retry must validate and reuse the promoted task rather than rescan it.
"$manager" run-task --run-root "$run_root" --task-index 0 \
    --scanner "$scanner" --duckdb duckdb --allow-local \
    --allow-scanner-provenance-mismatch >/dev/null

real_duckdb=$(type -P duckdb)
cat > "$temporary/failing-duckdb" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
for argument in "$@"; do
    if [[ $argument == */genome_scan_schema.sql ]]; then
        echo "deliberate catalog-build failure" >&2
        exit 17
    fi
done
exec "$REAL_DUCKDB" "$@"
EOF
chmod +x "$temporary/failing-duckdb"
if REAL_DUCKDB="$real_duckdb" "$manager" finalize --run-root "$run_root" \
    --duckdb "$temporary/failing-duckdb" \
    >"$temporary/finalize-failure.stdout" 2>"$temporary/finalize-failure.stderr"; then
    echo "E: Deliberately interrupted finalizer unexpectedly succeeded." >&2
    exit 1
fi
grep -Fq 'deliberate catalog-build failure' "$temporary/finalize-failure.stderr"
[[ ! -e $run_root/package/manifest.json ]] || {
    echo "E: Failed finalizer published its completion manifest." >&2
    exit 1
}

# Retrying after the failed catalog phase must publish one coherent package,
# with the manifest still appearing last.
"$manager" finalize --run-root "$run_root" --duckdb duckdb >/dev/null
"$manager" finalize --run-root "$run_root" --duckdb duckdb >/dev/null

(
    cd "$run_root/package"
    duckdb jaspar_genome_scan.duckdb -bail -c "
        SELECT CASE WHEN (SELECT COUNT(*) FROM scan_task) <> 2
            THEN error('scan_task inventory is incomplete') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM scan_file_inventory) <> 4
            THEN error('scan_file_inventory is incomplete') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM scan_task_completeness WHERE NOT complete
        ) THEN error('a promoted task is incomplete') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM scan_threshold_policy
            WHERE policy_id = 'tp73_calibrated' AND minimum_score = -5
        ) OR NOT EXISTS (
            SELECT 1 FROM scan_threshold_policy
            WHERE policy_id = 'default_uncalibrated' AND minimum_score = -1
        ) THEN error('threshold policies were not preserved') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM scan_motif_threshold
            WHERE motif_id = 'MA0861.2' AND final_minimum_score = -5
        ) OR NOT EXISTS (
            SELECT 1 FROM scan_motif_threshold
            WHERE motif_id = 'MA0001.1' AND final_minimum_score = -1
        ) THEN error('legacy per-motif thresholds were not materialized') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM scan_file_inventory
            WHERE state <> 'complete'
               OR expected_windows < emitted_hits
               OR expected_windows <> skipped_n_windows
                    + sentinel_score_windows + below_minimum_score_windows
                    + pwm_filtered_windows + emitted_hits
               OR bytes <= 0 OR length(sha256) <> 64
               OR mtime_ns IS NULL
               OR source_commit <> '$source_commit'
               OR slurm_job_id IS NULL OR slurm_array_task_id IS NULL
               OR minimum_pwm_relative_score IS NOT NULL
               OR maximum_pwm_relative_score IS NOT NULL
               OR coordinate_mode <> 'bed' OR n_policy <> 'skip'
               OR matched_sequence_policy <> 'omitted'
        ) THEN error('file inventory lost validation or provenance fields') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM scan_run
            WHERE minimum_pwm_relative_score IS NOT NULL
               OR maximum_pwm_relative_score IS NOT NULL
               OR emitted_hit_count <= 0
               OR parquet_bytes <= 0
               OR length(scanner_sha256) <> 64
               OR scanner_build_json IS NULL
               OR duckdb_version IS NULL
               OR finalization_validation_mode <>
                    'exact_inventory_size_mtime_when_available'
        ) THEN error('scan_run lost explicit score configuration') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM genome
            WHERE genome_id <> 'synthetic_genome_v1'
               OR taxon_id <> 9606
               OR assembly_accession <> 'GCA_SYNTHETIC'
               OR ensembl_release <> 113
               OR length(fasta_sha256) <> 64
               OR length(sequence_set_sha256) <> 64
        ) THEN error('genome identity or checksums are incomplete') END;
    " >/dev/null
)

query_tool="$repository_root/scripts/query_genome_scan.py"
"$query_tool" summary --package "$run_root/package" --format json \
    --output "$temporary/summary.json"
python3 - "$temporary/summary.json" <<'PY'
import json
import sys

rows = json.load(open(sys.argv[1], encoding="utf-8"))
assert len(rows) == 1
assert rows[0]["run_id"] == "synthetic_jaspar_scan_v1"
assert rows[0]["file_count"] == 4
PY

hits_json=$(cd "$temporary" && "$query_tool" hits \
    --package "$run_root/package" --motif MA0861.2 --chrom 1 \
    --format json --limit 5)
python3 - "$hits_json" <<'PY'
import json
import sys

rows = json.loads(sys.argv[1])
assert rows
for row in rows:
    assert row["run_id"] == "synthetic_jaspar_scan_v1"
    assert row["genome_id"] == "synthetic_genome_v1"
    assert row["motif_set_id"] == "synthetic_jaspar_v1"
    assert row["motif_id"] == "MA0861.2"
    assert row["chrom"] == "1"
PY

# Catalog construction must use metadata only and therefore still work while
# the payload directory is temporarily unavailable.
mv "$run_root/package/task_data" "$temporary/task_data.hidden"
"$manager" build-catalog --run-root "$run_root" --duckdb duckdb \
    --output "$temporary/rebuilt.duckdb" >/dev/null
mv "$temporary/task_data.hidden" "$run_root/package/task_data"
[[ -s $temporary/rebuilt.duckdb ]] || {
    echo "E: Metadata-only catalog rebuild produced no database." >&2
    exit 1
}

layout_json=$(
    "$repository_root/scripts/benchmark_sparse_layout.py" \
        --package "$run_root/package" --motif MA0861.2 --chrom 1 \
        --output-dir "$temporary/layout-benchmark" --repetitions 1
)
python3 - "$layout_json" <<'PY'
import json
import sys

result = json.loads(sys.argv[1])
assert result["split_file_count"] == 2
assert result["combined_file_count"] == 1
assert result["split_bytes"] > 0
assert result["combined_bytes"] > 0
PY

# The full audit is deliberately separate and resumable. SIGUSR1 must produce
# exactly one manager progress line without terminating the audit.
"$manager" verify --run-root "$run_root" --checksums --max-files 1 \
    --max-read-mib-per-second 0.005 \
    >"$temporary/verify.stdout" 2>"$temporary/verify.stderr" &
verify_pid=$!
audit_started=0
for _ in {1..200}; do
    if grep -Fq 'I: Checksum audit started' "$temporary/verify.stderr"; then
        audit_started=1
        break
    fi
    kill -0 "$verify_pid" 2>/dev/null || break
    sleep 0.01
done
[[ $audit_started -eq 1 ]] || {
    echo "E: Checksum audit did not reach its signal-safe phase." >&2
    exit 1
}
kill -0 "$verify_pid" 2>/dev/null || {
    echo "E: Checksum audit finished before its progress signal test." >&2
    exit 1
}
kill -USR1 "$verify_pid"
kill -0 "$verify_pid" 2>/dev/null || {
    echo "E: Checksum audit did not survive SIGUSR1." >&2
    exit 1
}
wait "$verify_pid"
progress_lines=$(grep -c 'I: manager progress request' "$temporary/verify.stderr")
[[ $progress_lines -eq 1 ]] || {
    echo "E: Expected one manager progress line, got $progress_lines." >&2
    exit 1
}
grep -Eq 'phase=checksum_audit.*files_total=4.*bytes_total=[1-9][0-9]*' \
    "$temporary/verify.stderr" || {
    echo "E: Manager progress line lacks checksum totals." >&2
    exit 1
}
partial_status=$("$manager" status --run-root "$run_root")
grep -Fq $'checksum_verification_state\tpartial' <<< "$partial_status"
"$manager" verify --run-root "$run_root" --checksums >/dev/null
"$manager" verify --run-root "$run_root" --checksums >/dev/null

# A same-size payload mutation with its mtime restored must still be caught by
# the explicit SHA-256 audit, demonstrating that fast finalization and full
# integrity verification are genuinely separate checks.
corrupt_run="$temporary/corrupt-run"
mkdir -p "$corrupt_run/package"
cp -pR "$run_root/plan" "$corrupt_run/plan"
cp -pR "$run_root/package/task_data" "$corrupt_run/package/task_data"
corrupt_file=$(find "$corrupt_run/package/task_data" -name '*.parquet' -type f -print -quit)
python3 - "$corrupt_file" <<'PY'
import os
import sys

path = sys.argv[1]
metadata = os.stat(path)
with open(path, "r+b") as stream:
    original = stream.read(1)
    stream.seek(0)
    stream.write(bytes([original[0] ^ 0x01]))
os.utime(path, ns=(metadata.st_atime_ns, metadata.st_mtime_ns))
PY
if "$manager" verify --run-root "$corrupt_run" --checksums \
    >"$temporary/corrupt.stdout" 2>"$temporary/corrupt.stderr"; then
    echo "E: Checksum audit accepted a same-size payload mutation." >&2
    exit 1
fi
grep -Fq 'checksum mismatch' "$temporary/corrupt.stderr"

status=$($manager status --run-root "$run_root")
grep -Fq $'tasks_complete\t2' <<< "$status"
grep -Fq $'finalized\ttrue' <<< "$status"
grep -Fq $'checksum_verification_state\tcomplete' <<< "$status"

# A failed scanner attempt must remain in its unique staging directory and
# must never create a promoted task.
failed_run_root="$temporary/failed-run"
"$manager" prepare \
    --run-root "$failed_run_root" --run-id synthetic_failed_scan_v1 \
    --source-commit "$source_commit" \
    --genome "$temporary/genome.fna" \
    --fasta-index "$temporary/genome.fna.fai" \
    --jaspar "$temporary/motifs.jaspar" \
    --genome-id synthetic_genome_v1 --motif-set-id synthetic_jaspar_v1 \
    --taxon-id 9606 --species 'Homo sapiens synthetic' \
    --assembly-name SYN1 --assembly-accession GCA_SYNTHETIC \
    --ensembl-release 113 \
    --fasta-url https://example.invalid/genome.fna \
    --jaspar-url https://example.invalid/motifs.jaspar \
    --chrom 1 --motif-batch-size 1 --minimum-free-gib 0 >/dev/null
if "$manager" run-task --run-root "$failed_run_root" --task-index 0 \
    --scanner "$(type -P false)" --duckdb duckdb --allow-local \
    --allow-scanner-provenance-mismatch \
    >/dev/null 2>&1
then
    echo "E: Deliberately failing scanner unexpectedly succeeded." >&2
    exit 1
fi
failure_count=$(find "$failed_run_root/staging" -name failure.json -type f | wc -l)
[[ $failure_count -eq 1 ]] || {
    echo "E: Failed attempt was not retained with exactly one failure marker." >&2
    exit 1
}
if find "$failed_run_root/package/task_data" -name task_result.json \
    -type f 2>/dev/null | grep -q .
then
    echo "E: Failed task was promoted." >&2
    exit 1
fi

# A rich per-motif registry must group motifs only with identical thresholds,
# survive task execution, and become a queryable provenance table.
threshold_registry="$temporary/motif-thresholds.tsv"
cat > "$threshold_registry" <<'EOF'
motif_id	informative_threshold	informative_source	default_minimum_score	candidate_minimum_score	density_minimum_spacing_bp	density_maximum_loci	density_threshold	final_minimum_score	density_limited	density_chrom	valid_locus_starts	skipped_locus_starts	loci_at_candidate_threshold	loci_at_final_threshold	mean_spacing_bp_at_final_threshold	distribution_sha256
MA0001.1	-2	synthetic_context	-1	-2	4	2	0	0	true	1	8	0	4	1	8	1111111111111111111111111111111111111111111111111111111111111111
MA0861.2	-5	synthetic_cutrun	-1	-5	4	1	-4	-4	true	1	7	0	3	1	7	2222222222222222222222222222222222222222222222222222222222222222
EOF
python3 - "$threshold_registry" "$temporary/motif-thresholds.json" \
    "$temporary/genome.fna" "$temporary/motifs.jaspar" "$source_commit" <<'PY'
import hashlib
import json
import pathlib
import sys

thresholds = pathlib.Path(sys.argv[1])
genome = pathlib.Path(sys.argv[3])
jaspar = pathlib.Path(sys.argv[4])
metadata = {
    "schema_version": 1,
    "threshold_set_id": "synthetic_density4_v1",
    "genome_id": "synthetic_genome_v1",
    "motif_set_id": "synthetic_jaspar_v1",
    "score_mode": "log2_relative_risk",
    "pseudocount": 1,
    "genome_fasta_sha256": hashlib.sha256(genome.read_bytes()).hexdigest(),
    "density_chrom": "1",
    "density_chrom_sequence_sha256": hashlib.sha256(b"ACGTNACGT").hexdigest(),
    "jaspar_sha256": hashlib.sha256(jaspar.read_bytes()).hexdigest(),
    "source_commit": sys.argv[5],
    "orientation_aggregation": "max_score_per_alignment_span",
    "distribution_bin_width": 1,
    "candidate_formula": "min(informative_threshold, default_minimum_score)",
    "final_formula": "max(candidate_minimum_score, density_threshold)",
    "density_minimum_spacing_bp": 4,
    "default_minimum_score": -1,
    "informative_thresholds_sha256": "3" * 64,
    "negative_sensitivity_sha256": None,
    "override_thresholds_sha256": "4" * 64,
    "distribution_manifest_sha256": "5" * 64,
    "motifs": 2,
    "threshold_tsv_sha256": hashlib.sha256(thresholds.read_bytes()).hexdigest(),
}
pathlib.Path(sys.argv[2]).write_text(
    json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="ascii"
)
PY

threshold_run="$temporary/threshold-run"
threshold_task_count=$(
    "$manager" prepare \
        --run-root "$threshold_run" --run-id synthetic_threshold_scan_v1 \
        --source-commit "$source_commit" \
        --genome "$temporary/genome.fna" \
        --fasta-index "$temporary/genome.fna.fai" \
        --jaspar "$temporary/motifs.jaspar" \
        --genome-id synthetic_genome_v1 --motif-set-id synthetic_jaspar_v1 \
        --taxon-id 9606 --species 'Homo sapiens synthetic' \
        --assembly-name SYN1 --assembly-accession GCA_SYNTHETIC \
        --ensembl-release 113 \
        --fasta-url https://example.invalid/genome.fna \
        --jaspar-url https://example.invalid/motifs.jaspar \
        --chrom 1 --motif-batch-size 64 --minimum-free-gib 0 \
        --motif-thresholds "$threshold_registry" \
        --motif-threshold-metadata "$temporary/motif-thresholds.json" \
        | awk -F= '$1 == "TASK_COUNT" { print $2 }'
)
[[ $threshold_task_count == 2 ]] || {
    echo "E: Distinct thresholds were not split into two tasks." >&2
    exit 1
}
python3 - "$threshold_run/plan/scan_plan.json" <<'PY'
import json
import sys

plan = json.load(open(sys.argv[1], encoding="utf-8"))
assert plan["schema_version"] == 2
assert len(plan["threshold_policies"]) == 1
assert plan["threshold_policies"][0]["selector_type"] == "per_motif_registry"
assert {task["minimum_score"] for task in plan["tasks"]} == {-4.0, 0.0}
for task in plan["tasks"]:
    assert len(task["motif_ids"]) == 1
PY

"$manager" run-chromosome --run-root "$threshold_run" --chrom 1 \
    --scanner "$scanner" --duckdb duckdb --allow-local \
    --scratch-directory "$temporary/threshold-scratch" \
    --minimum-scratch-free-gib 0 --batch-workers 2 \
    --allow-scanner-provenance-mismatch >/dev/null
"$manager" finalize --run-root "$threshold_run" --duckdb duckdb >/dev/null
(
    cd "$threshold_run/package"
    duckdb jaspar_genome_scan.duckdb -bail -c "
        SELECT CASE WHEN (SELECT count(*) FROM scan_motif_threshold) <> 2
            THEN error('rich motif threshold table has the wrong size') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM scan_motif_threshold
            WHERE motif_id = 'MA0861.2'
              AND threshold_set_id = 'synthetic_density4_v1'
              AND candidate_minimum_score = -5
              AND density_threshold = -4
              AND final_minimum_score = -4
              AND density_limited
        ) THEN error('rich TP73 threshold provenance was lost') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM scan_motif_threshold
            WHERE motif_id = 'MA0001.1' AND final_minimum_score = 0
        ) THEN error('rich cofactor threshold was lost') END;
    " >/dev/null
)

echo "Genome-scan manager tests passed."
