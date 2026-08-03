#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-density-calibration.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

scanner="$repository_root/pssm_scan"
[[ -x $scanner ]] || {
    echo "E: Build pssm_scan before this test." >&2
    exit 1
}

cat > "$temporary/genome.fna" <<'EOF'
>1 synthetic chromosome
ACGTACGTACGTACGTACGT
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

cat > "$temporary/informative.tsv" <<'EOF'
motif_id	recommended_threshold
MA0001.1	-2
EOF
cat > "$temporary/overrides.tsv" <<'EOF'
motif_id	informative_threshold	informative_source
MA0861.2	-5	synthetic_tp73_override
EOF

manager="$repository_root/scripts/manage_motif_density_calibration.py"
run_root="$temporary/run"
source_commit=0123456789abcdef0123456789abcdef01234567
batch_count=$(
    "$manager" prepare \
        --run-root "$run_root" --run-id synthetic_density_v1 \
        --source-commit "$source_commit" \
        --genome "$temporary/genome.fna" \
        --fasta-index "$temporary/genome.fna.fai" \
        --jaspar "$temporary/motifs.jaspar" \
        --informative-thresholds "$temporary/informative.tsv" \
        --override-thresholds "$temporary/overrides.tsv" \
        --threshold-set-id synthetic_density200_v1 \
        --genome-id synthetic_genome_v1 \
        --motif-set-id synthetic_jaspar_v1 \
        --chrom 1 --motif-batch-size 1 --minimum-spacing-bp 4 \
        --minimum-free-gib 0 \
        | awk -F= '$1 == "BATCH_COUNT" { print $2 }'
)
[[ $batch_count == 2 ]] || {
    echo "E: Expected two density batches, got $batch_count." >&2
    exit 1
}

"$manager" run --run-root "$run_root" --scanner "$scanner" \
    --scratch-directory "$temporary/scratch" --batch-workers 2 \
    --minimum-scratch-free-gib 0 --allow-local \
    --allow-scanner-provenance-mismatch >/dev/null
"$manager" status --run-root "$run_root" > "$temporary/status.tsv"
grep -Fq $'batches_complete\t2' "$temporary/status.tsv"
grep -Fq $'batches_pending\t0' "$temporary/status.tsv"

"$manager" finalize --run-root "$run_root" >/dev/null
"$manager" finalize --run-root "$run_root" >/dev/null

python3 - "$run_root" <<'PY'
import csv
import hashlib
import json
import pathlib
import sys

run = pathlib.Path(sys.argv[1])
final = run / "final"
manifest = json.loads((final / "manifest.json").read_text(encoding="utf-8"))
threshold_path = final / manifest["threshold_tsv"]
metadata_path = final / manifest["threshold_json"]
assert hashlib.sha256(threshold_path.read_bytes()).hexdigest() == manifest["threshold_tsv_sha256"]
assert hashlib.sha256(metadata_path.read_bytes()).hexdigest() == manifest["threshold_json_sha256"]

with threshold_path.open(newline="") as stream:
    rows = list(csv.DictReader(stream, delimiter="\t"))
assert {row["motif_id"] for row in rows} == {"MA0861.2", "MA0001.1"}
for row in rows:
    assert float(row["final_minimum_score"]) >= float(row["candidate_minimum_score"])
    assert int(row["loci_at_final_threshold"]) <= int(row["density_maximum_loci"])
    if int(row["loci_at_final_threshold"]):
        assert float(row["mean_spacing_bp_at_final_threshold"]) >= 4
assert next(row for row in rows if row["motif_id"] == "MA0861.2")["informative_source"] == "synthetic_tp73_override"

metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
assert metadata["threshold_set_id"] == "synthetic_density200_v1"
assert metadata["orientation_aggregation"] == "max_score_per_alignment_span"
assert metadata["density_minimum_spacing_bp"] == 4
assert len(metadata["genome_fasta_sha256"]) == 64
assert len(metadata["density_chrom_sequence_sha256"]) == 64
assert len(metadata["jaspar_sha256"]) == 64
assert metadata["source_commit"] == "0123456789abcdef0123456789abcdef01234567"
assert metadata["total_loci_at_final_threshold_density_chrom"] == sum(
    int(row["loci_at_final_threshold"]) for row in rows
)
assert metadata["maximum_orientation_rows_at_final_threshold_density_chrom"] == 2 * metadata["total_loci_at_final_threshold_density_chrom"]

with (final / "distribution_manifest.tsv").open(newline="") as stream:
    distributions = list(csv.DictReader(stream, delimiter="\t"))
assert len(distributions) == 2
for row in distributions:
    path = (final / row["path"]).resolve()
    assert path.is_file()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
PY

echo "Motif-density calibration tests passed."
