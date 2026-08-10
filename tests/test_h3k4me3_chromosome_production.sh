#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-h3-chromosome.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping H3 chromosome-production test." >&2
    exit 0
}

annotation="$temporary/annotation"
scan="$temporary/scan"
stage="$temporary/stage"
mkdir -p "$annotation/final" "$scan" "$stage"
touch "$temporary/anchor-chr1.parquet" "$temporary/anchor-chr2.parquet"
printf '{}\n' > "$annotation/final/manifest.json"
printf '{}\n' > "$scan/manifest.json"
printf '{}\n' > "$stage/input_manifest.json"
printf 'tracks\n' > "$temporary/tracks.tsv"
printf 'thresholds\n' > "$temporary/thresholds.tsv"

"$duckdb" -light-mode -batch "$annotation/final/context.duckdb" >/dev/null <<SQL
CREATE TABLE context_file_inventory AS
SELECT * FROM (VALUES
  ('$temporary/anchor-chr1.parquet', 'tp73_context_anchor', '1', true),
  ('$temporary/anchor-chr2.parquet', 'tp73_context_anchor', '2', true)
) AS t(absolute_path, dataset, chrom, is_parquet);
SQL

python3 - "$repository_root" "$duckdb" "$annotation" "$scan" "$stage" \
    "$temporary" <<'PY'
import argparse
import importlib.util
import pathlib
import sys

root, duckdb, annotation, scan, stage, temporary = map(pathlib.Path, sys.argv[1:])
spec = importlib.util.spec_from_file_location(
    "h3_production", root / "scripts/run_h3k4me3_chr1_production.py"
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

anchor = module.resolve_anchor(duckdb, annotation, "2")
assert anchor == (temporary / "anchor-chr2.parquet").resolve()

captured = []
module.run = lambda command: captured.append(command)
arguments = argparse.Namespace(
    run_root=temporary / "run", source=root, scan_package=scan,
    chrom="2", duckdb=duckdb,
)
output = module.stage_scan_inputs(arguments, ["MA0024.3"])
assert output == temporary / "run/input/chrom-2_cofactor_panel"
assert captured and captured[0][captured[0].index("--chrom") + 1] == "2"

common = dict(
    source=root, annotation_run=annotation, scan_package=scan,
    track_manifest=temporary / "tracks.tsv", track_root=temporary,
    thresholds=temporary / "thresholds.tsv", minimum_anchor_score=-1.0,
    source_score_floor=-1.0, context_flank=150,
    primary_window="flank_150_1000", pseudocount=1.0, block_size=5_000_000,
    include_occupancy_analysis=True,
)
chr1 = argparse.Namespace(**common, chrom="1", chrom_length=248_956_422)
chr2 = argparse.Namespace(**common, chrom="2", chrom_length=242_193_529)
contract1 = module.build_contract(
    chr1, temporary / "anchor-chr1.parquet", stage, ["MA0024.3"], "0" * 40
)
contract2 = module.build_contract(
    chr2, temporary / "anchor-chr2.parquet", stage, ["MA0024.3"], "0" * 40
)
assert contract1["analysis"] == "tp73_h3k4me3_change_chr1"
assert contract1["chrom"] == "1" and contract1["chrom_length"] == 248_956_422
assert contract1["analysis_role"] == "chromosome_1_development_analysis"
assert contract2["analysis"] == "tp73_h3k4me3_change_chr2"
assert contract2["chrom"] == "2" and contract2["chrom_length"] == 242_193_529
assert contract2["analysis_role"] == (
    "held_out_chromosome_validation_thresholds_fixed_on_chr1"
)
assert contract2["include_occupancy_analysis"] is True
assert "held-out chromosome 2" in module.occupancy_inference_status("2")
PY

if bash "$repository_root/scripts/submit_h3k4me3_chr1_production_slurm.sh" \
    --run-root /data/sm718/run --annotation-run /data/sm718/annotation \
    --scan-package /data/sm718/scan --track-root /data/sm718/tracks \
    --runtime-prefix /data/sm718/runtime --chrom 2 --dry-run \
    > "$temporary/missing-length.out" 2> "$temporary/missing-length.err"; then
    echo "E: chromosome-2 submission accepted no chromosome length" >&2
    exit 1
fi
grep -Fq -- '--chrom-length is required' "$temporary/missing-length.err"

echo "I: H3 chromosome-production tests passed." >&2
