#!/usr/bin/env bash

set -euo pipefail

if ! command -v duckdb >/dev/null 2>&1; then
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
fi

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-metadata-catalog.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

metadata="$temporary/metadata.tsv"
pssm="$temporary/motifs.jaspar"
catalog="$temporary/catalog"

cat > "$metadata" <<'EOF'
collection	tax_group	matrix_id	base_id	version	name	class	family	uniprot_ids	validation	comment	source	type	tax_id	species
CORE	vertebrates	MA0001.1	MA0001	1	ONE	Class 1	Family 1	P1	PMID1		Source 1	SELEX	9606::10090	Homo sapiens::Mus musculus
CORE	plants	MA0002.1	MA0002	1	TWO	Class 2	Family 2	P2	PMID2		Source 2	ChIP-seq	3702	Arabidopsis thaliana
CORE	vertebrates	MA0003.1	MA0003	1	THREE	Class 3	Family 3				Source 3	PBM		
EOF

cat > "$pssm" <<'EOF'
>MA0001.1	ONE
A  [ 1 0 ]
C  [ 0 1 ]
G  [ 0 0 ]
T  [ 0 0 ]
>MA0003.1	THREE
A  [ 0 1 ]
C  [ 1 0 ]
G  [ 0 0 ]
T  [ 0 0 ]
EOF

metadata_sha256=$(shasum -a 256 "$metadata" | awk '{print $1}')
"$repository_root/scripts/build_jaspar_metadata_catalog.py" \
    --metadata "$metadata" --pssm "$pssm" --output "$catalog" \
    --expected-metadata-sha256 "$metadata_sha256" \
    --expected-motif-count 2 --motif-set-id synthetic_core

database="$catalog/jaspar_metadata.duckdb"
duckdb -light-mode -readonly -batch "$database" <<'SQL'
SELECT CASE WHEN (SELECT count(*) FROM jaspar_matrix) <> 3
    THEN error('wrong matrix count') END;
SELECT CASE WHEN (SELECT count(*) FROM jaspar_matrix_species) <> 3
    THEN error('wrong exploded species count') END;
SELECT CASE WHEN (SELECT count(*) FROM jaspar_motif_set_matrix) <> 2
    THEN error('wrong motif-set membership count') END;
SELECT CASE WHEN NOT (SELECT includes_homo_sapiens FROM jaspar_matrix
                       WHERE matrix_id = 'MA0001.1')
    THEN error('human source flag was lost') END;
SELECT CASE WHEN (SELECT species_count FROM jaspar_matrix
                  WHERE matrix_id = 'MA0003.1') <> 0
    THEN error('unknown source species was not represented') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM jaspar_motif_set_matrix WHERE NOT name_matches_metadata
) THEN error('matching PSSM names were rejected') END;
SQL

python3 - "$catalog/manifest.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    manifest = json.load(handle)
assert manifest["schema_version"] == 1
assert manifest["metadata_matrix_count"] == 3
assert manifest["metadata_species_row_count"] == 3
assert manifest["motif_set_matrix_count"] == 2
assert manifest["motif_set_human_source_count"] == 1
assert manifest["motif_set_tax_group_counts"] == {"vertebrates": 2}
assert manifest["pssm_metadata_name_mismatch_count"] == 0
PY

echo "JASPAR metadata catalog tests passed."
