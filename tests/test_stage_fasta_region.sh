#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-stage-fasta.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

cat > "$temporary/genome.fa" <<'EOF'
>chr1 first
AAAA
CCCC
>chr2 second
AACCGG
TT
EOF
python3 "$repository_root/scripts/build_fasta_index.py" \
    "$temporary/genome.fa" --output "$temporary/genome.fa.fai" >/dev/null
expected_sha=$(python3 -c 'import hashlib; print(hashlib.sha256(b"AACCGGTT").hexdigest())')

metadata=$(
    "$repository_root/scripts/stage_fasta_region.py" \
        --genome "$temporary/genome.fa" \
        --fasta-index "$temporary/genome.fa.fai" \
        --sequence chr2 --output "$temporary/chr2.fa" \
        --expected-length 8 --expected-sha256 "$expected_sha"
)
[[ $(cat "$temporary/chr2.fa") == $'>chr2\nAACCGGTT' ]] || {
    echo "E: Staged FASTA sequence is incorrect." >&2
    exit 1
}
[[ $(cat "$temporary/chr2.fa.fai") == $'chr2\t8\t6\t8\t9' ]] || {
    echo "E: Staged FASTA index is incorrect." >&2
    exit 1
}
grep -Fq '"sequence": "chr2"' <<< "$metadata"
grep -Fq '"length": 8' <<< "$metadata"

python3 "$repository_root/scripts/build_fasta_index.py" \
    "$temporary/chr2.fa" --output "$temporary/rebuilt.fai" >/dev/null
cmp "$temporary/chr2.fa.fai" "$temporary/rebuilt.fai"

if "$repository_root/scripts/stage_fasta_region.py" \
    --genome "$temporary/genome.fa" --fasta-index "$temporary/genome.fa.fai" \
    --sequence chr1 --output "$temporary/bad.fa" \
    --expected-sha256 "$expected_sha" >/dev/null 2>&1; then
    echo "E: Scratch staging accepted the wrong planned sequence checksum." >&2
    exit 1
fi
[[ ! -e $temporary/bad.fa ]] || {
    echo "E: Failed scratch staging published its final FASTA path." >&2
    exit 1
}

echo "FASTA scratch-staging tests passed."
