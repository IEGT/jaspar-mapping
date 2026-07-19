#!/usr/bin/env bash

set -euo pipefail

case $(uname -s) in
    Linux|Darwin|FreeBSD) ;;
    *)
        echo "Progress-signal integration test skipped on non-POSIX platform."
        exit 0
        ;;
esac

repository_root=$(cd "$(dirname "$0")/.." && pwd)
scanner=${PSSM_SCAN_BIN:-"$repository_root/pssm_scan"}
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-progress-signal.XXXXXX")
scanner_pid=""

cleanup() {
    if [[ -n $scanner_pid ]] && kill -0 "$scanner_pid" 2>/dev/null; then
        kill -TERM "$scanner_pid" 2>/dev/null || true
        wait "$scanner_pid" 2>/dev/null || true
    fi
    rm -rf "$temporary"
}
trap cleanup EXIT HUP INT TERM

fail() {
    echo "E: $*" >&2
    if [[ -s $temporary/stderr.log ]]; then
        tail -n 20 "$temporary/stderr.log" >&2
    fi
    exit 1
}

[[ -x $scanner ]] || fail "Scanner is not executable: $scanner"

genome="$temporary/genome.fa"
fasta_index="$genome.fai"
pssm="$temporary/motif.pfm"
output_dir="$temporary/output"

sequence_line=ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
awk -v line="$sequence_line" 'BEGIN {
    for (chromosome = 1; chromosome <= 3; ++chromosome) {
        print ">chr" chromosome
        for (line_number = 0; line_number < 12000; ++line_number) {
            print line
        }
    }
}' > "$genome"

python3 "$repository_root/scripts/build_fasta_index.py" \
    "$genome" --output "$fasta_index" >/dev/null

cat > "$pssm" <<'EOF'
>SIGNAL.1 SIGNAL_TEST
A [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]
C [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]
G [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]
T [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]
EOF

"$scanner" -v \
    --genome "$genome" \
    --fasta-index "$fasta_index" \
    --pssm "$pssm" \
    --motif SIGNAL.1 \
    --outdir "$output_dir" \
    --threshold 100 \
    --strand both \
    --coordinate-mode bed \
    > "$temporary/stdout.log" 2> "$temporary/stderr.log" &
scanner_pid=$!

scan_started=0
for ((attempt = 0; attempt < 1000; ++attempt)); do
    if grep -Fq 'I: Scanning chromosome chr1 for motif SIGNAL_TEST' \
        "$temporary/stderr.log"; then
        scan_started=1
        break
    fi
    kill -0 "$scanner_pid" 2>/dev/null || fail "Scanner exited before its scan loop started."
    sleep 0.01
done
[[ $scan_started -eq 1 ]] || fail "Timed out waiting for the scan loop."

kill -USR1 "$scanner_pid" || fail "Could not send SIGUSR1 to scanner PID $scanner_pid."
sleep 0.02
kill -0 "$scanner_pid" 2>/dev/null || fail "Scanner did not remain alive after SIGUSR1."

progress_seen=0
for ((attempt = 0; attempt < 1000; ++attempt)); do
    progress_count=$(awk 'index($0, "I: progress request") { count++ }
        END { print count + 0 }' "$temporary/stderr.log")
    if [[ $progress_count -ge 1 ]]; then
        progress_seen=1
        break
    fi
    kill -0 "$scanner_pid" 2>/dev/null || fail "Scanner exited without reporting requested progress."
    sleep 0.01
done
[[ $progress_seen -eq 1 ]] || fail "Timed out waiting for requested progress."

if ! wait "$scanner_pid"; then
    scanner_pid=""
    fail "Scanner failed after handling SIGUSR1."
fi
scanner_pid=""

progress_count=$(awk 'index($0, "I: progress request") { count++ }
    END { print count + 0 }' "$temporary/stderr.log")
[[ $progress_count -eq 1 ]] || fail "Expected one progress line, found $progress_count."

progress_line=$(grep -o 'I: progress request.*' "$temporary/stderr.log")
grep -Eq 'phase=scan operation=scan' <<< "$progress_line" ||
    fail "Progress line lacks its scan phase."
grep -Eq 'motif_index=1 motif_total=1' <<< "$progress_line" ||
    fail "Progress line lacks global motif indices."
grep -Eq 'chromosome_index=[123] chromosome_total=3' <<< "$progress_line" ||
    fail "Progress line lacks global chromosome indices."
grep -Eq 'emitted_hits=0' <<< "$progress_line" ||
    fail "Progress line lacks the cumulative emitted-hit count."
grep -Eq 'elapsed_seconds=[0-9]+\.[0-9]{3}' <<< "$progress_line" ||
    fail "Progress line lacks elapsed time."

echo "Progress-signal integration test passed."
