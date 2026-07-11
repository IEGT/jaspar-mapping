#!/usr/bin/env python3

# Fetch a human genomic region and features from Ensembl REST and write an EMBL flatfile.
# The output EMBL file can be imported into SnapGene or similar tools.
#
# Created by chatGPT-4 based on user specifications and edited by Steffen Möller.

import argparse
import datetime as _dt
import json
import re
import sys
from collections import Counter
from typing import Dict, List, Optional, Tuple



def requests_module():
    """Load the optional network dependency only after CLI parsing."""
    try:
        import requests
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "The requests package is required for Ensembl retrieval. "
            "Install the project's Python requirements first."
        ) from error
    return requests


def ensembl_host(assembly: str) -> str:
    # Ensembl provides a stable GRCh37 REST host; otherwise use current REST.
    if assembly.lower() in {"grch37", "hg19"}:
        return "https://grch37.rest.ensembl.org"
    return "https://rest.ensembl.org"


def http_get_json(url: str, headers: Dict[str, str], params: Optional[Dict[str, str]] = None):
    r = requests_module().get(url, headers=headers, params=params, timeout=60)
    if not r.ok:
        raise RuntimeError(f"GET {url} failed: HTTP {r.status_code}: {r.text[:500]}")
    return r.json()


def http_get_text(url: str, headers: Dict[str, str], params: Optional[Dict[str, str]] = None) -> str:
    r = requests_module().get(url, headers=headers, params=params, timeout=60)
    if not r.ok:
        raise RuntimeError(f"GET {url} failed: HTTP {r.status_code}: {r.text[:500]}")
    return r.text


def fetch_sequence(host: str, species: str, region: str) -> str:
    # GET /sequence/region/:species/:region returns sequence.
    url = f"{host}/sequence/region/{species}/{region}"
    headers = {"Content-Type": "text/plain", "Accept": "text/plain"}
    seq = http_get_text(url, headers=headers)
    seq = re.sub(r"\s+", "", seq).upper()
    if not re.fullmatch(r"[ACGTN]*", seq):
        raise RuntimeError("Sequence contains unexpected characters.")
    return seq


def fetch_overlap(host: str, species: str, region: str, features: List[str]) -> List[dict]:
    # GET /overlap/region/:species/:region?feature=gene&feature=exon...
    url = f"{host}/overlap/region/{species}/{region}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = []
    # requests wants dict; for repeated params use list of tuples via "params".
    # We'll build manually in URL-like style by passing list of tuples.
    params_tuples = [("feature", f) for f in features]
    r = requests_module().get(
        url, headers=headers, params=params_tuples, timeout=60
    )
    if not r.ok:
        raise RuntimeError(f"GET {url} failed: HTTP {r.status_code}: {r.text[:500]}")
    return r.json()


def revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(tr)[::-1]


def clamp_feature_to_region(f_start: int, f_end: int, region_start: int, region_end: int) -> Optional[Tuple[int, int]]:
    s = max(f_start, region_start)
    e = min(f_end, region_end)
    if s > e:
        return None
    return s, e


def embl_location(region_start: int, region_end: int, f_start: int, f_end: int, strand: int) -> str:
    # Convert genomic coordinates to 1-based coordinates relative to extracted sequence.
    # EMBL locations are 1..len, and for reverse strand use complement(a..b).
    rel_start = f_start - region_start + 1
    rel_end = f_end - region_start + 1
    loc = f"{rel_start}..{rel_end}"
    if strand == -1:
        return f"complement({loc})"
    return loc


def wrap_seq_for_embl(seq: str) -> List[str]:
    # EMBL SQ lines: 60 bases per line, grouped in blocks of 10, with a right-justified position.
    lines = []
    n = len(seq)
    for i in range(0, n, 60):
        chunk = seq[i:i+60].lower()
        blocks = [chunk[j:j+10] for j in range(0, len(chunk), 10)]
        left = " ".join(blocks)
        pos = i + len(chunk)
        # 6 blocks of 10 bases => 5 spaces between blocks => 60 + 5 = 65 chars
        lines.append(f"     {left:<65} {pos:>9d}")
    return lines


def embl_sq_header(seq: str) -> str:
    n = len(seq)
    c = Counter(seq.lower())
    # Provide counts for a,c,g,t and "other" as typical flatfiles do.
    a = c.get("a", 0)
    c_ = c.get("c", 0)
    g = c.get("g", 0)
    t = c.get("t", 0)
    other = n - (a + c_ + g + t)
    return f"SQ   Sequence {n} BP; {a} A; {c_} C; {g} G; {t} T; {other} other;"


def read_custom_tsv(path: str, separator: str,
                    from_col: str, to_col: str, name_col: Optional[str], type_col: str) -> List[dict]:
    """
    TSV columns (header required):
      type  start  end  strand  name  qualifiers_json
    - start/end are genomic coordinates (same coordinate system as your region)
    - strand is +1 or -1 (optional; default +1)
    - qualifiers_json is optional JSON object, e.g. {"note":"my feature","gene":"XYZ"}
    """
    feats = []
    with open(path, "r", encoding="utf-8") as fh:
        header_line = ""
        while True:
            header_line = fh.readline()
            if not header_line:
                raise RuntimeError("Custom TSV is empty or has no header.")
            if header_line.strip() and not header_line.lstrip().startswith("#"):
                break
        header = header_line.rstrip("\n").split(separator)
        idx = {k: i for i, k in enumerate(header)}
        required = {type_col, from_col, to_col}
        missing = required - set(idx)
        if missing:
            raise RuntimeError(f"Custom TSV missing required columns: {', '.join(sorted(missing))}")
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(separator)
            f = {
                "feature_type": parts[idx[type_col]],
                "start": int(parts[idx[from_col]]),
                "end": int(parts[idx[to_col]]),
                "strand": int(parts[idx["strand"]]) if "strand" in idx and parts[idx["strand"]] else 1,
                "name": parts[idx[name_col]] if name_col in idx and parts[idx[name_col]] else None,
                "qualifiers": {},
            }
            if "qualifiers_json" in idx and parts[idx["qualifiers_json"]]:
                f["qualifiers"] = json.loads(parts[idx["qualifiers_json"]])
            feats.append(f)
    return feats


def ft_lines(feature_key: str, location: str, qualifiers: Dict[str, str]) -> List[str]:
    lines = [f"FT   {feature_key:<15}{location}"]
    # Qualifiers: /key="value" (quote values; escape quotes)
    for k, v in qualifiers.items():
        if v is None:
            continue
        v = str(v).replace('"', '\\"')
        lines.append(f'FT                   /{k}="{v}"')
    return lines


def main():
    ap = argparse.ArgumentParser(
        description="Fetch human genomic region and features and write an EMBL flatfile (SnapGene-importable)."
    )
    ap.add_argument("--chr", required=True, help="Chromosome name, e.g. 1, 2, X, Y, MT")
    ap.add_argument("--from", dest="start", type=int, required=True, help="Start coordinate (1-based, inclusive)")
    ap.add_argument("--to", dest="end", type=int, required=True, help="End coordinate (1-based, inclusive)")
    ap.add_argument("--assembly", default="GRCh38", help="GRCh38 (default) or GRCh37")
    ap.add_argument("--species", default="human", help="Ensembl species name (default: human)")
    ap.add_argument("--out", required=True, help="Output EMBL file")
    ap.add_argument("--name", default=None, help="Record name/ID (default derived from region)")
    ap.add_argument("--accession", default="NA000000", help="Accession to place in AC line")
    ap.add_argument("--include", default="gene,exon", help="Comma-separated Ensembl overlap feature types")
    ap.add_argument("--custom-tsv", default=None, help="Optional TSV with custom features (see docstring)")
    ap.add_argument("--custom-tsv-column-separator", default="\t", help="Custom TSV column separator (default: tab)")
    ap.add_argument("--custom-tsv-column-from", default="from", help="Name of 'from' column in custom TSV (default: 'from')")
    ap.add_argument("--custom-tsv-column-to", default="to", help="Name of 'to' column in custom TSV (default: 'to')")
    ap.add_argument("--custom-tsv-column-name", default="name", help="Name of 'name' column in custom TSV (default: 'name')")
    ap.add_argument("--custom-tsv-column-type", default="type", help="Name of 'type' column in custom TSV (default: 'type')")
    ap.add_argument("--extra-basepairs", type=int, default=0, help="Extra base pairs to fetch on each side")
    ap.add_argument("--rc", action="store_true", help="Write the extracted sequence as reverse-complement")
    args = ap.parse_args()

    if args.start < 1 or args.end < 1:
        raise SystemExit("Coordinates must be positive (1-based).")
    if args.end < args.start:
        raise SystemExit("--to must be >= --from.")

    host = ensembl_host(args.assembly)
    chrom = args.chr
    region_start = args.start - args.extra_basepairs
    region_end = args.end + args.extra_basepairs

 
    custom_feats = read_custom_tsv(args.custom_tsv, separator=args.custom_tsv_column_separator,
                                   from_col=args.custom_tsv_column_from, to_col=args.custom_tsv_column_to, name_col=args.custom_tsv_column_name, type_col=args.custom_tsv_column_type) if args.custom_tsv else []
    # iterate over custom_features and extend region if needed
    for cf in custom_feats:
        fstart = int(cf["start"])
        fend = int(cf["end"])
        if fstart < region_start:
            print("I: Extending region start to include custom feature at", fstart, file=sys.stderr)
            region_start = fstart
        if fend > region_end:
            print("I: Extending region end to include custom feature at", fend, file=sys.stderr)
            region_end = fend


    if region_start < 1:
        raise SystemExit("E: After applying extra base pairs or custom features, region start is < 1.")

    if region_end - region_start > 10_000_000:
        raise SystemExit("E: Region too large (>10 Mbp) after applying extra base pairs or custom features.")

    # Ensembl region format: {chrom}:{start}..{end}:{strand}
    # We'll request strand 1, and optionally reverse-complement after retrieval.
    region = f"{chrom}:{region_start}..{region_end}:1"

    # Overlap features
    include = [x.strip() for x in args.include.split(",") if x.strip()]
    overlaps = fetch_overlap(host, args.species, region, include) if include else []

    # Fetch sequence
    seq = fetch_sequence(host, args.species, region)
    if len(seq) != (region_end - region_start + 1):
        # Ensembl usually returns exact length; if not, fail loudly.
        raise RuntimeError(f"E: Returned sequence length {len(seq)} does not match requested inclusive interval length {(region_end - region_start + 1)}.")


    # Prepare EMBL
    record_name = args.name or f"hs_chr{chrom}_{region_start}_{region_end}"
    n = len(seq)
    today = _dt.date.today().isoformat()

    # If user wants reverse complement record, flip sequence and invert feature strands.
    record_strand_flip = False
    if args.rc:
        seq = revcomp(seq)
        record_strand_flip = True

    embl = []
    embl.append(f"ID   {record_name}; SV 1; linear; genomic DNA; STD; HUM; {n} BP.")
    embl.append("XX")
    embl.append(f"AC   {args.accession};")
    embl.append("XX")
    embl.append(f"DE   Homo sapiens genomic region chr{chrom}:{region_start}-{region_end} ({args.assembly}); extracted {today}.")
    embl.append("XX")
    embl.append("FH   Key             Location/Qualifiers")
    embl.append("FH")
    embl.extend(ft_lines(
        "source",
        f"1..{n}",
        {
            "organism": "Homo sapiens",
            "mol_type": "genomic DNA",
            "db_xref": "taxon:9606",
            "note": f"Source coordinates: chr{chrom}:{region_start}-{region_end} ({args.assembly}); Ensembl REST sequence/region",
        },
    ))

    def add_overlap_feature(obj: dict):
        fstart = int(obj.get("start"))
        fend = int(obj.get("end"))
        strand = int(obj.get("strand", 1) or 1)

        clamped = clamp_feature_to_region(fstart, fend, region_start, region_end)
        if not clamped:
            return
        s, e = clamped

        # If the whole record is reverse-complemented, invert strands relative to record orientation.
        eff_strand = -strand if record_strand_flip else strand
        # If record is RC, coordinates also mirror: position p becomes (n - p + 1).
        if record_strand_flip:
            # Convert original relative positions, then mirror.
            rel_s = s - region_start + 1
            rel_e = e - region_start + 1
            # Mirror within 1..n:
            mir_s = n - rel_e + 1
            mir_e = n - rel_s + 1
            loc = f"{mir_s}..{mir_e}"
            if eff_strand == -1:
                loc = f"complement({loc})"
        else:
            loc = embl_location(region_start, region_end, s, e, eff_strand)

        ftype = obj.get("feature_type") or obj.get("type") or "misc_feature"
        # Keep gene as "gene", but map exon to misc_feature for broad compatibility.
        if ftype == "exon":
            key = "misc_feature"
        elif ftype == "gene":
            key = "gene"
        else:
            key = "misc_feature"

        qualifiers = {}
        if obj.get("external_name"):
            qualifiers["gene"] = obj["external_name"] if key == "gene" else None
            if key != "gene":
                qualifiers["note"] = f'{ftype}: {obj["external_name"]}'
        if obj.get("id"):
            qualifiers["db_xref"] = f'Ensembl:{obj["id"]}'
        if obj.get("biotype"):
            qualifiers["note"] = (qualifiers.get("note", "") + ("; " if qualifiers.get("note") else "") + f'biotype={obj["biotype"]}').strip()
        if obj.get("description"):
            qualifiers["note"] = (qualifiers.get("note", "") + ("; " if qualifiers.get("note") else "") + str(obj["description"])).strip()

        # Ensure exon is at least labeled
        if ftype == "exon" and not qualifiers.get("note"):
            qualifiers["note"] = "exon"

        embl.extend(ft_lines(key, loc, qualifiers))

    # Add overlaps (genes/exons/etc.)
    for o in overlaps:
        add_overlap_feature(o)

    # Add custom features
    for cf in custom_feats:
        fstart = int(cf["start"])
        fend = int(cf["end"])
        strand = int(cf.get("strand", 1))
        clamped = clamp_feature_to_region(fstart, fend, region_start, region_end)
        if not clamped:
            continue
        s, e = clamped

        eff_strand = -strand if record_strand_flip else strand
        if record_strand_flip:
            rel_s = s - region_start + 1
            rel_e = e - region_start + 1
            mir_s = n - rel_e + 1
            mir_e = n - rel_s + 1
            loc = f"{mir_s}..{mir_e}"
            if eff_strand == -1:
                loc = f"complement({loc})"
        else:
            loc = embl_location(region_start, region_end, s, e, eff_strand)

        key = cf.get("feature_type") or "misc_feature"
        qualifiers = dict(cf.get("qualifiers") or {})
        if cf.get("name"):
            # Put a human-visible label in /note if none provided
            qualifiers.setdefault("note", cf["name"])
        embl.extend(ft_lines(key, loc, qualifiers))

    embl.append("XX")
    embl.append(embl_sq_header(seq))
    embl.extend(wrap_seq_for_embl(seq))
    embl.append("//")

    with open(args.out, "w", encoding="utf-8") as out:
        out.write("\n".join(embl) + "\n")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
