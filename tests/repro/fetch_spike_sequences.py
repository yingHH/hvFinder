#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_spike_sequences.py
========================
Download the reference sequences required for C3 (in-silico near-neighbour spike-in).

The C3 rows of `09_public_negative_controls.csv` have `source_type=in_silico_mix`:
there is no wet-lab sample and no SRA run, so what is needed are **viral references**
(GenBank/RefSeq), spiked into a negative background as designed in `02b_insilico_mix_manifest.csv`.

Usage:
    python fetch_spike_sequences.py                 # five near-neighbour viruses, all RefSeq segments
    python fetch_spike_sequences.py --dry-run       # list what would be downloaded
    python fetch_spike_sequences.py --max-per-virus 1
    python fetch_spike_sequences.py --outdir derived/spike_sequences

Output:
    derived/spike_sequences/<slug>_<accession>.fasta
    derived/spike_sequences_manifest.csv    virus / taxid / accession (version kept) / length / definition / download time

Conventions:
    * read-only NCBI E-utilities only; without an API key the rate is about 3 req/s;
    * RefSeq records only; the accession **keeps its version** in the manifest (traceability);
    * these sequences are only the parents for the in-silico spike-in and do not enter the FPR denominator;
    * like `resolve_public_controls.py`, this script never modifies existing data and only writes to the output directory.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#: default near-neighbour viruses (taxids checked against NCBI Taxonomy)
DEFAULT_SOURCES = [
    {"virus": "Dugbe virus", "acronym": "DUGV", "taxid": "3052514"},
    {"virus": "Nairobi sheep disease virus", "acronym": "NSDV", "taxid": "194540"},
    {"virus": "Hazara virus", "acronym": "HAZV", "taxid": "11596"},
    {"virus": "Bhanja virus", "acronym": "BHAV", "taxid": "1213620"},
    {"virus": "Dabie bandavirus", "acronym": "SFTSV", "taxid": "2748958"},
]

MANIFEST_COLUMNS = ["virus", "acronym", "taxid", "accession", "length",
                    "definition", "fetched_at", "fasta_path"]


def _get(url: str, retries: int = 3, timeout: int = 90) -> str:
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"request failed: {url} ({last})")


def esearch(db: str, term: str, retmax: int = 20) -> list[str]:
    q = urllib.parse.urlencode({"db": db, "term": term, "retmax": retmax,
                                "retmode": "json"})
    import json
    data = json.loads(_get(f"{EUTILS}/esearch.fcgi?{q}"))
    return data["esearchresult"].get("idlist", [])


def esummary(db: str, ids: list[str]) -> dict:
    import json
    q = urllib.parse.urlencode({"db": db, "id": ",".join(ids), "retmode": "json"})
    return json.loads(_get(f"{EUTILS}/esummary.fcgi?{q}"))["result"]


def efetch_fasta(db: str, ids: list[str]) -> str:
    q = urllib.parse.urlencode({"db": db, "id": ",".join(ids),
                                "rettype": "fasta", "retmode": "text"})
    return _get(f"{EUTILS}/efetch.fcgi?{q}")


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def load_sources(path: str | None) -> list[dict]:
    if not path:
        return [dict(s) for s in DEFAULT_SOURCES]
    with open(path, "r", encoding="utf-8-sig") as fh:
        rows = [r for r in csv.DictReader(fh) if (r.get("taxid") or "").strip()]
    if not rows:
        raise SystemExit(f"[error] {path}: no line with a taxid")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fetch near-neighbour virus reference sequences (C3)")
    ap.add_argument("--outdir", default="derived/spike_sequences",
                    help="output directory (default derived/spike_sequences)")
    ap.add_argument("--manifest-out", default="derived/spike_sequences_manifest.csv")
    ap.add_argument("--sources", default=None,
                    help="custom manifest CSV (columns: virus,acronym,taxid)")
    ap.add_argument("--max-per-virus", type=int, default=0,
                    help="maximum records per virus (0 = all RefSeq segments)")
    ap.add_argument("--include-non-refseq", action="store_true",
                    help="do not add refseq[filter] (RefSeq only by default)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    sources = load_sources(args.sources)
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.dirname(args.manifest_out) or ".", exist_ok=True)

    rows: list[dict] = []
    for src in sources:
        virus = (src.get("virus") or src.get("name") or "").strip()
        acronym = (src.get("acronym") or "").strip()
        taxid = str(src.get("taxid") or "").strip()
        if not (virus and taxid.isdigit()):
            print(f"[skip] {virus or '<unnamed>'}: missing virus/taxid")
            continue
        term = f"txid{taxid}[Organism]"
        if not args.include_non_refseq:
            term += " AND refseq[filter]"
        try:
            ids = esearch("nuccore", term, retmax=50)
        except Exception as exc:  # noqa: BLE001
            print(f"[warn] {virus}: esearch failed ({exc})")
            continue
        if not ids:
            print(f"[warn] {virus}: no nucleotide record for taxid {taxid}")
            continue
        if args.max_per_virus > 0:
            ids = ids[:args.max_per_virus]
        try:
            summary = esummary("nuccore", ids)
        except Exception as exc:  # noqa: BLE001
            print(f"[warn] {virus}: esummary failed ({exc})")
            continue

        records = []
        for uid in ids:
            rec = summary.get(uid) or {}
            acc = rec.get("accessionversion") or ""
            if not acc:
                continue
            records.append((acc, rec))
        if not records:
            print(f"[warn] {virus}: no usable accession")
            continue

        if args.dry_run:
            for acc, rec in records:
                print(f"[dry] {virus:32s} {acc:12s} len={rec.get('slen')} "
                      f"{(rec.get('title') or '')[:60]}")
            time.sleep(0.4)
            continue

        try:
            fasta = efetch_fasta("nuccore", [u for u in ids])
        except Exception as exc:  # noqa: BLE001
            print(f"[warn] {virus}: efetch failed ({exc})")
            continue

        # efetch returns several FASTA records; split by accession
        blocks: dict[str, list[str]] = {}
        current = None
        for line in fasta.splitlines():
            if line.startswith(">"):
                current = line[1:].split()[0]
                blocks[current] = [line]
            elif current:
                blocks[current].append(line)

        for acc, rec in records:
            body = None
            for key, lines in blocks.items():
                if key.split(".")[0] == acc.split(".")[0]:
                    body = lines
                    break
            if not body:
                print(f"[warn] {virus}: FASTA for {acc} not found in efetch output")
                continue
            name = f"{slugify(acronym or virus)}_{acc}.fasta"
            path = os.path.join(args.outdir, name)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("\n".join(body) + "\n")
            rows.append({
                "virus": virus,
                "acronym": acronym,
                "taxid": taxid,
                "accession": acc,
                "length": rec.get("slen", ""),
                "definition": (rec.get("title") or "")[:120],
                "fetched_at": datetime.now().strftime("%Y-%m-%d"),
                "fasta_path": path,
            })
        print(f"[ok] {virus}: {len(records)} record(s)")
        time.sleep(0.4)

    if args.dry_run:
        return 0

    if not rows:
        print("\nNo sequence was downloaded.")
        return 1

    with open(args.manifest_out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"\n[ok] {args.manifest_out} ({len(rows)} rows)")
    print(f"[ok] sequences in {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
