#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_dilution_series.py
========================
Expand a baseline table (one row per parent library) into the per-run rows of a
dilution series (default 3 parents x 10 fractions x 3 seeds = 90 rows).

Usage:
    python make_dilution_series.py --dry-run      # preview the table
    python make_dilution_series.py --baseline baseline.csv --out manifest.csv

Key conventions:
    * `Total Reads` of a public PE run equals the SRA `spots`, i.e. read **pairs**.
      Hence target_total_pairs = round(spots x fraction); do not divide by 2 again.
    * expected_target_reads is left empty: the target read count can only be measured
      from the fraction=1.0 run (Reads_Mapped in detection_virus_summary.tsv); other points scale with it.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

COLUMNS = [
    "series_id", "parent_library_id", "library_id", "track", "gradient_index",
    "gradient_type", "downsample_fraction", "target_total_pairs", "replicate_index",
    "seed", "expected_target_reads", "actual_total_pairs", "actual_target_reads",
    "hvfinder_project_dir", "status", "notes",
]

DEFAULT_FRACTIONS = [1.0, 0.5, 0.25, 0.10, 0.05, 0.025, 0.010, 0.005, 0.0025, 0.001]


def _fmt_fraction(x: float) -> str:
    if x == 1:
        return "1.0"
    return f"{x:g}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", default="03b_baseline_three_samples.csv")
    ap.add_argument("--out", default="02_dilution_series_manifest.csv")
    ap.add_argument("--fractions", default=",".join(str(f) for f in DEFAULT_FRACTIONS))
    ap.add_argument("--seeds", default="11,22,33")
    ap.add_argument("--max-pairs", type=int, default=20_000_000,
                    help="cap on total pairs per point (inactive when the parent is smaller)")
    ap.add_argument("--project-root", default="projects")
    ap.add_argument("--min-pairs", type=int, default=100,
                    help="mark underpowered points (below this many pairs) in notes")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.baseline):
        print(f"[error] baseline not found: {args.baseline}", file=sys.stderr)
        return 2

    fractions = [float(x) for x in args.fractions.split(",") if x.strip()]
    seeds = [int(x) for x in args.seeds.split(",") if x.strip()]

    with open(args.baseline, "r", encoding="utf-8-sig") as fh:
        samples = [r for r in csv.DictReader(fh) if (r.get("sample_id") or "").strip()]

    if not samples:
        print("[error] baseline has no rows", file=sys.stderr)
        return 2

    rows = []
    for s in samples:
        sid = s["sample_id"].strip()
        srr = s["srr_id"].strip()
        spots = int(float(s.get("sra_spots_pairs") or s.get("total_reads") or 0))
        if spots <= 0:
            print(f"[skip] {sid}: no spot count")
            continue
        for gi, frac in enumerate(fractions, start=1):
            pairs = min(max(1, int(round(spots * frac))), args.max_pairs)
            for ri, seed in enumerate(seeds, start=1):
                lib = f"{sid}_G{gi:02d}_R{ri}"
                note = "seed=%d" % seed
                if frac == 1.0:
                    note += "; baseline run: fill actual_target_reads from detection_virus_summary.tsv"
                else:
                    note += "; expected_target_reads ~= baseline_target_reads x %s" % _fmt_fraction(frac)
                if pairs < args.min_pairs:
                    note += "; underpowered (<%d pairs)" % args.min_pairs
                rows.append({
                    "series_id": f"DIL_{srr}",
                    "parent_library_id": srr,
                    "library_id": lib,
                    "track": "in_silico",
                    "gradient_index": gi,
                    "gradient_type": "fraction_of_library",
                    "downsample_fraction": _fmt_fraction(frac),
                    "target_total_pairs": pairs,
                    "replicate_index": ri,
                    "seed": seed,
                    "expected_target_reads": "",
                    "actual_total_pairs": "",
                    "actual_target_reads": "",
                    "hvfinder_project_dir": f"{args.project_root.rstrip('/')}/{lib}",
                    "status": "planned",
                    "notes": note,
                })

    # print an overview
    print(f"{'sample':6s} {'fraction':>10s} {'pairs':>10s}   (x{len(seeds)} seeds)")
    for s in samples:
        sid = s["sample_id"].strip()
        spots = int(float(s.get("sra_spots_pairs") or s.get("total_reads") or 0))
        for gi, frac in enumerate(fractions, start=1):
            pairs = min(max(1, int(round(spots * frac))), args.max_pairs)
            print(f"{sid:6s} {frac:>10g} {pairs:>10d}")
        print()

    print(f"total rows: {len(rows)}  (fractions={len(fractions)}, seeds={len(seeds)}, samples={len(samples)})")

    if args.dry_run:
        print("[dry-run] nothing written")
        return 0

    with open(args.out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"[ok] {args.out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
