#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_dilution_fastq.py
======================
Down-sample the 90 dilution points of `02_dilution_series_manifest.csv`
(3 series x 10 levels x 3 replicates) from the parent libraries with `seqtk sample` into paired FASTQ.

Usage:
    # preview the commands (seqtk not needed)
    python3 make_dilution_fastq.py --dry-run \
        --parent SRR29741326=/data/CCHFV/SRR29741326_1.fastq.gz,/data/CCHFV/SRR29741326_2.fastq.gz \
        ...

    # only the first three points of one series
    python3 make_dilution_fastq.py --series DIL_SRR25133377 --limit 3 --dry-run

    # actually generate (4 parallel jobs)
    python3 make_dilution_fastq.py --outdir derived/dilution_fastq --jobs 4 --parents-csv parents.csv

Key conventions:
--------
* **pairs stay paired**: R1 and R2 of the same library must use the **same seed**.
  `seqtk sample` draws one random number per record; equal record counts plus the same seed pick the same records.
* `S1_G01_R1/R2/R3` in the manifest are **replicate numbers**, not read1/read2;
  read1/read2 is given by the output file name `_1.fastq.gz` / `_2.fastq.gz`.
* `downsample_fraction` is passed to seqtk directly; `target_total_pairs = round(spots x fraction)`
  is precomputed by `expand_dilution_series.py`; this script only samples.
* output is written to `.tmp` first and then renamed, so an interrupted run leaves no half file.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MANIFEST = os.path.join(HERE, "02_dilution_series_manifest.csv")
MATES = (1, 2)


def _series_srr(series_id: str) -> str:
    """`DIL_SRR29741326` → `SRR29741326`。"""
    return series_id.split("_", 1)[1] if "_" in series_id else series_id


def read_manifest(path: str) -> list[dict]:
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        return [r for r in csv.DictReader(fh) if (r.get("library_id") or "").strip()]


def read_parents_csv(path: str) -> dict:
    """`srr_id,r1,r2` → {srr: (r1, r2)}。"""
    out = {}
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            srr = (row.get("srr_id") or row.get("parent_library_id") or "").strip()
            r1 = (row.get("r1") or "").strip()
            r2 = (row.get("r2") or "").strip()
            if srr and r1 and r2:
                out[srr] = (r1, r2)
    return out


def parse_parents(items: list[str]) -> dict:
    """`--parent SRR1=/a_1.fq.gz,/a_2.fq.gz` → {SRR1: (r1, r2)}。"""
    out = {}
    for item in items or []:
        key, _, value = item.partition("=")
        parts = [p.strip() for p in value.split(",") if p.strip()]
        if not key.strip() or len(parts) != 2:
            raise SystemExit(f"ERROR: bad --parent {item!r} (want SRR=/r1.fastq.gz,/r2.fastq.gz)")
        out[key.strip()] = (parts[0], parts[1])
    return out


def out_paths(outdir: str, library_id: str) -> tuple[str, str]:
    return (os.path.join(outdir, f"{library_id}_1.fastq.gz"),
            os.path.join(outdir, f"{library_id}_2.fastq.gz"))


def build_command(seqtk: str, seed: str, src: str, fraction: str,
                  dst: str) -> str:
    """Build one down-sampling command.

    Pitfall: `seqtk sample` guesses whether its third argument is a **fraction or a count** --
    `1.0` is taken as **1 read** (one record in practice), not 100%.
    Therefore `fraction >= 1.0` always copies the whole file and never calls seqtk.
    """
    tmp = dst + ".tmp"
    if float(fraction) >= 1.0:
        return (f"cp -f {shlex.quote(src)} {shlex.quote(tmp)} "
                f"&& mv {shlex.quote(tmp)} {shlex.quote(dst)}")
    return (f"set -o pipefail; {shlex.quote(seqtk)} sample -s {int(seed)} "
            f"{shlex.quote(src)} {fraction} | gzip -c > {shlex.quote(tmp)} "
            f"&& mv {shlex.quote(tmp)} {shlex.quote(dst)}")


def count_records(path: str) -> int:
    """Count FASTQ records (.gz and plain text)."""
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for _ in fh:
            n += 1
    return n // 4


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=DEFAULT_MANIFEST)
    ap.add_argument("--outdir", default=os.path.join(HERE, "derived", "dilution_fastq"))
    ap.add_argument("--parent", action="append", default=[],
                    help="SRR=/r1.fastq.gz,/r2.fastq.gz (repeatable)")
    ap.add_argument("--parents-csv", default=None,
                    help="CSV with columns srr_id,r1,r2")
    ap.add_argument("--seqtk", default="seqtk")
    ap.add_argument("--series", default="",
                    help="comma-separated series_id or SRR id to keep")
    ap.add_argument("--limit", type=int, default=0, help="only the first N points")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true", help="overwrite existing outputs")
    ap.add_argument("--skip-parent-check", action="store_true",
                    help="do not compare R1/R2 record counts of the parents")
    ap.add_argument("--allow-unequal-parents", action="store_true",
                    help="proceed even if a parent pair has unequal record counts")
    ap.add_argument("--skip-pair-check", action="store_true",
                    help="do not verify that generated R1/R2 have equal counts")
    ap.add_argument("--manifest-out", default=None,
                    help="record generated files (default: <outdir>/dilution_fastq_manifest.tsv)")
    args = ap.parse_args(argv)

    rows = read_manifest(args.manifest)
    if args.series:
        wanted = {s.strip() for s in args.series.split(",") if s.strip()}
        rows = [r for r in rows
                if r["series_id"] in wanted or _series_srr(r["series_id"]) in wanted]
    if args.limit > 0:
        rows = rows[:args.limit]
    if not rows:
        print("ERROR: no manifest row selected", file=sys.stderr)
        return 2

    parents = dict(read_parents_csv(args.parents_csv)) if args.parents_csv else {}
    parents.update(parse_parents(args.parent))

    todo, missing_parent = [], []
    for r in rows:
        srr = _series_srr(r["series_id"])
        if srr not in parents:
            missing_parent.append(srr)
            continue
        r1, r2 = parents[srr]
        dst1, dst2 = out_paths(args.outdir, r["library_id"])
        if not args.force and os.path.exists(dst1) and os.path.exists(dst2):
            continue
        todo.append((r, r1, r2, dst1, dst2))

    if missing_parent:
        uniq = sorted(set(missing_parent))
        print(f"ERROR: no --parent given for {', '.join(uniq)}", file=sys.stderr)
        print("       pass --parent SRR=/r1.fastq.gz,/r2.fastq.gz or --parents-csv",
              file=sys.stderr)
        return 2

    # paired length: take the shorter side. Extra trailing R1 records have no mate and are
    # dropped by the first fastp step (warning only); R2 longer than R1 means truncation -> error.
    if not args.skip_parent_check and not args.dry_run:
        bad = []
        for srr in sorted({_series_srr(r["series_id"]) for r in rows}):
            if srr not in parents:
                continue
            r1, r2 = parents[srr]
            if not (os.path.exists(r1) and os.path.exists(r2)):
                continue
            n1, n2 = count_records(r1), count_records(r2)
            print(f"parent {srr}: R1={n1:,} pairs  R2={n2:,} pairs")
            if n1 > n2:
                # trailing R1 records have no mate; pairing follows R2, the first fastp step
                # drops them -- not data loss, so warn instead of failing.
                print(f"WARN: {srr}: R1 has {n1 - n2:,} extra record(s) after R2 "
                      f"ends; pairs are limited "
                      f"to {n2:,}", file=sys.stderr)
            elif n2 > n1:
                bad.append((srr, n1, n2))
        if bad and not args.allow_unequal_parents:
            print("ERROR: parent R2 has more records than R1 (R1 looks "
                  "truncated):", file=sys.stderr)
            for srr, n1, n2 in bad:
                print(f"  {srr}: R1={n1} R2={n2} (R2-R1={n2 - n1})",
                      file=sys.stderr)
            print("       re-download the mate first, or pass "
                  "--allow-unequal-parents (pairs will be truncated to the "
                  "shorter mate downstream)", file=sys.stderr)
            return 2

    total_pairs = sum(int(r["target_total_pairs"]) for r, *_ in todo)
    print(f"manifest : {args.manifest}")
    print(f"selected : {len(rows)} point(s) | to generate: {len(todo)} "
          f"| skip existing: {len(rows) - len(todo)}")
    print(f"outdir   : {args.outdir}")
    print(f"pairs    : {total_pairs:,} (~{total_pairs / 1e6:.2f} M pairs, PE)")

    cmds = []
    for r, r1, r2, dst1, dst2 in todo:
        for src, dst in ((r1, dst1), (r2, dst2)):
            cmds.append((r, build_command(args.seqtk, r["seed"], src,
                                          r["downsample_fraction"], dst)))

    if args.dry_run:
        for i, (r, cmd) in enumerate(cmds):
            if i >= 6:
                print(f"... and {len(cmds) - 6} more command(s)")
                break
            print(f"DRY  {r['library_id']}: {cmd}")
        return 0

    os.makedirs(args.outdir, exist_ok=True)
    failures = []

    def run_one(item):
        r, cmd = item
        res = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        if res.returncode != 0:
            failures.append((r["library_id"], (res.stderr or "")[-400:]))
        return res.returncode == 0

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        list(pool.map(run_one, cmds))

    if not args.skip_pair_check:
        mism = []
        for r, _r1, _r2, dst1, dst2 in todo:
            if not (os.path.exists(dst1) and os.path.exists(dst2)):
                continue
            n1, n2 = count_records(dst1), count_records(dst2)
            if n1 != n2:
                mism.append((r["library_id"], n1, n2))
        if mism:
            print(f"WARN: {len(mism)} point(s) have unequal R1/R2 record counts:",
                  file=sys.stderr)
            for lib, n1, n2 in mism[:8]:
                print(f"  {lib}: R1={n1} R2={n2} (diff={n1 - n2})", file=sys.stderr)

    man_out = args.manifest_out or os.path.join(args.outdir, "dilution_fastq_manifest.tsv")
    with open(man_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["series_id", "library_id", "gradient_index", "fraction",
                    "seed", "target_pairs", "r1", "r2"])
        for r, _r1, _r2, dst1, dst2 in todo:
            w.writerow([r["series_id"], r["library_id"], r["gradient_index"],
                        r["downsample_fraction"], r["seed"],
                        r["target_total_pairs"], dst1, dst2])
    print(f"manifest written: {man_out}")

    if failures:
        print(f"FAILED {len(failures)} file(s):", file=sys.stderr)
        for lib, err in failures[:5]:
            print(f"  {lib}: {err}", file=sys.stderr)
        return 1
    print(f"OK: generated {len(cmds)} file(s) for {len(todo)} point(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
