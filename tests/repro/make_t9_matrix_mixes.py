#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Build the T9 environmental-matrix in-silico mixes.

T9 asks whether a complex background (tick virome, soil metagenome) hides a
true signal or produces false calls.  Ten mixes are declared here: in-panel
viruses (CCHFV, HTNV, SEOV) spiked into the RNA/DNA virome of *Haemaphysalis
longicornis* or into a soil metagenome, at 10^2 / 10^3 / 10^4 read pairs, plus
two unspiked background controls.

The heavy lifting is done by ``make_in_silico_mixes.py`` (same directory); this
script only enumerates the matrix and calls it once per mix, so a single
command reproduces every sample of the experiment.

Inputs
------
--bg-dir      directory holding the background FASTQ files, named
              ``<accession>_1.fastq.gz`` / ``<accession>_2.fastq.gz``.
              If they are missing and ``<accession>.sra`` is available in
              ``--sra-dir``, the first ``--spots`` spots are extracted with
              ``fastq-dump`` (PE150 libraries, no full conversion).
--spike-root  directory with one sub-directory per virus, named after the
              acronym in lower case (``cchfv/``, ``htnv/``, ``seov/``), each
              holding the RefSeq segments as ``*.fasta``.

Outputs (``<outdir>/``)
-----------------------
``<mix_id>_R1.fastq.gz``, ``<mix_id>_R2.fastq.gz`` and ``<mix_id>.report.json``
per mix, plus ``t9_matrix_mixes.csv`` summarising what was produced.

Example
-------
    python3 make_t9_matrix_mixes.py \
        --bg-dir  /data/t9_bg \
        --sra-dir /data/sra \
        --spike-root /data/spike_sequences \
        --outdir  derived/t9_mixes
"""
import argparse
import csv
import os
import pathlib
import shutil
import subprocess
import sys

THIS_DIR = pathlib.Path(__file__).resolve().parent
GENERATOR = THIS_DIR / "make_in_silico_mixes.py"

#: (mix_id, background accession, spike acronym (lower case) or "", pairs, expected call)
MIXES = [
    ("T9-01a", "SRR40755179", "cchfv", 100, "POSITIVE"),
    ("T9-01b", "SRR40755179", "cchfv", 1000, "POSITIVE"),
    ("T9-01c", "SRR40755179", "cchfv", 10000, "POSITIVE"),
    ("T9-02", "SRR40754971", "cchfv", 1000, "POSITIVE"),
    ("T9-03a", "SRR41515514", "htnv", 100, "POSITIVE"),
    ("T9-03b", "SRR41515514", "htnv", 1000, "POSITIVE"),
    ("T9-03c", "SRR41515514", "htnv", 10000, "POSITIVE"),
    ("T9-04", "SRR41515514", "seov", 1000, "POSITIVE"),
    ("T9-05", "SRR40755179", "", 0, "NEGATIVE"),
    ("T9-06", "SRR41515514", "", 0, "NEGATIVE"),
]


def extract_spots(accession, sra_dir, bg_dir, spots, threads=8):
    """Extract ``spots`` spots of one SRA run as gzipped paired FASTQ."""
    r1 = pathlib.Path(bg_dir) / f"{accession}_1.fastq.gz"
    r2 = pathlib.Path(bg_dir) / f"{accession}_2.fastq.gz"
    if r1.exists() and r2.exists():
        return r1, r2
    if not sra_dir:
        raise SystemExit(f"missing background FASTQ {r1} (no --sra-dir given)")
    sra = pathlib.Path(sra_dir) / accession / f"{accession}.sra"
    if not sra.exists():
        sra = pathlib.Path(sra_dir) / f"{accession}.sra"
    if not sra.exists():
        raise SystemExit(f"missing SRA file for {accession} under {sra_dir}")
    if not shutil.which("fastq-dump"):
        raise SystemExit("fastq-dump not found; install sra-tools or provide the FASTQ files")
    pathlib.Path(bg_dir).mkdir(parents=True, exist_ok=True)
    cmd = ["fastq-dump", "--split-files", "--skip-technical", "--gzip",
           "-X", str(spots), "-O", str(bg_dir), str(sra)]
    print("[extract]", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return r1, r2


def build_mix(mix_id, bg_r1, bg_r2, spike_dir, pairs, expected, args):
    cmd = [sys.executable, str(GENERATOR),
           "--bg-r1", str(bg_r1), "--bg-r2", str(bg_r2),
           "--spike-pairs", str(pairs),
           "--bg-pairs", str(args.bg_pairs),
           "--read-len", str(args.read_len),
           "--insert-size", str(args.insert_size),
           "--seed", str(args.seed),
           "--mix-id", mix_id,
           "--outdir", args.outdir,
           "--expected-call", expected]
    if spike_dir:
        # The generator needs at least one reference even for a pure background
        # control; with --spike-pairs 0 nothing is inserted.
        cmd += ["--spike-dir", str(spike_dir)]
    else:
        cmd += ["--spike-dir", str(pathlib.Path(args.spike_root) / "cchfv")]
    print("[mix]", mix_id, "->", args.outdir)
    subprocess.run(cmd, check=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bg-dir", required=True,
                    help="directory with <accession>_1/_2.fastq.gz backgrounds")
    ap.add_argument("--sra-dir", default=None,
                    help="optional directory with <accession>.sra (used when the FASTQ is missing)")
    ap.add_argument("--spike-root", required=True,
                    help="directory with cchfv/ htnv/ seov/ sub-directories of *.fasta")
    ap.add_argument("--outdir", default="derived/t9_mixes")
    ap.add_argument("--bg-pairs", type=int, default=100000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--insert-size", type=int, default=300)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--spots", type=int, default=300000,
                    help="spots to extract per SRA run when converting (default 300000)")
    ap.add_argument("--only", default=None,
                    help="comma-separated mix ids to build (default: all)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the commands without writing anything")
    args = ap.parse_args(argv)

    wanted = {m.strip() for m in args.only.split(",")} if args.only else None
    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    rows = []
    for mix_id, accession, spike, pairs, expected in MIXES:
        if wanted and mix_id not in wanted:
            continue
        spike_dir = pathlib.Path(args.spike_root) / spike if spike else None
        if spike_dir and not spike_dir.is_dir() and not args.dry_run:
            raise SystemExit(f"missing spike directory: {spike_dir}")
        if args.dry_run:
            print(f"[dry-run] {mix_id}: {accession} + {pairs} pairs of "
                  f"{spike or 'no spike'} -> {args.outdir}")
        else:
            bg_r1, bg_r2 = extract_spots(accession, args.sra_dir, args.bg_dir,
                                         args.spots)
            build_mix(mix_id, bg_r1, bg_r2, spike_dir, pairs, expected, args)
        rows.append({
            "mix_id": mix_id,
            "background_library_id": accession,
            "spike_acronym": spike.upper(),
            "target_spike_reads": pairs,
            "total_pairs": args.bg_pairs + pairs,
            "expected_call": expected,
            "seed": args.seed,
            "read_len": args.read_len,
            "insert_size": args.insert_size,
            "outdir": args.outdir,
        })

    summary = pathlib.Path(args.outdir) / "t9_matrix_mixes.csv"
    if not args.dry_run and rows:
        with open(summary, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"[ok] {len(rows)} mix(es) -> {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
