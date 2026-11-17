#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spike_in_silico.py
==================
C3 (in-silico cross-reactivity): simulate reads from near-neighbour viral references,
spike them into a negative background, and produce paired FASTQ for hvFinder.

Design points (reviewer Q1: no sensitivity floor / cross-reactivity):
    * the spike-in amount is controlled exactly in **read pairs** (1e2 / 1e3 / 1e4);
    * each mix changes one variable only (virus or amount); the seed is fixed;
    * every simulated read carries its accession, e.g. `@spike|NC_004159.1|17`;
      the count can be checked with `grep -c '^@spike|NC_004159.1'`;
    * the expected call is `NEGATIVE` (these viruses are outside the target taxids):
      a POSITIVE call would indicate cross-reactivity or a false positive.

Usage:
    # one virus, one spike-in amount
    python spike_in_silico.py --bg-r1 bg_R1.fastq.gz --bg-r2 bg_R2.fastq.gz \
        --spike-fasta derived/spike_sequences/dugv_NC_004159.1.fasta \
        --spike-pairs 1000 --mix-id C3_DUGV_1e3 --outdir derived/mixes

    # spike every reference in a directory (reads split in proportion)
    python spike_in_silico.py ... --spike-dir derived/spike_sequences --spike-pairs 1000

    # gradient (loop, or call several times)
    for n in 100 1000 10000; do ... --spike-pairs $n --mix-id C3_DUGV_$n ; done

Output:
    <outdir>/<mix_id>_R1.fastq.gz
    <outdir>/<mix_id>_R2.fastq.gz
    <outdir>/<mix_id>.report.json     counts / seed / reference versions / md5 / expected call

Notes:
    * run this script **locally**; the FASTQ files are uploaded to the server afterwards;
    * background reads are reservoir-sampled (first 5,000,000 pairs by default, for speed);
      sampled reads are shuffled; simulated reads match the background read length.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import json
import os
import random
import sys
from datetime import datetime

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def open_maybe_gz(path: str):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b" or path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def read_fastq_pairs(path: str):
    """Yield (name, seq, qual) tuples; four lines per record; paired files must match."""
    with open_maybe_gz(path) as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            yield head[1:].rstrip("\n").split()[0], seq, qual


def reservoir_pairs(r1: str, r2: str, want: int, max_scan: int, rng: random.Random):
    """Reservoir-sample paired FASTQ; returns [(n, s1, q1, s2, q2), ...]."""
    if want <= 0:
        want = 10 ** 12                      # 0 = take everything
    pool: list[tuple] = []
    seen = 0
    for (n1, s1, q1), (n2, s2, q2) in zip(read_fastq_pairs(r1), read_fastq_pairs(r2)):
        seen += 1
        rec = (n1, s1, q1, s2, q2)
        if len(pool) < want:
            pool.append(rec)
        else:
            j = rng.randrange(seen)
            if j < want:
                pool[j] = rec
        if max_scan and seen >= max_scan:
            break
    rng.shuffle(pool)
    return pool, seen


def read_fasta(path: str):
    name, chunks = None, []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line.strip().upper())
    if name:
        yield name, "".join(chunks)


def collect_spikes(spike_files: list[str]) -> list[tuple[str, str]]:
    seqs = []
    for path in spike_files:
        for name, seq in read_fasta(path):
            if seq:
                seqs.append((name, seq))
    return seqs


def simulate_reads(spikes: list[tuple[str, str]], n_pairs: int, read_len: int,
                   insert_size: int, rng: random.Random, quality: str,
                   error_rate: float) -> list[tuple]:
    """Take random fragments from the references; yield (name, s1, q1, s2, q2)."""
    if not spikes:
        return []
    out = []
    # rotate over the references so that every fragment gets a chance
    idx = 0
    for i in range(n_pairs):
        acc, seq = spikes[idx % len(spikes)]
        idx += 1
        ins = min(insert_size, len(seq))
        if ins <= read_len:
            ins = min(len(seq), read_len * 2)
        start = rng.randrange(0, max(1, len(seq) - ins + 1))
        frag = seq[start:start + ins]
        if len(frag) <= read_len:
            r1 = frag
            r2 = revcomp(frag)
        else:
            r1 = frag[:read_len]
            r2 = revcomp(frag[-read_len:])
        if error_rate > 0:
            r1 = mutate(r1, error_rate, rng)
            r2 = mutate(r2, error_rate, rng)
        name = f"spike|{acc}|{i}"
        out.append((name, r1, quality * len(r1), r2, quality * len(r2)))
    return out


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    bases = "ACGT"
    out = []
    for b in seq:
        if rng.random() < rate:
            out.append(rng.choice([x for x in bases if x != b]))
        else:
            out.append(b)
    return "".join(out)


def write_pairs(path: str, recs: list[tuple], mate: int, gz: bool):
    opener = gzip.open if gz else open
    kw = {"compresslevel": 6} if gz else {}
    with opener(path, "wt", **kw) as fh:
        for rec in recs:
            name = rec[0]
            s = rec[1] if mate == 1 else rec[3]
            q = rec[2] if mate == 1 else rec[4]
            fh.write(f"@{name}/{mate}\n{s}\n+\n{q}\n")


def md5_of(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="C3: in-silico spike-in producing paired FASTQ")
    ap.add_argument("--bg-r1", required=True)
    ap.add_argument("--bg-r2", required=True)
    ap.add_argument("--spike-fasta", action="append", default=[])
    ap.add_argument("--spike-dir", default=None, help="spike every *.fasta in this directory")
    ap.add_argument("--spike-pairs", type=int, default=1000)
    ap.add_argument("--bg-pairs", type=int, default=100000,
                    help="how many pairs to draw from the background (0 = all)")
    ap.add_argument("--max-scan-pairs", type=int, default=5_000_000,
                    help="maximum pairs scanned by the reservoir (0 = unlimited)")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--insert-size", type=int, default=300)
    ap.add_argument("--error-rate", type=float, default=0.0)
    ap.add_argument("--quality", default="I")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--mix-id", required=True)
    ap.add_argument("--outdir", default="derived/mixes")
    ap.add_argument("--no-gzip", action="store_true")
    ap.add_argument("--expected-call", default="NEGATIVE")
    args = ap.parse_args(argv)

    spike_files = list(args.spike_fasta)
    if args.spike_dir:
        spike_files += sorted(glob.glob(os.path.join(args.spike_dir, "*.fasta")))
    if not spike_files:
        print("[error] --spike-fasta or --spike-dir is required", file=sys.stderr)
        return 2
    spikes = collect_spikes(spike_files)
    if not spikes:
        print("[error] reference sequences are empty", file=sys.stderr)
        return 2

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    bg, scanned = reservoir_pairs(args.bg_r1, args.bg_r2, args.bg_pairs,
                                  args.max_scan_pairs, rng)
    if not bg:
        print("[error] background FASTQ is empty or unreadable", file=sys.stderr)
        return 2

    spikes_recs = simulate_reads(spikes, args.spike_pairs, args.read_len,
                                 args.insert_size, rng, args.quality,
                                 args.error_rate)
    recs = bg + spikes_recs
    rng.shuffle(recs)

    gz = not args.no_gzip
    suffix = ".fastq.gz" if gz else ".fastq"
    r1 = os.path.join(args.outdir, f"{args.mix_id}_R1{suffix}")
    r2 = os.path.join(args.outdir, f"{args.mix_id}_R2{suffix}")
    write_pairs(r1, recs, 1, gz)
    write_pairs(r2, recs, 2, gz)

    report = {
        "mix_id": args.mix_id,
        "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "background": {"r1": os.path.abspath(args.bg_r1),
                       "r2": os.path.abspath(args.bg_r2),
                       "pairs_scanned": scanned,
                       "pairs_used": len(bg)},
        "spike": {"pairs": len(spikes_recs),
                  "read_len": args.read_len,
                  "insert_size": args.insert_size,
                  "error_rate": args.error_rate,
                  "references": [name for name, _ in spikes],
                  "files": [os.path.abspath(p) for p in spike_files]},
        "seed": args.seed,
        "total_pairs": len(recs),
        "expected_call": args.expected_call,
        "outputs": {"r1": os.path.abspath(r1), "r2": os.path.abspath(r2),
                    "r1_md5": md5_of(r1), "r2_md5": md5_of(r2)},
    }
    rep = os.path.join(args.outdir, f"{args.mix_id}.report.json")
    with open(rep, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    print(f"[ok] {args.mix_id}: bg={len(bg)} pairs (scanned {scanned}), "
          f"spike={len(spikes_recs)} pairs, total={len(recs)}")
    print(f"     {r1}")
    print(f"     {r2}")
    print(f"     {rep}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
