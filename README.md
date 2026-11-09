# hvFinder

[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)]()
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](pyproject.toml)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-brightgreen.svg)]()
[![](https://img.shields.io/badge/build-passing-brightgreen)]()

**hvFinder** — a two-stage mNGS pipeline for detecting highly pathogenic viruses,
with a desktop GUI for setting up, submitting and monitoring runs on a Linux
server.

Given paired-end reads it performs QC, host and rRNA removal, assembly, a
protein-level pre-screen against a curated high-consequence virus database, then a
nucleotide-level confirmation against **hvn** (a curated high pathogenic virus
nucleotide database — not NCBI `nt`), followed by a comprehensive nr DIAMOND search
on the *confirmed* contigs only. The two-tier design keeps DIAMOND's
sensitivity for divergent viruses while BLASTN controls specificity, which is what
prevents near-neighbour cross-reactivity from being reported.

- **Server side** (`server/`) — plain Python 3.9+, **standard library + PyYAML
  only**. No GUI or SSH dependency, so it runs on compute nodes as-is.
- **Client side** (`client/`) — PySide6 desktop app (Windows / macOS / Linux) with
  three tabs: **1. SSH & Files** (host, credentials, work directory, reference
  data check), **2. Pre-processing** and **3. Pathogen Detection** (all tool paths
  and parameters). It uploads the scripts and config, submits the job over SSH —
  directly or through Slurm — then polls it and streams the remote log back into
  the window, showing progress per step. Every run is checkpointed, so a run that
  was stopped or lost can be resumed from the GUI without redoing finished steps,
  and results are pulled back for local inspection once the job finishes.

---

## Pipeline overview

```
                     Tab 2 — pre-processing                    Tab 3 — detection
  R1/R2 FASTQ ──► fastp ──► bowtie2 ──► ribodetector ──► MEGAHIT ──► CD-HIT-EST
                  (QC)     (host-)      (rRNA)          (assembly)  (dedup)
                                                    │
                                                    ▼
                                        PRESCREEN (DIAMOND, hvp)
                                                    │
                                     prescreen_candidates.fasta
                                                    ▼
                               BLASTN (hvn) ──► nr DIAMOND (confirmed contigs)
                                                    │
                                        coverage (bowtie2/samtools/bedtools)
                                                    │
                                     detection_results.tsv + bigWig + summary
```

Both stages are checkpointed (`preprocess.checkpoints.json`,
`pathogen_detect.checkpoints.json`). Re-running the **same command** after an
interruption resumes from the last completed step; `--no-resume` starts over.

---

## Requirements

| | Server (analysis) | Client (desktop GUI) |
|---|---|---|
| Python | 3.9+ | 3.9+ |
| Python packages | PyYAML | PySide6, paramiko (cryptography for key-based auth) |
| External programs | 11 command-line tools — listed under *Installation* and checked by `tests/selfcheck.sh` | — |
| Reference data | hvn, hvp, nr DIAMOND database, bowtie2 host index, taxid list — see [Required data](#required-data) | — |
| Scheduler | optional: a job runs on the host directly, `sbatch` is only used when the GUI is set to Slurm mode | — |
| Network | SSH from the client, outbound not required | SSH to the server |

## Installation

### Server (Linux)

```bash
python -m pip install -e .          # installs PyYAML only
# or simply run the scripts directly, no install needed:
python server/mNGS_pip.py --help
python server/pathogen_detect.py --help
```

External tools must be on `PATH` (or given as absolute paths in the config):
`fastp`, `bowtie2`, `samtools`, `ribodetector_cpu`, `megahit`, `cdhit-est`,
`diamond`, `blastn`, `blastdbcmd`, `bedtools`, `bedGraphToBigWig`.

### Client (desktop GUI)

The pre-built standalone executables for Windows and Linux platforms are
available in `dist/` — no Python installation is needed for those. macOS users
either run from source as shown below or build the binary themselves.

Running from source:

```bash
python -m pip install -e ".[client]"      # PySide6 + paramiko
python client/gui.py
```

To build a standalone executable (PyInstaller cannot cross-compile, so build on
the target platform):

```powershell
powershell -ExecutionPolicy Bypass -File build\build_windows.ps1   # → dist\hvFinder.exe
```

```bash
bash build/build_unix.sh                  # macOS / Linux → dist/hvFinder
bash build/build_unix.sh --refresh-env    # update the conda env first, then build
```

Both scripts create (or reuse) the `hvFinder` conda environment from
`build/environment-win.yml` / `build/environment-unix.yml`.
See `build/README_BUILD.md` for details.

---

## Self-check

Run this **once on the server** after installation — it is the quickest way to see
what is still missing:

```bash
bash tests/selfcheck.sh                                   # check everything
bash tests/selfcheck.sh --stage preprocess --outdir /data/run
bash tests/selfcheck.sh --deep                            # also run blastdbcmd -info / diamond dbinfo
```

It runs in the foreground and is **read-only** (no job submitted, nothing written):
Python version and PyYAML, code syntax, every tool declared in `config2.yaml`, the
reference data (taxid list, hvn, hvp, bowtie2 host index, the nr DIAMOND
database), output-directory writability and free space, and whether `sbatch`
exists. Exit code `0` = ready (`[WARN]` items are optional), `1` = at least one
`[FAIL]`.

## Required data

The repository does **not** ship reference databases. See [`data/README.md`](data/README.md)
for the expected layout; the short version:

| Item | Purpose |
|---|---|
| `data/virus_taxids_high_pathogenic.txt` | target taxids, one per line (`#` comments allowed) |
| `data/hvn/hvn.fa` (+ BLAST DB) | **hvn** — high pathogenic virus nucleotide DB, used by BLASTN for confirmation |
| `data/hvp/hvp.dmnd` | **hvp** — high pathogenic virus protein DB, used by the DIAMOND pre-screen |
| `data/databases/nr.dmnd` | NCBI nr, for the comprehensive DIAMOND pass |
| `host_idx/host_idx.*` | bowtie2 host index |

Taxids can be looked up in NCBI Taxonomy (web search or the E-utilities API);
the file itself is plain text, one taxid per line.

---

## Configuration

Everything the pipeline reads lives in **one** file, `server/configs/config2.yaml`:

- each tool has a `bin` key — a **plain command name** by default, so nothing is
  pinned to a specific machine; set an absolute path only if the tool is not on
  `PATH`;
- database paths (`diamond.db_path`, `prescreen.db_path`, `bowtie2.params.index_path`,
  ...) are empty or relative placeholders and must be filled in;
- thresholds and DIAMOND memory parameters are documented inline.

The GUI offers the same knobs: **Tab 2** owns the pre-processing stage
(including the small pre-screen DB), **Tab 3** owns the detection stage
(BLASTN parameters and the nr DIAMOND DB/parameters). Saving from the GUI writes
both `client/settings.ini` and the server-side `config2.yaml`.

DIAMOND memory model: `-c 1` is fixed, `-b = max(2, ⌊(node_mem − min(30, 0.3·node_mem)) / 6⌋)`
(e.g. 231 GB → 33). Nothing is auto-detected — if `node_mem_gb` is 0, no `-b/-c`
is passed and DIAMOND keeps its own defaults.

---

## Credentials

`client/settings.ini` ships as a **blank template** — host, user, password and all
paths are empty, and every tool binary is a plain command name. Fill it in on your
own machine; if you keep it under version control, avoid committing real
credentials (`git update-index --assume-unchanged client/settings.ini`).

## Usage

### Command line (server)

```bash
# stage 1 — pre-processing, writes <out>/prescreen_candidates.fasta
python server/mNGS_pip.py -1 R1.fq.gz -2 R2.fq.gz -o out -t 64 \
    -c server/configs/config2.yaml

# stage 2 — detection (BLASTN vs hvn → nr DIAMOND → coverage → report)
python server/pathogen_detect.py --fq1 R1.fq.gz --fq2 R2.fq.gz \
    --megahit out --diamond out/nr_hits.tsv --output out/pathogen \
    --taxids data/virus_taxids_high_pathogenic.txt \
    --virus-index data/hvn/hvn.fa \
    --threads 64 -c server/configs/config2.yaml
```

### GUI

1. Fill in the SSH host / port / user / password and the remote working directory.
2. Pick or create a project, set R1/R2, choose `Direct Python` or `Slurm (srun)`.
3. **RUN PIPELINE** (Tab 2) → wait for `prescreen_candidates.fasta` →
   **RUN DETECTION** (Tab 3).
4. The log box reports `✅ … finished` or `❌ … failed (exit code N)`; the tab
   status shows `● Finished` even for runs started outside the GUI.

---

## Outputs

| File | Content |
|---|---|
| `preprocess/prescreen_candidates.fasta` | contigs kept by the protein pre-screen (input of stage 2) |
| `pathogen/detection_results.tsv` | per-contig evidence (23 columns, stable order) |
| `pathogen/detection_virus_summary.tsv` | per-virus call, margin and coverage |
| `pathogen/detection_thresholds.json` | every threshold and parameter actually used |
| `pathogen/contig_coverage.bw` | coverage track for genome browsers |
| `pathogen/detection_summary.txt` | human-readable verdict (`No pathogen detected.` when empty) |

Verdicts are `POSITIVE` / `INDETERMINATE` / `NEGATIVE`. `INDETERMINATE` means the
evidence passed every gate but the margin over the thresholds is below
`margin_positive`; it is not a weak positive and should not be counted as one.
CSV/TSV column order is treated as a contract: columns are appended, never
renamed or reordered.

---

## Repository layout

```
client/            PySide6 GUI: SSH, submit, poll, view results
server/            pipeline: mNGS_pip.py (stage 1), pathogen_detect.py (stage 2)
  configs/         config2.yaml (all tool paths & parameters), species_aliases.csv
tests/             server-side self-check (selfcheck.sh) — no test suite is shipped
build/             PyInstaller build scripts (Windows / macOS / Linux) + conda env
hvFinder.spec      PyInstaller spec (working-directory independent)
pyproject.toml     packaging metadata and optional extras (client / build)
data/              (placeholder) reference databases — see data/README.md
dist/              pre-built standalone client (Windows .exe, Linux binary)
```

---

## License

MIT.

## Citation

If you use hvFinder in published work, please cite the accompanying manuscript
(reference to be added upon publication).
