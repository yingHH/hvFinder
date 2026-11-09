# -*- coding: utf-8 -*-
"""
Author: Ying Huang
Date: 2026-01-27
Description: Pathogen Detection Pipeline
- Filter diamond results by taxid list
- Extract corresponding contigs from megahit results
- Run bowtie2 alignment to target virus index
- Generate detection report
"""

import os
import sys
import re
import time
import argparse
import subprocess
import shutil
import json
from collections import defaultdict
from datetime import datetime

try:
    import yaml
except ImportError:
    yaml = None

try:
    from .logger import setup_logger
except ImportError:
    try:
        from logger import setup_logger
    except ImportError:
        import logging

        def setup_logger(name, log_file):
            l = logging.getLogger(name)
            l.setLevel(logging.INFO)
            fh = logging.FileHandler(log_file)
            sh = logging.StreamHandler()
            fmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            fh.setFormatter(fmt)
            sh.setFormatter(fmt)
            l.addHandler(fh)
            l.addHandler(sh)
            return l


def load_checkpoints(checkpoint_file):
    """
    Load the checkpoint; an unreadable file is treated as absent.
    """
    if not os.path.exists(checkpoint_file):
        return {"completed_steps": []}
    try:
        with open(checkpoint_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        if not isinstance(data, dict):
            raise ValueError("unexpected checkpoint structure")
        data.setdefault("completed_steps", [])
        return data
    except (OSError, ValueError) as exc:
        try:
            os.replace(checkpoint_file, checkpoint_file + ".corrupt")
        except OSError:
            pass
        logger = logging.getLogger("PathogenDetect")
        logger.warning(
            f"Checkpoint file is unreadable or corrupt ({exc}); "
            "ignoring it and recomputing from the first step.")
        return {"completed_steps": []}

def save_checkpoints(checkpoints, checkpoint_file):
    """
    Persist the checkpoint.
    """
    with open(checkpoint_file, 'w', encoding='utf-8') as f:
        json.dump(checkpoints, f, indent=2)

def should_skip_step(step_name, checkpoints, output_files):
    """
    Whether a step can be skipped.
    """
    if step_name not in checkpoints.get("completed_steps", []):
        return False
    for f in output_files:
        if not os.path.exists(f):
            return False
    return True


def load_tool_config(config_path):
    """
    Read config2.yaml (the only source of tool paths); {} when unavailable.
    """
    if not config_path or yaml is None or not os.path.exists(config_path):
        return {}
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def _config_bin(cfg, section_name):
    """
    Bin of one section (accepts `bin` and `*bin` spellings); '' when absent.
    """
    if not cfg:
        return ''
    section = cfg.get(section_name)
    if not isinstance(section, dict):
        return ''
    value = section.get('bin') or section.get(f'{section_name}_bin') or ''
    return str(value).strip()


BIN_SPECS = {
    'bowtie2_bin': ('bowtie2', 'bowtie2'),
    'samtools_bin': ('samtools', 'samtools'),
    'bedtools_bin': ('bedtools', 'bedtools'),
    'blastn_bin': ('blastn', 'blastn'),
    'bedgraph_to_bigwig_bin': ('bedgraph_to_bigwig', 'bedGraphToBigWig'),
    'blastdbcmd_bin': ('blastdbcmd', 'blastdbcmd'),
}


def resolve_tool_bins(config_data, cli_values):
    """
    Tool path priority: explicit CLI value > config section bin > built-in name.
    """
    resolved = {}
    for attr, (section, fallback) in BIN_SPECS.items():
        cli_value = (cli_values or {}).get(attr)
        resolved[attr] = cli_value or _config_bin(config_data, section) or fallback
    return resolved


def resolve_db_path(cli_value, config_data, key):
    """
    Database path priority: explicit CLI value > `databases.<key>` in
    config2.yaml. None when neither is given - no path is hard-coded.
    """
    if cli_value:
        return os.path.expanduser(str(cli_value))
    dbs = (config_data or {}).get('databases')
    if not isinstance(dbs, dict):
        return None
    value = dbs.get(key)
    if not value:
        return None
    return os.path.expanduser(str(value))


def pre_check(args, logger):
    """
    Pre-flight check of the environment and of the input files.
    """
    logger.info("--- [Phase: Pre-Check] Validating Environment ---")
    validation_failed = False
    
    input_files = {
        'diamond': args.diamond,
        'candidates': args.candidates,
        'fq1': args.fq1,
        'fq2': args.fq2,
        'taxids': args.taxids,
        'megahit': args.megahit,
    }
    for name, path in input_files.items():
        if not path:
            continue
        if not os.path.exists(path):
            logger.error(f"Input file missing: [{name}] -> {path}")
            validation_failed = True
        else:
            logger.info(f"Verified input: {name:15} at {path}")
    
    required_tools = {
        'bowtie2': args.bowtie2_bin,
        'bowtie2-build': f"{args.bowtie2_bin}-build",
        'samtools': args.samtools_bin,
        'bedtools': args.bedtools_bin,
        'bedGraphToBigWig': args.bedgraph_to_bigwig_bin,
        'blastn': args.blastn_bin,
    }
    for tool, path in required_tools.items():
        expanded_path = os.path.expanduser(path)
        if not (shutil.which(expanded_path) or os.path.exists(expanded_path)):
            logger.error(f"Software path invalid: [{tool}] -> {path}")
            validation_failed = True
        else:
            logger.info(f"Verified tool: {tool:20} at {expanded_path}")
    
    virus_db = os.path.expanduser(args.virus_index)
    required_index_files = ['.nhr', '.nin', '.nsq']
    for ext in required_index_files:
        index_file = f"{virus_db}{ext}"
        if not os.path.exists(index_file):
            logger.error(f"BLASTN index missing: {index_file}")
            validation_failed = True
    if not validation_failed:
        logger.info(f"Verified BLASTN database: {virus_db}")
    
    if validation_failed:
        logger.critical("Pre-check failed. Please rectify tool paths or input files.")
        sys.exit(1)
    
    logger.info("--- [Phase: Pre-Check] All systems GO ---\n")


def load_taxids(taxid_file):
    """Load taxid list from file.
    
    Support formats:
    - Single taxid per line: "544571"
    - Taxid with name: "544571  California mosquito pool virus"
    - Comments starting with #
    
    Returns:
        taxids: set of taxid strings
        taxid_names: dict mapping taxid to virus name
    """
    taxids = set()
    taxid_names = {}
    
    with open(taxid_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Split by whitespace (tab or space)
                parts = line.split(None, 1)  # Split into max 2 parts
                taxid = parts[0]
                name = parts[1] if len(parts) > 1 else ''
                
                taxids.add(taxid)
                if name:
                    taxid_names[taxid] = name
    
    return taxids, taxid_names


# ---------------------------------------------------------------------------
# Species composition table and annotation
# ---------------------------------------------------------------------------
DEFAULT_SPECIES_ALIASES = {
    "zaire ebolavirus": "EBOV", "ebola virus": "EBOV", "zaire virus": "EBOV",
    "sudan ebolavirus": "SUDV", "sudan virus": "SUDV",
    "bundibugyo ebolavirus": "BDBV", "bundibugyo virus": "BDBV",
    "tai forest ebolavirus": "TAFV", "taï forest ebolavirus": "TAFV",
    "reston ebolavirus": "RESTV", "reston virus": "RESTV",
    "bombali ebolavirus": "BOMV", "bombali virus": "BOMV",
    "marburg virus": "MARV", "marburg marburgvirus": "MARV",
    "hendra virus": "HeV", "hendra henipavirus": "HeV",
    "nipah virus": "NiV", "nipah henipavirus": "NiV",
    "lassa virus": "LASV", "lassa mammarenavirus": "LASV",
    "machupo virus": "MACV", "machupo mammarenavirus": "MACV",
    "hantaan virus": "HTNV", "hantaan orthohantavirus": "HTNV",
    "seoul virus": "SEOV", "seoul orthohantavirus": "SEOV",
    "crimean congo hemorrhagic fever virus": "CCHFV",
    "crimean congo haemorrhagic fever virus": "CCHFV",
    "congo hemorrhagic fever virus": "CCHFV",
    "crimean-congo hemorrhagic fever virus": "CCHFV",
    "monkeypox virus": "MPXV", "mpox virus": "MPXV",
    "variola virus": "VARV", "smallpox virus": "VARV",
}


def _norm_name(text):
    """
    Lowercase and collapse non-alphanumerics to single spaces, for matching.
    """
    return re.sub(r'[^a-z0-9]+', ' ', str(text).lower()).strip()


def load_species_table(path):
    """
    Read the species table (Family, Genus, Species, Acronym, sequence count).

    Conventions: `#` starts a comment line; the first row is the header; empty
    Family/Genus cells inherit the row above; a thousands separator may split a
    number across columns (e.g. 8,749).
    """
    rows = []
    family = genus = ""
    with open(path, 'r', encoding='utf-8-sig') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 4:
                continue
            if parts[0].lower() == 'family':
                continue
            fam, gen, sp, acr = parts[0], parts[1], parts[2], parts[3]
            if fam:
                family = fam
            if gen:
                genus = gen
            tail = [p for p in parts[4:] if p != '']
            n_seq = 0
            if tail:
                try:
                    n_seq = int(''.join(tail).replace(',', ''))
                except ValueError:
                    n_seq = 0
            rows.append({
                'family': family, 'genus': genus,
                'species': sp, 'acronym': acr, 'n_sequences': n_seq,
            })
    return rows


def build_species_index(rows, alias_path=None):
    """
    Build the normalised-name -> species-row index (scientific name, acronym,
    synonyms).
    """
    alias_map = {_norm_name(k): v for k, v in DEFAULT_SPECIES_ALIASES.items()}
    if alias_path and os.path.exists(alias_path):
        with open(alias_path, 'r', encoding='utf-8-sig') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = [p.strip() for p in line.split(',')]
                if len(parts) < 2 or parts[0].lower() == 'alias':
                    continue
                alias_map[_norm_name(parts[0])] = parts[1]

    by_acronym = {r['acronym']: r for r in rows}
    index = {}
    for r in rows:
        index[_norm_name(r['species'])] = r
        index[_norm_name(r['acronym'])] = r
    for alias_norm, acronym in alias_map.items():
        row = by_acronym.get(acronym)
        if row:
            index.setdefault(alias_norm, row)
    return index


def match_species(name, index, rows):
    """
    Match any name to a species row; None when nothing matches.
    """
    if not name:
        return None
    norm = _norm_name(name)
    if not norm:
        return None
    hit = index.get(norm)
    if isinstance(hit, dict):
        return hit

    tokens = set(norm.split())
    for key, row in index.items():
        if not isinstance(row, dict):
            continue
        kt = set(key.split())
        if kt and kt <= tokens:
            return row
    for key, row in index.items():
        if isinstance(row, dict) and len(key) >= 4 and ' ' not in key and key in norm:
            return row
    return None


def count_db_sequences(db_path, blastdbcmd_bin='blastdbcmd'):
    """
    Count sequences in a BLAST database with blastdbcmd; None when unavailable.
    """
    bin_path = os.path.expanduser(blastdbcmd_bin)
    if not (shutil.which(bin_path) or os.path.exists(bin_path)):
        return None
    cmd = f"{bin_path} -db {os.path.expanduser(db_path)} -entry all -outfmt '%i'"
    result = subprocess.run(cmd, shell=True, capture_output=True, encoding='utf-8', errors='replace')
    if result.returncode != 0:
        return None
    return len([line for line in result.stdout.splitlines() if line.strip()])


def filter_diamond(diamond_tsv, taxids, pident, evalue, logger=None):
    """Filter diamond results by taxid and parameters
    
    Args:
        diamond_tsv: Path to diamond TSV file
        taxids: Set of target taxids
        pident: Minimum identity threshold
        evalue: Maximum e-value threshold
        logger: Logger instance for output
    
    Returns:
        filtered_contigs: Dict of filtered contig information
    """
    filtered_contigs = defaultdict(dict)
    
    # Statistics for pre-filter summary
    total_hits = 0
    all_pidents = []
    all_qcovs = []
    all_evalues = []
    
    with open(diamond_tsv, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 15:
                continue
            
            qseqid = fields[0]  # contig ID
            sseqid = fields[1]  # reference ID
            qcovhsp = float(fields[2]) if fields[2] else 0  # qcovhsp is 3rd column in config2.yaml outfmt
            pident_val = float(fields[3]) if fields[3] else 0  # pident is 4th column
            evalue_val = float(fields[11]) if fields[11] else 1  # evalue is 12th column
            staxids = fields[12]  # taxonomy ID is 13th column
            sscinames = fields[13]  # scientific name is 14th column
            
            # Parse multiple taxids if present
            contig_taxids = staxids.split(';') if staxids else []

            # Check if any taxid matches our target list
            if not (set(contig_taxids) & taxids):
                continue
            
            total_hits += 1
            all_pidents.append(pident_val)
            all_qcovs.append(qcovhsp)
            all_evalues.append(evalue_val)
            
            if pident_val < pident:
                continue
            if evalue_val > evalue:
                continue
            
            # Store best match for each contig
            if qseqid not in filtered_contigs:
                filtered_contigs[qseqid] = {
                    'sseqid': sseqid,
                    'sscinames': sscinames,
                    'staxids': staxids,
                    'pident': pident_val,
                    'qcov': qcovhsp,
                    'evalue': evalue_val
                }
            else:
                # Keep the best match (highest pident)
                if pident_val > filtered_contigs[qseqid]['pident']:
                    filtered_contigs[qseqid] = {
                        'sseqid': sseqid,
                        'sscinames': sscinames,
                        'staxids': staxids,
                        'pident': pident_val,
                        'qcov': qcovhsp,
                        'evalue': evalue_val
                    }
    
    # Log pre-filter statistics
    if logger and total_hits > 0:
        logger.info("-" * 50)
        logger.info("Pre-filter statistics (taxid-matched hits only):")
        logger.info(f"  Total hits before filtering: {total_hits}")
        logger.info(f"  Identity: min={min(all_pidents):.1f}%, max={max(all_pidents):.1f}%, avg={sum(all_pidents)/len(all_pidents):.1f}%")
        logger.info(f"  Coverage: min={min(all_qcovs):.1f}%, max={max(all_qcovs):.1f}%, avg={sum(all_qcovs)/len(all_qcovs):.1f}%")
        logger.info(f"  E-value: min={min(all_evalues):.2e}, max={max(all_evalues):.2e}, median={sorted(all_evalues)[len(all_evalues)//2]:.2e}")
        logger.info("-" * 50)
    
    return filtered_contigs


def diagnose_diamond_filter(diamond_tsv, taxids, pident, evalue, logger=None):
    """
    Explain how many contigs the DIAMOND pre-filter dropped and how many a
    relaxed threshold would recover.
    """
    best = {}      # contig -> [max pident, min evalue]
    hits = 0
    if os.path.exists(diamond_tsv):
        with open(diamond_tsv, 'r', encoding='utf-8') as f:
            for line in f:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 15:
                    continue
                try:
                    pident_val = float(fields[3]) if fields[3] else 0.0
                    evalue_val = float(fields[11]) if fields[11] else 1.0
                except ValueError:
                    continue
                contig_taxids = set(fields[12].split(';')) if fields[12] else set()
                if not (contig_taxids & taxids):
                    continue
                hits += 1
                cur = best.get(fields[0])
                if cur is None:
                    best[fields[0]] = [pident_val, evalue_val]
                else:
                    cur[0] = max(cur[0], pident_val)
                    cur[1] = min(cur[1], evalue_val)

    n_contigs = len(best)
    kept = sum(1 for p, e in best.values() if p >= pident and e <= evalue)
    pident_sweep = {t: sum(1 for p, e in best.values()
                           if p >= t and e <= evalue)
                    for t in (0, 30, 40, 50, 60, 70, 80, 90)}
    evalue_sweep = {t: sum(1 for p, e in best.values()
                           if p >= pident and e <= t)
                    for t in (1e-2, 1e-5, 1e-10, 1e-20)}

    stats = {
        'hits': hits, 'contigs': n_contigs, 'kept': kept,
        'dropped': n_contigs - kept,
        'pident_sweep': pident_sweep, 'evalue_sweep': evalue_sweep,
    }

    if logger:
        logger.info("DIAMOND prefilter: %d taxid-matched hit(s) on %d contig(s); "
                    "%d contig(s) kept at pident>=%g and evalue<=%g",
                    hits, n_contigs, kept, pident, evalue)
        if kept < n_contigs:
            logger.warning(
                "  DIAMOND prefilter dropped %d/%d contig(s). "
                "Contigs kept vs. pident threshold: %s",
                n_contigs - kept, n_contigs,
                ", ".join(f">={t}:{c}" for t, c in pident_sweep.items()))
            logger.warning(
                "  Contigs kept vs. evalue threshold: %s",
                ", ".join(f"<={t:g}:{c}" for t, c in evalue_sweep.items()))
            logger.warning(
                "  The archived run used pident>=50 and evalue<=1e-02; "
                "BLASTN (identity/coverage) does the real confirmation, so a "
                "looser prefilter is usually safe.")

    return stats


def run_blastn(query_fa, db_path, out_file, threads, evalue, blastn_bin, logger):
    """Run BLASTN alignment (contigs vs virus reference)"""
    blastn_path = os.path.expanduser(blastn_bin)
    cmd = f"""
    {blastn_path} -query {query_fa} \
           -db {db_path} \
           -out {out_file} \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs slen" \
           -evalue {evalue} \
           -num_threads {threads}
    """
    logger.info("Running BLASTN alignment...")
    result = subprocess.run(cmd, shell=True, capture_output=True, encoding='utf-8', errors='replace')
    if result.returncode != 0:
        logger.error(f"BLASTN failed: {result.stderr}")
        return False
    logger.info(f"BLASTN completed. Results saved to {out_file}")
    return True


def filter_blastn_results(input_file, output_file, identity, coverage, logger):
    """Filter BLASTN results by identity and coverage"""
    count = 0
    with open(input_file, 'r', encoding='utf-8') as f_in, open(output_file, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue
            pident = float(fields[2])
            qcov = float(fields[12])
            if pident >= identity and qcov >= coverage:
                f_out.write(line)
                count += 1
    if logger:
        logger.info(f"Filtered BLASTN results: {count} hits saved to {output_file}")
    return count


def diagnose_blastn_filter(blastn_raw, identity, coverage, logger):
    """
    Explain why BLASTN kept nothing: hit distribution and pass counts per
    criterion.
    """
    total = pass_id = pass_cov = pass_both = 0
    p_min = q_min = None
    p_max = q_max = None
    if os.path.exists(blastn_raw):
        with open(blastn_raw, 'r', encoding='utf-8') as f:
            for line in f:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 13:
                    continue
                try:
                    pident = float(fields[2])
                    qcov = float(fields[12])
                except ValueError:
                    continue
                total += 1
                p_min = pident if p_min is None else min(p_min, pident)
                p_max = pident if p_max is None else max(p_max, pident)
                q_min = qcov if q_min is None else min(q_min, qcov)
                q_max = qcov if q_max is None else max(q_max, qcov)
                if pident >= identity:
                    pass_id += 1
                if qcov >= coverage:
                    pass_cov += 1
                if pident >= identity and qcov >= coverage:
                    pass_both += 1

    stats = {'total': total, 'pass_pident': pass_id, 'pass_qcov': pass_cov,
             'pass_both': pass_both, 'min_pident': p_min, 'max_pident': p_max,
             'min_qcov': q_min, 'max_qcov': q_max}

    if logger:
        if total == 0:
            logger.warning("BLASTN raw results are empty: none of the extracted "
                           "contigs hits the viral database")
            return stats
        logger.warning(f"0 hits left after BLASTN filtering (raw hits: {total})")
        logger.warning(f"  pident: threshold {identity}%  observed {p_min:.1f}-{p_max:.1f}%  "
                       f"pass {pass_id}/{total}")
        logger.warning(f"  qcovs : threshold {coverage}%  observed {q_min:.1f}-{q_max:.1f}%  "
                       f"pass {pass_cov}/{total}")
        blockers = []
        if pass_id == 0:
            blockers.append(f"--blastn-identity {identity}")
        if pass_cov == 0:
            blockers.append(f"--blastn-coverage {coverage}")
        if blockers:
            logger.warning("  blocker(s): " + " , ".join(blockers))
    return stats


def calculate_contig_coverage_with_bigwig(fq1, fq2, contigs_fa, output_bam, 
                                          output_summary, output_bw, 
                                          threads, bowtie2_bin, samtools_bin, 
                                          bedtools_bin, bedgraph_to_bigwig_bin, logger):
    """
    Map reads to contigs and calculate coverage
    Output: BAM, summary TSV, BigWig
    """
    bowtie2_path = os.path.expanduser(bowtie2_bin)
    samtools_path = os.path.expanduser(samtools_bin)
    bedtools_path = os.path.expanduser(bedtools_bin)
    bedgraph_to_bigwig_path = os.path.expanduser(bedgraph_to_bigwig_bin)
    
    logger.info("Building contigs index...")
    cmd_index = f"{bowtie2_path}-build {contigs_fa} {output_bam}.contigs_idx"
    result = subprocess.run(cmd_index, shell=True, capture_output=True, encoding='utf-8', errors='replace')
    if result.returncode != 0:
        logger.error(f"Failed to build contigs index: {result.stderr}")
        return False
    
    logger.info("Mapping reads to contigs...")
    sam_file = output_bam.replace('.bam', '.sam')
    cmd_map = f"""
    {bowtie2_path} -x {output_bam}.contigs_idx \
        -1 {fq1} -2 {fq2} \
        -S {sam_file} \
        --threads {threads}
    """
    result = subprocess.run(cmd_map, shell=True, capture_output=True, encoding='utf-8', errors='replace')
    if result.returncode != 0:
        logger.error(f"Failed to map reads: {result.stderr}")
        return False
    
    logger.info("Converting SAM to BAM and sorting...")
    cmd_bam = f"{samtools_path} view -bS {sam_file} | {samtools_path} sort -o {output_bam}"
    result = subprocess.run(cmd_bam, shell=True, capture_output=True, encoding='utf-8', errors='replace')
    if result.returncode != 0:
        logger.error(f"Failed to create BAM: {result.stderr}")
        return False
    
    cmd_index = f"{samtools_path} index {output_bam}"
    subprocess.run(cmd_index, shell=True, capture_output=True, encoding='utf-8', errors='replace')
    
    logger.info("Calculating coverage summary...")
    idxstats_file = output_summary + ".idxstats.tsv"
    depth_file = output_summary + ".depth.tsv"

    cmd_idxstats = f"{samtools_path} idxstats {output_bam} > {idxstats_file}"
    result = subprocess.run(cmd_idxstats, shell=True, capture_output=True, encoding='utf-8', errors='replace')
    if result.returncode != 0:
        logger.error(f"Failed to compute idxstats: {result.stderr}")
        return False

    cmd_summary = f"""
    {samtools_path} depth -a {output_bam} | awk '{{ 
        sum[$1]+=$3; n[$1]++; 
        if($3>max[$1]) max[$1]=$3; 
        if(!($1 in min) || $3<min[$1]) min[$1]=$3;
        if($3>0) cov[$1]++
    }} END {{
        for(c in sum) printf "%s\\t%.4f\\t%d\\t%d\\t%d\\t%d\\n", c, sum[c]/n[c], min[c], max[c], cov[c], n[c]
    }}' > {depth_file}
    """
    result = subprocess.run(cmd_summary, shell=True, capture_output=True, encoding='utf-8', errors='replace')
    if result.returncode != 0:
        logger.error(f"Failed to calculate coverage: {result.stderr}")
        return False

    write_coverage_summary(idxstats_file, depth_file, output_summary, logger)
    for tmp in (idxstats_file, depth_file):
        if os.path.exists(tmp):
            os.remove(tmp)
    
    logger.info("Generating BigWig...")
    if shutil.which(bedgraph_to_bigwig_path) or os.path.exists(bedgraph_to_bigwig_path):
        cmd_bedgraph = f"{bedtools_path} genomecov -ibam {output_bam} -bg > {output_bam}.coverage.bedgraph"
        result = subprocess.run(cmd_bedgraph, shell=True, capture_output=True, encoding='utf-8', errors='replace')
        if result.returncode != 0:
            logger.error(f"Failed to create bedGraph: {result.stderr}")
            return False
        
        cmd_sizes = f"{samtools_path} faidx {contigs_fa} && cut -f1,2 {contigs_fa}.fai > {contigs_fa}.sizes"
        result = subprocess.run(cmd_sizes, shell=True, capture_output=True, encoding='utf-8', errors='replace')
        if result.returncode != 0:
            logger.error(f"Failed to create sizes file: {result.stderr}")
            return False
        
        cmd_bw = f"{bedgraph_to_bigwig_path} {output_bam}.coverage.bedgraph {contigs_fa}.sizes {output_bw}"
        result = subprocess.run(cmd_bw, shell=True, capture_output=True, encoding='utf-8', errors='replace')
        if result.returncode != 0:
            logger.error(f"Failed to create BigWig: {result.stderr}")
            return False
        
        if os.path.exists(f"{output_bam}.coverage.bedgraph"):
            os.remove(f"{output_bam}.coverage.bedgraph")
        if os.path.exists(f"{contigs_fa}.sizes"):
            os.remove(f"{contigs_fa}.sizes")
        logger.info(f"BigWig saved to {output_bw}")
    else:
        logger.warning(f"bedGraphToBigWig not found at {bedgraph_to_bigwig_path}. Skipping BigWig generation.")
    
    if os.path.exists(sam_file):
        os.remove(sam_file)
    
    logger.info(f"Coverage calculation completed. BAM: {output_bam}")
    return True






# ---------------------------------------------------------------------------
# Detection thresholds and decision layer
# ---------------------------------------------------------------------------
COVERAGE_HEADER = ["Contig_ID", "Contig_Length", "Reads_Mapped", "Mean_Depth",
                   "Min_Depth", "Max_Depth", "Covered_Positions", "Breadth_Pct"]

VIRUS_SUMMARY_HEADER = (
    "Virus_Seq\tCall\tN_Contigs_Blastn_Pass\tN_Genome_Regions\t"
    "Best_Contig_ID\tBest_Contig_Accession\tBest_Contig_Len\t"
    "Best_Contig_Mean_Depth\t"
    "Best_Contig_Breadth_Pct\tBest_Contig_Max_Depth\t"
    "Supporting_Reads_Mapped\tSupporting_Contig_Sum_Len\t"
    "Blastn_Span_Union_Bp\tBlastn_Span_Overlap_Ratio\t"
    "Genome_Breadth_Pct_Blastn_Approx\tSubject_Length\tMargin_Fold\t"
    "Decision_Rule\tFlags\tReason\tVirus_Species\tVirus_Acronym\t"
    "Virus_Family\tVirus_Genus\n")

THRESHOLD_PROFILES = {
    "P0_current_default": {
        "min_contig_len": 0, "min_reads_mapped": 0, "min_mean_depth": 0,
        "min_contig_breadth": 0, "min_genome_breadth": 0,
        "margin_positive": 1.0, "min_supporting_contigs": 1,
    },
    "P1_two_tier": {
        "min_contig_len": 300, "min_reads_mapped": 10, "min_mean_depth": 1,
        "min_contig_breadth": 50, "min_genome_breadth": 10,
        "margin_positive": 2.0, "min_supporting_contigs": 1,
    },
    "P2_stringent": {
        "min_contig_len": 1000, "min_reads_mapped": 30, "min_mean_depth": 3,
        "min_contig_breadth": 70, "min_genome_breadth": 20,
        "margin_positive": 2.0, "min_supporting_contigs": 2,
    },
}


def write_coverage_summary(idxstats_file, depth_file, out_file, logger=None):
    """
    Merge idxstats (length, reads) with depth (mean depth, breadth) into one
    table.

    Breadth_Pct = bases with depth >= 1 / contig length * 100. The old
    Coverage(%) column stored Mean_Depth by mistake; both are reported now.
    """
    lengths, mapped = {}, {}
    if os.path.exists(idxstats_file):
        with open(idxstats_file, 'r', encoding='utf-8') as f:
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) < 4 or p[0] == '*':
                    continue
                try:
                    lengths[p[0]] = int(p[1])
                    mapped[p[0]] = int(p[2])
                except ValueError:
                    continue

    stats = {}
    if os.path.exists(depth_file):
        with open(depth_file, 'r', encoding='utf-8') as f:
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) < 6:
                    continue
                try:
                    stats[p[0]] = {
                        'mean': float(p[1]), 'min': int(float(p[2])),
                        'max': int(float(p[3])), 'cov': int(float(p[4])),
                        'n': int(float(p[5])),
                    }
                except ValueError:
                    continue

    n_written = 0
    with open(out_file, 'w', encoding='utf-8') as out:
        out.write('\t'.join(COVERAGE_HEADER) + '\n')
        for cid in sorted(set(lengths) | set(stats)):
            s = stats.get(cid)
            length = lengths.get(cid) or (s['n'] if s else 0)
            if not s or not length:
                mean = minv = maxv = cov = 0
                breadth = 0.0
            else:
                mean, minv, maxv, cov = s['mean'], s['min'], s['max'], s['cov']
                breadth = cov / length * 100.0
            out.write(f"{cid}\t{length}\t{mapped.get(cid, 0)}\t{mean:.4f}\t"
                      f"{minv}\t{maxv}\t{cov}\t{breadth:.2f}\n")
            n_written += 1

    if logger:
        logger.info(f"Coverage summary written for {n_written} contigs "
                    f"(contig length, mapped reads, depth, breadth)")


def parse_coverage_summary(path):
    """
    Parse contig_coverage_summary.tsv (the old headerless 5-column layout is
    also accepted).
    """
    stats = {}
    if not os.path.exists(path):
        return stats

    with open(path, 'r', encoding='utf-8') as f:
        lines = [ln.rstrip('\n') for ln in f if ln.strip()]

    has_header = bool(lines) and lines[0].split('\t')[0].strip() == 'Contig_ID'
    if has_header:
        for line in lines[1:]:
            p = line.split('\t')
            if len(p) < 8:
                continue
            try:
                stats[p[0]] = {
                    'length': int(float(p[1])), 'reads': int(float(p[2])),
                    'mean_depth': float(p[3]), 'max_depth': int(float(p[5])),
                    'breadth': float(p[7]),
                }
            except ValueError:
                continue
    else:
        for line in lines:
            p = line.split('\t')
            if len(p) < 5:
                continue
            try:
                stats[p[0]] = {
                    'length': int(float(p[4])), 'reads': 0,
                    'mean_depth': float(p[1]), 'max_depth': int(float(p[3])),
                    'breadth': 0.0,
                }
            except ValueError:
                continue
    return stats


def parse_blastn_best_hits(blastn_filtered):
    """
    Best BLASTN hit per contig (highest bitscore), keeping the subject span for
    the genome-breadth estimate.
    """
    hits = defaultdict(list)
    with open(blastn_filtered, 'r', encoding='utf-8') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 13:
                continue
            try:
                hits[p[0]].append({
                    'sseqid': p[1],
                    'pident': float(p[2]),
                    'evalue': float(p[10]),
                    'bitscore': float(p[11]),
                    'qcov': float(p[12]),
                    'sstart': int(p[8]),
                    'send': int(p[9]),
                    'slen': int(p[13]) if len(p) > 13 and p[13] else None,
                })
            except ValueError:
                continue

    best = {}
    for cid, lst in hits.items():
        b = max(lst, key=lambda x: x['bitscore'])
        b['interval'] = (min(b['sstart'], b['send']), max(b['sstart'], b['send']))
        best[cid] = b
    return best


def merge_intervals(intervals):
    """
    Merge closed intervals, to count non-overlapping genomic regions hit.
    """
    ivs = sorted((int(s), int(e)) for s, e in intervals if s is not None and e is not None)
    if not ivs:
        return []
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def get_subject_lengths(db_path, subject_ids, blastdbcmd_bin='blastdbcmd', logger=None):
    """
    Reference lengths via blastdbcmd - fallback only; BLASTN's slen is the
    normal source. Returns {} (breadth becomes NA) when unavailable.
    """
    lengths = {}
    if not subject_ids:
        return lengths

    bin_path = os.path.expanduser(blastdbcmd_bin)
    if not (shutil.which(bin_path) or os.path.exists(bin_path)):
        if logger:
            logger.warning("blastdbcmd not found; genome breadth will be reported as NA")
        return lengths

    tmp_path = None
    try:
        import tempfile
        with tempfile.NamedTemporaryFile('w', suffix='.ids', delete=False) as tmp:
            tmp_path = tmp.name
            for sid in sorted(subject_ids):
                tmp.write(f"{sid}\n")

        cmd = (f"{bin_path} -db {os.path.expanduser(db_path)} "
               f"-entry_batch {tmp_path} -outfmt '%i\\t%l'")
        result = subprocess.run(cmd, shell=True, capture_output=True, encoding='utf-8', errors='replace')
        if result.returncode != 0:
            if logger:
                logger.warning(f"blastdbcmd failed: {result.stderr.strip()}")
            return lengths
        for line in result.stdout.splitlines():
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                try:
                    lengths[p[0]] = int(p[1])
                except ValueError:
                    continue
    except Exception as exc:
        if logger:
            logger.warning(f"Could not resolve subject lengths: {exc}")
    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)
    return lengths


def evaluate_detection(per_contig, subject_lengths, thr, logger=None):
    """
    Contig-length gate, then reads/depth/breadth thresholds and the margin call.

    Returns {subject_id: {call, margin, genome_breadth_pct, ...}} where call is
    POSITIVE / INDETERMINATE / NEGATIVE.
    """
    groups = defaultdict(list)
    for row in per_contig:
        groups[row['Best_Hit_Seq']].append(row)

    decisions = {}
    for subject, contigs in groups.items():
        kept = [c for c in contigs if c['Contig_Length'] >= thr['min_contig_len']]
        dropped = len(contigs) - len(kept)
        n_contigs = len(kept)
        total_len = sum(c['Contig_Length'] for c in kept)
        reads = sum(c['Reads_Mapped'] for c in kept)
        mean_depth = (sum(c['Mean_Depth'] * c['Contig_Length'] for c in kept) / total_len
                      if total_len else 0.0)
        max_depth = max((c['Max_Depth'] for c in kept), default=0)
        max_breadth = max((c['Breadth_Pct'] for c in kept), default=0.0)

        best = max(kept, key=lambda c: (c['Breadth_Pct'], c['Mean_Depth'],
                                        c['Contig_Length'])) if kept else None

        merged = merge_intervals([c['interval'] for c in kept if c.get('interval')])
        span = sum(e - s + 1 for s, e in merged)
        span_sum = sum(e - s + 1 for c in kept if c.get('interval')
                       for s, e in [c['interval']])
        overlap_ratio = (1.0 - span / span_sum) if span_sum else None
        ref_len = subject_lengths.get(subject)
        genome_breadth = (span / ref_len * 100.0) if ref_len else None

        metrics = []
        if thr['min_reads_mapped'] > 0:
            metrics.append(('reads_mapped', reads, thr['min_reads_mapped']))
        if thr['min_mean_depth'] > 0:
            metrics.append(('mean_depth', mean_depth, thr['min_mean_depth']))
        if thr['min_contig_breadth'] > 0:
            metrics.append(('contig_breadth', max_breadth, thr['min_contig_breadth']))
        if thr['min_genome_breadth'] > 0 and genome_breadth is not None:
            metrics.append(('genome_breadth', genome_breadth, thr['min_genome_breadth']))

        reason = ''
        if n_contigs == 0:
            call, margin = 'NEGATIVE', 0.0
            reason = f"no contig >= {thr['min_contig_len']} bp"
        elif n_contigs < thr['min_supporting_contigs']:
            call, margin = 'NEGATIVE', 0.0
            reason = (f"n_contigs {n_contigs} < min_supporting_contigs "
                      f"{thr['min_supporting_contigs']}")
        elif not metrics:
            call, margin = 'POSITIVE', None
            reason = 'no numeric threshold active (legacy behaviour)'
        else:
            ratios = [(name, (value / t) if t > 0 else float('inf'))
                      for name, value, t in metrics]
            margin_candidates = [r for name, r in ratios
                                 if name in ('reads_mapped', 'mean_depth')]
            margin = (min(margin_candidates) if margin_candidates
                      else min(r for _, r in ratios))
            below = [name for name, r in ratios if r < 1.0]
            if below:
                call = 'NEGATIVE'
                reason = 'below threshold: ' + ','.join(below)
            elif margin >= thr['margin_positive']:
                call = 'POSITIVE'
            else:
                call = 'INDETERMINATE'
                reason = f"margin {margin:.2f} < {thr['margin_positive']}"

        flags = []
        if margin is not None and 1.0 <= margin < thr['margin_positive']:
            flags.append('near_threshold')
        if n_contigs == 1:
            flags.append('single_contig')
        if thr['min_contig_breadth'] > 0 and 1.0 <= max_breadth < 1.5 * thr['min_contig_breadth']:
            flags.append('low_breadth')
        if merged and len(merged) < 2:
            flags.append('single_locus')
        if dropped:
            flags.append(f'dropped_{dropped}_short_contig')

        decisions[subject] = {
            'call': call,
            'margin': margin,
            'genome_breadth_pct': genome_breadth,
            'n_contigs': n_contigs,
            'n_genome_regions': len(merged),
            'total_contig_len': total_len,
            'reads_mapped': reads,
            'mean_depth': mean_depth,
            'max_depth': max_depth,
            'contig_breadth_max': max_breadth,
            'subject_length': ref_len,
            'flags': ';'.join(flags) if flags else 'NA',
            'reason': reason,
            'best_contig_id': best['Contig_ID'] if best else None,
            'best_contig_accession': best['Best_Hit_Seq'] if best else None,
            'best_contig_len': best['Contig_Length'] if best else None,
            'best_contig_mean_depth': best['Mean_Depth'] if best else None,
            'best_contig_breadth_pct': best['Breadth_Pct'] if best else None,
            'best_contig_max_depth': best['Max_Depth'] if best else None,
            'span_union_bp': span,
            'span_sum_bp': span_sum,
            'span_overlap_ratio': overlap_ratio,
        }

        if logger:
            margin_txt = 'NA' if margin is None else f"{margin:.2f}"
            logger.info(f"  {subject}: {call} (contigs={n_contigs}, reads={reads}, "
                        f"depth={mean_depth:.2f}, breadth={max_breadth:.1f}%, "
                        f"margin={margin_txt})")
            if reason:
                logger.info(f"    reason: {reason}")

    return decisions



# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
DIAMOND_GB_PER_BLOCK = 6.0
DIAMOND_MEM_RESERVE_GB = 30.0
DIAMOND_MAX_RESERVE_RATIO = 0.3
DIAMOND_MIN_BLOCKS = 2
DIAMOND_INDEX_CHUNKS = 1
DIAMOND_OUTFMT = ("6 qseqid sseqid qcovhsp pident qlen slen length qstart qend "
                  "sstart send evalue staxids sscinames sskingdoms")


def plan_diamond_blocks(node_mem_gb, reserve_gb=DIAMOND_MEM_RESERVE_GB,
                        gb_per_block=DIAMOND_GB_PER_BLOCK,
                        chunks=DIAMOND_INDEX_CHUNKS,
                        min_blocks=DIAMOND_MIN_BLOCKS):
    """
    Derive DIAMOND's ``-b`` from ``memory(GB) ~ gb_per_block * (b / c)``.
    """
    try:
        total = float(node_mem_gb or 0)
    except (TypeError, ValueError):
        return None, None
    if total <= 0:
        return None, None
    chunks = max(1, int(chunks or 1))
    reserve = min(float(reserve_gb), DIAMOND_MAX_RESERVE_RATIO * total)
    budget = total - reserve
    blocks = max(int(min_blocks), int(budget * chunks / float(gb_per_block)))
    return blocks, budget


def plan_diamond_mem(config_data, node_mem_gb, logger=None):
    """
    Build the -b/-c arguments from config2.yaml diamond.params and
    --node-mem-gb.
    """
    params = ((config_data or {}).get('diamond') or {}).get('params') or {}
    total = node_mem_gb if node_mem_gb is not None else params.get('node_mem_gb')
    try:
        total = float(total or 0)
    except (TypeError, ValueError):
        total = 0.0
    if total <= 0:
        return "", None
    chunks = int(params.get('index_chunks') or DIAMOND_INDEX_CHUNKS)
    reserve = float(params.get('mem_reserve_gb', DIAMOND_MEM_RESERVE_GB))
    gb_per_block = float(params.get('gb_per_block', DIAMOND_GB_PER_BLOCK))
    blocks, budget = plan_diamond_blocks(total, reserve, gb_per_block,
                                         chunks=chunks)
    manual = params.get('blocks') or 0
    try:
        manual = int(manual)
    except (TypeError, ValueError):
        manual = 0
    if manual > 0:
        blocks = manual
    if not blocks:
        return "", None
    if logger:
        logger.info(
            "DIAMOND memory plan: %.0f GB total, reserve %.0f GB, budget %.0f GB "
            "-> -c %d -b %d (~%.0f GB per block)",
            total, total - budget, budget, chunks, blocks, gb_per_block)
    return (f"--index-chunks {chunks} --block-size {blocks}", blocks)



def write_empty_report(args, _profile, thr):
    """
    Write the empty (`no pathogen detected`) report files.
    """
    output_tsv = os.path.join(args.output, "detection_results.tsv")
    output_summary = os.path.join(args.output, "detection_summary.txt")
    with open(output_tsv, 'w', encoding='utf-8') as f:
        f.write("Virus_Seq\tBest_Hit_Seq\tCoverage(%)\tDepth\tContig_ID\tContig_Length\t"
                "Reads_Mapped\tMean_Depth\tMax_Depth\tBreadth_Pct\tPident\tQcov\tEvalue\t"
                "Genome_Breadth_Pct\tN_Contigs_Blastn_Pass\tCall\tMargin_Fold\tDecision_Rule\tFlags\t"
                "Virus_Species\tVirus_Acronym\tVirus_Family\tVirus_Genus\n")
    with open(output_summary, 'w', encoding='utf-8') as f:
        f.write("No pathogen detected.\n")
    with open(os.path.join(args.output, "detection_virus_summary.tsv"), 'w',
              encoding='utf-8') as f:
        f.write("%s" % VIRUS_SUMMARY_HEADER)
    with open(os.path.join(args.output, "detection_thresholds.json"), 'w',
              encoding='utf-8') as f:
        json.dump({'profile': _profile, 'thresholds': thr,
                   'call_levels': ['POSITIVE', 'INDETERMINATE', 'NEGATIVE']}, f, indent=2)
    return None


def write_fasta_subset(src_fa, ids, out_fa):
    """
    Write the records of src_fa whose id is in `ids` to out_fa; return the count.
    """
    n = 0
    write = False
    with open(src_fa, 'r', encoding='utf-8') as fh, \
            open(out_fa, 'w', encoding='utf-8') as out:
        for line in fh:
            if line.startswith('>'):
                cid = line[1:].strip().split()[0]
                write = cid in ids
                if write:
                    n += 1
            if write:
                out.write(line)
    return n


def run_nr_diamond(query_fa, db_path, out_tsv, threads, evalue,
                   max_target_seqs, sensitivity, mem_args, diamond_bin, logger):
    """
    Stage 2 nr blastx: runs only on the contigs BLASTN has confirmed.
    """
    cmd = f"""{os.path.expanduser(diamond_bin)} blastx \\
        --threads {threads} {mem_args} \\
        --query {query_fa} \\
        --db {os.path.expanduser(db_path)} \\
        --out {out_tsv} \\
        {sensitivity} \\
        --evalue {evalue} \\
        --max-target-seqs {max_target_seqs} \\
        --outfmt {DIAMOND_OUTFMT}
    """
    logger.info("Running nr DIAMOND: %s", out_tsv)
    result = subprocess.run(cmd, shell=True, capture_output=True,
                            encoding='utf-8', errors='replace')
    if result.returncode != 0:
        logger.error(f"nr DIAMOND failed: {result.stderr[:2000]}")
        return False
    logger.info(f"nr DIAMOND completed -> {out_tsv}")
    return True


def run_blastn_first_flow(args, taxids, logger, checkpoints, checkpoint_file,
                          _profile, thr, config_data):
    """
    Current order: BLASTN on every stage-1 candidate, then nr DIAMOND on the
    survivors.

    The final set is the same as with the old order (both gates are ANDed), but
    negative libraries often have zero BLASTN hits, so the whole nr step can be
    skipped. Returns filtered_contigs, or None when an empty result was written.
    """
    blastn_raw = os.path.join(args.output, "blastn_raw.tsv")
    blastn_filtered = os.path.join(args.output, "blastn_results.tsv")
    virus_contigs_fa = os.path.join(args.output, "virus_contigs.fa")
    diamond_tsv = os.path.join(args.output, "diamond_nr.tsv")

    logger.info("Step 2: BLASTN on %d prescreen candidate(s)...",
                _count_fasta(args.candidates))
    if not should_skip_step("step2_blastn", checkpoints, [blastn_filtered]):
        if not run_blastn(args.candidates, args.virus_index, blastn_raw,
                          args.threads, args.blastn_evalue, args.blastn_bin, logger):
            logger.error("BLASTN failed!")
            sys.exit(1)
        kept = filter_blastn_results(blastn_raw, blastn_filtered,
                                     args.blastn_identity, args.blastn_coverage, logger)
        checkpoints["completed_steps"].append("step2_blastn")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 2 (completed)")
        kept = _count_rows(blastn_filtered)
    logger.info("BLASTN pass contigs: %d", kept)
    if kept == 0:
        diagnose_blastn_filter(blastn_raw, args.blastn_identity,
                               args.blastn_coverage, logger)
        logger.info("DIAMOND prefilter: 0 taxid-matched hit(s) on 0 contig(s); "
                    "nr DIAMOND skipped (no BLASTN hit)")
        logger.info("Found 0 virus-related contigs after filtering")
        logger.info("nr pass contigs: 0")
        logger.warning("No contig passed BLASTN; nr DIAMOND will be skipped")
        return write_empty_report(args, _profile, thr)

    if not should_skip_step("step3_confirmed", checkpoints, [virus_contigs_fa]):
        confirmed = set()
        with open(blastn_filtered, 'r', encoding='utf-8') as fh:
            for line in fh:
                p = line.split('\t', 1)
                if p and p[0].strip():
                    confirmed.add(p[0].strip())
        n = write_fasta_subset(args.candidates, confirmed, virus_contigs_fa)
        logger.info(f"Step 3: {n} contig(s) confirmed by BLASTN -> {virus_contigs_fa}")
        checkpoints["completed_steps"].append("step3_confirmed")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 3 (completed)")

    step4_cache = os.path.join(args.output, "step4_nr_diamond.json")
    if not should_skip_step("step4_nr_diamond", checkpoints, [diamond_tsv, step4_cache]):
        if not args.diamond_db:
            logger.error("--diamond-db is required for the nr DIAMOND step")
            sys.exit(1)
        mem_args, blocks = plan_diamond_mem(config_data, args.node_mem_gb, logger)
        if not run_nr_diamond(virus_contigs_fa, args.diamond_db, diamond_tsv,
                              args.threads, args.diamond_evalue,
                              args.diamond_max_target_seqs,
                              args.diamond_sensitivity, mem_args,
                              args.diamond_bin, logger):
            sys.exit(1)
        logger.info("Step 4: filtering nr DIAMOND results by taxid/pident/evalue")
        filtered = filter_diamond(diamond_tsv, taxids, args.pident, args.evalue, logger)
        diagnose_diamond_filter(diamond_tsv, taxids, args.pident, args.evalue, logger)
        logger.info("nr pass contigs: %d", len(filtered))
        with open(step4_cache, 'w', encoding='utf-8') as f:
            json.dump({"filtered": filtered, "blocks": blocks}, f, indent=2)
        checkpoints["completed_steps"].append("step4_nr_diamond")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 4 (completed)")
        with open(step4_cache, 'r', encoding='utf-8') as f:
            filtered = json.load(f).get("filtered", {})

    logger.info(f"Found {len(filtered)} virus-related contigs after filtering")
    if not filtered:
        logger.warning("nr DIAMOND found no target-taxid hit "
                       "(BLASTN passed %d contig(s))", _count_fasta(virus_contigs_fa))
        return write_empty_report(args, _profile, thr)
    return filtered


def _count_fasta(path):
    n = 0
    try:
        with open(path, 'r', encoding='utf-8') as fh:
            for line in fh:
                if line.startswith('>'):
                    n += 1
    except OSError:
        return 0
    return n


def _count_rows(path):
    """
    Count non-empty rows of a headerless TSV (resume restores the BLASTN gate
    count).
    """
    n = 0
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.strip():
                    n += 1
    except OSError:
        return 0
    return n


def run_legacy_flow(args, taxids, logger, checkpoints, checkpoint_file,
                    _profile, thr):
    """
    Legacy order (only --diamond TSV given): nr filter -> extract contigs ->
    BLASTN.

    Returns filtered_contigs, or None when an empty result was written.
    """
    # Step 2: Filter diamond results
    step2_filtered_cache = os.path.join(args.output, "step2_filtered_contigs.json")
    
    if not should_skip_step("step2_filter", checkpoints, [step2_filtered_cache]):
        logger.info("Step 2: Filtering diamond results...")
        logger.info(f"Filters: pident>={args.pident}%, evalue<={args.evalue}")
        filtered_contigs = filter_diamond(args.diamond, taxids, args.pident,
                                          args.evalue, logger)
        logger.info(f"Found {len(filtered_contigs)} virus-related contigs after filtering")
        logger.info("nr pass contigs: %d", len(filtered_contigs))
        diagnose_diamond_filter(args.diamond, taxids, args.pident,
                                args.evalue, logger)
        
        # Save filtered_contigs for checkpoint resume
        with open(step2_filtered_cache, 'w', encoding='utf-8') as f:
            json.dump(filtered_contigs, f, indent=2)
        
        checkpoints["completed_steps"].append("step2_filter")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 2 (completed)")
        with open(step2_filtered_cache, 'r', encoding='utf-8') as f:
            filtered_contigs = json.load(f)
        logger.info(f"Loaded {len(filtered_contigs)} filtered contigs from cache")
    
    if filtered_contigs:
        logger.info("-" * 50)
        logger.info("Filtered contigs summary:")
        logger.info(f"  Total contigs: {len(filtered_contigs)}")
        
        pidents = [info['pident'] for info in filtered_contigs.values()]
        qcovs = [info['qcov'] for info in filtered_contigs.values()]
        
        logger.info(f"  Identity range: {min(pidents):.1f}% - {max(pidents):.1f}% (avg: {sum(pidents)/len(pidents):.1f}%)")
        logger.info(f"  Coverage range: {min(qcovs):.1f}% - {max(qcovs):.1f}% (avg: {sum(qcovs)/len(qcovs):.1f}%)")
        
        virus_groups = defaultdict(list)
        for contig_id, info in filtered_contigs.items():
            virus_groups[info['sscinames']].append(contig_id)
        
        logger.info(f"  Virus groups: {len(virus_groups)}")
        for virus_name, contig_list in sorted(virus_groups.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
            logger.info(f"    - {virus_name}: {len(contig_list)} contigs")
        
        logger.info("-" * 50)
    
    if not filtered_contigs:
        logger.warning("No virus-related contigs found!")
        return write_empty_report(args, _profile, thr)
    # Step 3: Extract virus-related contigs FASTA
    virus_contigs_fa = os.path.join(args.output, "virus_contigs.fa")
    
    if not should_skip_step("step3_extract", checkpoints, [virus_contigs_fa]):
        logger.info("Step 3: Extracting virus-related contigs FASTA...")
        
        # Get megahit contigs file
        megahit_dir_name = os.path.basename(args.megahit.rstrip('/'))
        if '_megahit_result' in megahit_dir_name:
            contig_name = megahit_dir_name.replace('_megahit_result', '')
        else:
            contig_name = "results"
        
        possible_contig_files = [
            os.path.join(args.megahit, f"{contig_name}.contigs.fa"),
            os.path.join(args.megahit, "results.contigs.fa"),
            os.path.join(args.megahit, "final.contigs.fa")
        ]
        megahit_contigs_fa = None
        for f in possible_contig_files:
            if os.path.exists(f):
                megahit_contigs_fa = f
                break
        
        if not megahit_contigs_fa:
            logger.error("No megahit contigs file found!")
            sys.exit(1)
        
        # Extract virus-related contigs to separate FASTA
        contig_ids_set = set(filtered_contigs.keys())
        extracted_count = 0
        with open(virus_contigs_fa, 'w', encoding='utf-8') as out:
            write_seq = False
            with open(megahit_contigs_fa, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('>'):
                        contig_id = line[1:].strip().split()[0]
                        if contig_id in contig_ids_set:
                            out.write(line)
                            write_seq = True
                            extracted_count += 1
                        else:
                            write_seq = False
                    elif write_seq:
                        out.write(line)
        logger.info(f"Extracted {extracted_count} virus-related contigs to {virus_contigs_fa}")
        
        checkpoints["completed_steps"].append("step3_extract")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 3 (completed)")
    
    # Step 4: BLASTN alignment (contigs vs virus reference)
    logger.info("Step 4: Running BLASTN alignment...")
    blastn_raw = os.path.join(args.output, "blastn_raw.tsv")
    blastn_filtered = os.path.join(args.output, "blastn_results.tsv")
    
    if not should_skip_step("step4_blastn", checkpoints, [blastn_filtered]):
        if run_blastn(virus_contigs_fa, args.virus_index, blastn_raw, 
                     args.threads, args.blastn_evalue, args.blastn_bin, logger):
            kept = filter_blastn_results(blastn_raw, blastn_filtered, 
                                 args.blastn_identity, args.blastn_coverage, logger)
            logger.info("BLASTN pass contigs: %d", kept)
            if kept == 0 and filtered_contigs:
                logger.warning(
                    f"Step 4: {len(filtered_contigs)} contig(s) entered BLASTN, "
                    f"but 0 hits passed the filters - the report will be empty")
                diagnose_blastn_filter(blastn_raw, args.blastn_identity,
                                       args.blastn_coverage, logger)
            checkpoints["completed_steps"].append("step4_blastn")
            save_checkpoints(checkpoints, checkpoint_file)
        else:
            logger.error("BLASTN failed!")
            sys.exit(1)
    else:
        logger.info("Skipping Step 4 (completed)")
    
    return filtered_contigs


def main():
    parser = argparse.ArgumentParser(description="Pathogen Detection Pipeline")
    parser.add_argument("--diamond", default=None,
                        help="Legacy: pre-computed DIAMOND-vs-nr TSV "
                             "(only used when --candidates is not given)")
    parser.add_argument("--candidates", default=None,
                        help="pre-screen candidate FASTA from the pre-processing stage "
                             "(preprocess/prescreen_candidates.fasta). When set, "
                             "BLASTN runs FIRST and nr DIAMOND only on the contigs "
                             "that passed BLASTN")
    parser.add_argument("--diamond-db", default=None,
                        help="DIAMOND nr database for the post-BLASTN stage "
                             "(default: config diamond.db_path)")
    parser.add_argument("--diamond-bin", default=None,
                        help="DIAMOND executable (default: config diamond.bin)")
    parser.add_argument("--diamond-evalue", default="1e-5",
                        help="nr DIAMOND search e-value (default 1e-5)")
    parser.add_argument("--diamond-max-target-seqs", type=int, default=10,
                        help="nr DIAMOND --max-target-seqs (default 10)")
    parser.add_argument("--diamond-sensitivity", default="--sensitive",
                        help="nr DIAMOND sensitivity flag (default --sensitive)")
    parser.add_argument("--node-mem-gb", type=float, default=None,
                        help="Node memory in GB, used to derive nr DIAMOND -b "
                             "(moved here from Tab 2)")
    parser.add_argument("--fq1", required=True, help="Original forward reads (FASTQ)")
    parser.add_argument("--fq2", required=True, help="Original reverse reads (FASTQ)")
    parser.add_argument("--megahit", required=True, help="Megahit results directory")
    parser.add_argument("--taxids", default=None,
                        help="Target taxid list file (default: config "
                             "databases.virus_taxids)")
    parser.add_argument("--virus-index", default=None,
                        help="BLASTN virus database prefix (default: config "
                             "databases.blastn_index)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("-c", "--config", default=None,
                        help="config2.yaml; explicit --*-bin options take precedence, "
                             "otherwise each tool's `bin` is read from this file")
    
    parser.add_argument("--pident", type=float, default=50,
                        help="Diamond identity threshold (%%); default matches the "
                             "archived run (50)")
    parser.add_argument("--evalue", type=float, default=1e-2,
                        help="Diamond e-value threshold; default matches the "
                             "archived run (1e-02)")
    
    parser.add_argument("--blastn-identity", type=float, default=70, help="BLASTN identity threshold (%%)")
    parser.add_argument("--blastn-evalue", type=float, default=1e-4, help="BLASTN e-value threshold")
    parser.add_argument("--blastn-coverage", type=float, default=70, help="BLASTN coverage threshold (%%)")

    parser.add_argument("--threshold-profile", default="P0_current_default",
                        choices=sorted(THRESHOLD_PROFILES.keys()),
                        help="Threshold preset: P0_current_default (no gating, legacy) "
                             "/ P1_two_tier / P2_stringent")
    parser.add_argument("--min-contig-len", type=int, default=None,
                        help="Minimum contig length (bp) to be considered")
    parser.add_argument("--min-reads-mapped", type=int, default=None,
                        help="Minimum reads mapped to the virus contigs (per virus)")
    parser.add_argument("--min-mean-depth", type=float, default=None,
                        help="Minimum mean depth over the virus contigs")
    parser.add_argument("--min-contig-breadth", type=float, default=None,
                        help="Minimum breadth (%%) of the best virus contig")
    parser.add_argument("--min-genome-breadth", type=float, default=None,
                        help="Minimum breadth (%%) over the reference genome (from BLASTN slen)")
    parser.add_argument("--margin-positive", type=float, default=None,
                        help="A metric must exceed its threshold by this fold to call POSITIVE")
    parser.add_argument("--min-supporting-contigs", type=int, default=None,
                        help="Minimum number of supporting contigs per virus")
    parser.add_argument("--blastdbcmd-bin", default=None,
                        help="blastdbcmd executable (only used if BLASTN output has no slen)")

    parser.add_argument("--species-table", default=None,
                        help="Optional species-composition CSV "
                             "(Family,Genus,Species,Acronym,N sequences); "
                             "lines starting with '#' are comments")
    parser.add_argument("--species-alias", default=None,
                        help="CSV of alias,acronym used to map taxid names to the species table")
    parser.add_argument("--verify-db", action="store_true",
                        help="Count BLAST DB entries with blastdbcmd and compare to the species-table total")
    
    parser.add_argument("--threads", type=int, default=64, help="Number of threads")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose/debug logging")
    parser.add_argument("--no-resume", action="store_true", help="Ignore checkpoint and start fresh")
    
    parser.add_argument("--bowtie2-bin", default=None, help="bowtie2 executable path")
    parser.add_argument("--samtools-bin", default=None, help="samtools executable path")
    parser.add_argument("--bedtools-bin", default=None, help="bedtools executable path")
    parser.add_argument("--bedgraph-to-bigwig-bin", default=None,
                        help="bedGraphToBigWig executable path")
    parser.add_argument("--blastn-bin", default=None, help="blastn executable path")
    
    args = parser.parse_args()

    config_data = load_tool_config(args.config)
    args.taxids = resolve_db_path(args.taxids, config_data, "virus_taxids")
    args.virus_index = resolve_db_path(args.virus_index, config_data, "blastn_index")
    _missing = [flag for flag, val in (("--taxids", args.taxids),
                                       ("--virus-index", args.virus_index))
                if not val]
    if _missing:
        parser.error("missing database path(s): " + ", ".join(_missing)
                     + " — pass them on the command line or set "
                       "databases.virus_taxids / databases.blastn_index in the config")

    _dm = ((config_data or {}).get('diamond') or {})
    args.diamond_db = args.diamond_db or _dm.get('db_path')
    args.diamond_bin = args.diamond_bin or _dm.get('bin') or 'diamond'
    if not args.candidates and not args.diamond:
        parser.error("either --candidates (current order) or --diamond "
                     "(legacy pre-computed TSV) is required")

    _resolved = resolve_tool_bins(
        config_data, {attr: getattr(args, attr, None) for attr in BIN_SPECS})
    for _attr, _path in _resolved.items():
        setattr(args, _attr, _path)

    _profile = args.threshold_profile
    thr = dict(THRESHOLD_PROFILES[_profile])
    for _key in list(thr.keys()):
        _val = getattr(args, _key, None)
        if _val is not None:
            thr[_key] = _val

    _negative = sorted(k for k, v in thr.items()
                       if isinstance(v, (int, float)) and v < 0)
    if _negative:
        parser.error("negative threshold value(s): " + ", ".join(_negative))
    if thr["margin_positive"] <= 0:
        parser.error("--margin-positive must be > 0")

    # Create output directory FIRST (before logger tries to create log file)
    os.makedirs(args.output, exist_ok=True)
    
    # Setup logger (log file can now be created in existing directory)
    log_file = os.path.join(args.output, f"pathogen_detect_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    logger = setup_logger("PathogenDetect", log_file)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.info("Verbose logging enabled")

    logger.info(f"Threshold profile: {_profile}")
    for _k, _v in thr.items():
        logger.info(f"  {_k:24s} = {_v}")

    if args.config:
        if config_data:
            logger.info(f"Tool paths from config: {args.config}")
        else:
            logger.warning(f"Could not read config (or PyYAML missing): {args.config}; "
                           f"falling back to built-in tool names")
    else:
        logger.warning("No --config given; tool paths fall back to built-in names "
                       "(pass -c server/configs/config2.yaml to read them from config)")
    for _attr in ('bowtie2_bin', 'samtools_bin', 'bedtools_bin', 'blastn_bin',
                  'bedgraph_to_bigwig_bin', 'blastdbcmd_bin'):
        logger.info(f"  {_attr:24s} = {getattr(args, _attr)}")

    species_rows, species_index = [], {}
    default_alias = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "configs", "species_aliases.csv")
    alias_path = args.species_alias or default_alias
    if args.species_table:
        if not os.path.exists(args.species_table):
            logger.warning(f"Species table not found: {args.species_table}")
        else:
            try:
                species_rows = load_species_table(args.species_table)
                species_index = build_species_index(species_rows, alias_path)
                total_seq = sum(r['n_sequences'] for r in species_rows)
                logger.info(f"Species table: {len(species_rows)} species, "
                            f"{total_seq} sequences (from --species-table)")
                for _r in species_rows:
                    logger.info(f"  {_r['acronym']:6s} {_r['species']:45s} "
                                f"{_r['n_sequences']}")
                if args.verify_db:
                    n_db = count_db_sequences(args.virus_index, args.blastdbcmd_bin)
                    if n_db is None:
                        logger.warning("Could not count BLAST DB entries "
                                       "(blastdbcmd unavailable)")
                    elif n_db != total_seq:
                        logger.warning(f"DB entry count {n_db} != species-table total {total_seq}")
                    else:
                        logger.info(f"DB entry count matches species-table total: {n_db}")
            except Exception as exc:
                logger.warning(f"Failed to load species table: {exc}")
    
    checkpoint_file = os.path.join(args.output, "pathogen_detect.checkpoints.json")
    
    # Auto-detect checkpoint file
    if args.no_resume and os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)
        checkpoints = {"completed_steps": []}
        logger.info("Checkpoint file removed. Starting fresh.")
    elif os.path.exists(checkpoint_file):
        checkpoints = load_checkpoints(checkpoint_file)
        logger.info(f"Found checkpoint file: {checkpoint_file}")
        logger.info(f"Completed steps: {checkpoints.get('completed_steps', [])}")
        logger.info("Use --no-resume to ignore checkpoint and start fresh")
    else:
        checkpoints = {"completed_steps": []}
    
    logger.info(f"Output directory ready: {args.output}")
    logger.info("=" * 60)
    logger.info("Pathogen Detection Pipeline Started")
    logger.info(f"Diamond results: {args.diamond}")
    logger.info(f"Original reads: {args.fq1}, {args.fq2}")
    logger.info(f"Megahit results: {args.megahit}")
    logger.info(f"Taxid list: {args.taxids}")
    logger.info(f"Virus index (BLASTN db): {args.virus_index}")
    logger.info(f"Output: {args.output}")
    logger.info("=" * 60)
    
    start_time = time.time()
    
    pre_check(args, logger)
    
    # Step 1: Load taxids
    step1_taxids_cache = os.path.join(args.output, "step1_taxids_cache.json")
    
    if not should_skip_step("step1_load_taxids", checkpoints, [step1_taxids_cache]):
        logger.info("Step 1: Loading taxid list...")
        taxids, taxid_names = load_taxids(args.taxids)
        logger.info(f"Loaded {len(taxids)} target taxids")
        
        # Save taxids result for checkpoint resume
        with open(step1_taxids_cache, 'w', encoding='utf-8') as f:
            json.dump({"taxids": list(taxids), "taxid_names": taxid_names}, f, indent=2)
        
        checkpoints["completed_steps"].append("step1_load_taxids")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 1 (completed)")
        with open(step1_taxids_cache, 'r', encoding='utf-8') as f:
            cache = json.load(f)
            taxids = set(cache["taxids"])
            taxid_names = cache["taxid_names"]
        logger.info(f"Loaded {len(taxids)} target taxids from cache")
    
    # Log taxid names if available (show first 5 with ellipsis if more)
    if taxid_names:
        logger.info("Target viruses:")
        display_items = list(taxid_names.items())[:5]
        for taxid, name in display_items:
            logger.info(f"  - {taxid}: {name}")
        if len(taxid_names) > 5:
            logger.info(f"  ... and {len(taxid_names) - 5} more")
    elif taxids:
        logger.info("Target taxids:")
        display_taxids = list(taxids)[:5]
        for taxid in display_taxids:
            logger.info(f"  - {taxid}")
        if len(taxids) > 5:
            logger.info(f"  ... and {len(taxids) - 5} more")
    
    taxid_map = taxid_names
    
    blastn_raw = os.path.join(args.output, "blastn_raw.tsv")
    blastn_filtered = os.path.join(args.output, "blastn_results.tsv")
    virus_contigs_fa = os.path.join(args.output, "virus_contigs.fa")
    if args.candidates:
        filtered_contigs = run_blastn_first_flow(
            args, taxids, logger, checkpoints, checkpoint_file,
            _profile, thr, config_data)
    else:
        filtered_contigs = run_legacy_flow(
            args, taxids, logger, checkpoints, checkpoint_file,
            _profile, thr)
    if filtered_contigs is None:
        return

    # Step 5: Calculate contig coverage
    # Use clean reads from Tab 2 if available, otherwise use original reads
    preprocess_dir = os.path.dirname(args.megahit.rstrip('/'))
    clean_reads_dir = os.path.join(preprocess_dir, "clean_reads")
    clean_fq1 = os.path.join(clean_reads_dir, "clean_R1.fq.gz")
    clean_fq2 = os.path.join(clean_reads_dir, "clean_R2.fq.gz")
    
    if os.path.exists(clean_fq1) and os.path.exists(clean_fq2):
        fq1 = clean_fq1
        fq2 = clean_fq2
        logger.info(f"Using clean reads from: {clean_reads_dir}")
    else:
        fq1 = args.fq1
        fq2 = args.fq2
        logger.info(f"Using original reads: {args.fq1}, {args.fq2}")
    
    logger.info("Step 5: Calculating contig coverage...")
    contigs_bam = os.path.join(args.output, "contigs_mapping.bam")
    coverage_summary = os.path.join(args.output, "contig_coverage_summary.tsv")
    coverage_bw = os.path.join(args.output, "contig_coverage.bw")
    
    if not should_skip_step("step5_coverage", checkpoints, [contigs_bam, coverage_summary]):
        if not calculate_contig_coverage_with_bigwig(
            fq1, fq2, virus_contigs_fa,
            contigs_bam, coverage_summary, coverage_bw,
            args.threads, args.bowtie2_bin, args.samtools_bin, 
            args.bedtools_bin, args.bedgraph_to_bigwig_bin, logger
        ):
            logger.error("Coverage calculation failed!")
            sys.exit(1)
        checkpoints["completed_steps"].append("step5_coverage")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 5 (completed)")
    
    # Step 6: Generate report
    output_tsv = os.path.join(args.output, "detection_results.tsv")
    output_summary = os.path.join(args.output, "detection_summary.txt")
    output_virus_tsv = os.path.join(args.output, "detection_virus_summary.tsv")
    output_thresholds = os.path.join(args.output, "detection_thresholds.json")

    results = []

    if not should_skip_step("step6_report", checkpoints,
                            [output_tsv, output_summary, output_virus_tsv]):
        logger.info("Step 6: Generating detection report...")

        best_hits = parse_blastn_best_hits(blastn_filtered)
        coverage_stats = parse_coverage_summary(coverage_summary)

        if args.candidates:
            _keep = set(filtered_contigs)
            best_hits = {c: h for c, h in best_hits.items() if c in _keep}
            logger.info("Step 6: %d contig(s) passed both BLASTN and nr DIAMOND",
                        len(best_hits))

        if not best_hits and filtered_contigs:
            logger.warning(
                f"Step 6: blastn_results.tsv is empty, so the report has header "
                f"only (Step 2 produced {len(filtered_contigs)} contig(s)). "
                f"Report rows come from BLASTN results only; please check "
                f"--blastn-identity / --blastn-coverage")

        per_contig = []
        for contig_id, hit in best_hits.items():
            cov = coverage_stats.get(contig_id, {})
            per_contig.append({
                'Contig_ID': contig_id,
                'Best_Hit_Seq': hit['sseqid'],
                'Contig_Length': cov.get('length', 0),
                'Reads_Mapped': cov.get('reads', 0),
                'Mean_Depth': cov.get('mean_depth', 0.0),
                'Max_Depth': cov.get('max_depth', 0),
                'Breadth_Pct': cov.get('breadth', 0.0),
                'Pident': hit['pident'],
                'Qcov': hit['qcov'],
                'Evalue': hit['evalue'],
                'interval': hit.get('interval'),
            })

        subject_lengths = {}
        for cid, hit in best_hits.items():
            if hit.get('slen'):
                subject_lengths[hit['sseqid']] = hit['slen']
        missing = {c['Best_Hit_Seq'] for c in per_contig} - set(subject_lengths)
        if missing:
            logger.warning(
                f"BLASTN output lacks slen for {len(missing)} subject(s); "
                f"falling back to blastdbcmd")
            subject_lengths.update(get_subject_lengths(
                args.virus_index, missing, args.blastdbcmd_bin, logger))
        decisions = evaluate_detection(per_contig, subject_lengths, thr, logger)

        contig_taxon = {}
        if species_rows:
            for cid, info in filtered_contigs.items():
                candidates = [info.get('sscinames', '')]
                for t in str(info.get('staxids', '')).split(';'):
                    t = t.strip()
                    if t in taxid_map:
                        candidates.append(taxid_map[t])
                for nm in candidates:
                    row = match_species(nm, species_index, species_rows)
                    if row:
                        contig_taxon[cid] = row
                        break
            missing = [c['Contig_ID'] for c in per_contig
                       if c['Contig_ID'] not in contig_taxon]
            if missing:
                logger.warning(f"Species annotation missing for {len(missing)} contig(s); "
                               f"extend aliases via --species-alias")

        subject_taxon = {}
        for c in per_contig:
            if c['Best_Hit_Seq'] not in subject_taxon and c['Contig_ID'] in contig_taxon:
                subject_taxon[c['Best_Hit_Seq']] = contig_taxon[c['Contig_ID']]

        header = ['Virus_Seq', 'Best_Hit_Seq', 'Coverage(%)', 'Depth', 'Contig_ID',
                  'Contig_Length', 'Reads_Mapped', 'Mean_Depth', 'Max_Depth',
                  'Breadth_Pct', 'Pident', 'Qcov', 'Evalue',
                  'Genome_Breadth_Pct', 'N_Contigs_Blastn_Pass', 'Call', 'Margin_Fold',
                  'Decision_Rule', 'Flags',
                  'Virus_Species', 'Virus_Acronym', 'Virus_Family', 'Virus_Genus']
        results = []
        with open(output_tsv, 'w', encoding='utf-8') as f:
            f.write('\t'.join(header) + '\n')
            for c in sorted(per_contig, key=lambda x: (x['Best_Hit_Seq'], x['Contig_ID'])):
                d = decisions.get(c['Best_Hit_Seq'], {})
                gb = d.get('genome_breadth_pct')
                mg = d.get('margin')
                tax = contig_taxon.get(c['Contig_ID'])
                row = {
                    'Virus_Seq': c['Best_Hit_Seq'],
                    'Best_Hit_Seq': c['Best_Hit_Seq'],
                    'Coverage(%)': f"{c['Breadth_Pct']:.2f}",
                    'Depth': f"{c['Mean_Depth']:.2f}",
                    'Contig_ID': c['Contig_ID'],
                    'Contig_Length': c['Contig_Length'],
                    'Reads_Mapped': c['Reads_Mapped'],
                    'Mean_Depth': f"{c['Mean_Depth']:.2f}",
                    'Max_Depth': c['Max_Depth'],
                    'Breadth_Pct': f"{c['Breadth_Pct']:.2f}",
                    'Pident': f"{c['Pident']:.2f}",
                    'Qcov': f"{c['Qcov']:.2f}",
                    'Evalue': f"{c['Evalue']:.3g}",
                    'Genome_Breadth_Pct': f"{gb:.2f}" if gb is not None else 'NA',
                    'N_Contigs_Blastn_Pass': d.get('n_contigs', 0),
                    'Call': d.get('call', 'NEGATIVE'),
                    'Margin_Fold': f"{mg:.2f}" if mg is not None else 'NA',
                    'Decision_Rule': _profile,
                    'Flags': d.get('flags', 'NA'),
                    'Virus_Species': tax['species'] if tax else 'NA',
                    'Virus_Acronym': tax['acronym'] if tax else 'NA',
                    'Virus_Family': tax['family'] if tax else 'NA',
                    'Virus_Genus': tax['genus'] if tax else 'NA',
                }
                results.append(row)
                f.write('\t'.join(str(row[k]) for k in header) + '\n')

        with open(output_virus_tsv, 'w', encoding='utf-8') as f:
            f.write(VIRUS_SUMMARY_HEADER)
            for subject in sorted(decisions):
                d = decisions[subject]
                tax = subject_taxon.get(subject)
                f.write('\t'.join(str(x) for x in [
                    subject, d['call'], d['n_contigs'], d['n_genome_regions'],
                    d['best_contig_id'] or 'NA',
                    d['best_contig_accession'] or 'NA',
                    d['best_contig_len'] if d['best_contig_len'] is not None else 'NA',
                    f"{d['best_contig_mean_depth']:.2f}"
                    if d['best_contig_mean_depth'] is not None else 'NA',
                    f"{d['best_contig_breadth_pct']:.2f}"
                    if d['best_contig_breadth_pct'] is not None else 'NA',
                    d['best_contig_max_depth'] if d['best_contig_max_depth'] is not None else 'NA',
                    d['reads_mapped'], d['total_contig_len'],
                    d['span_union_bp'],
                    f"{d['span_overlap_ratio']:.4f}"
                    if d['span_overlap_ratio'] is not None else 'NA',
                    f"{d['genome_breadth_pct']:.2f}" if d['genome_breadth_pct'] is not None else 'NA',
                    d['subject_length'] if d['subject_length'] else 'NA',
                    f"{d['margin']:.2f}" if d['margin'] is not None else 'NA',
                    _profile, d['flags'], d['reason'],
                    tax['species'] if tax else 'NA',
                    tax['acronym'] if tax else 'NA',
                    tax['family'] if tax else 'NA',
                    tax['genus'] if tax else 'NA',
                ]) + '\n')

        with open(output_thresholds, 'w', encoding='utf-8') as f:
            json.dump({
                'profile': _profile,
                'thresholds': thr,
                'call_levels': ['POSITIVE', 'INDETERMINATE', 'NEGATIVE'],
                'note': 'Coverage(%) reports breadth; Depth reports mean depth',
            }, f, indent=2)

        by_call = defaultdict(int)
        for d in decisions.values():
            by_call[d['call']] += 1
        with open(output_summary, 'w', encoding='utf-8') as f:
            f.write("Pathogen Detection Summary\n")
            f.write("=" * 60 + "\n")
            f.write(f"Threshold profile : {_profile}\n")
            f.write(f"Contigs reported  : {len(results)}\n")
            f.write(f"Viruses by call   : POSITIVE={by_call['POSITIVE']}, "
                    f"INDETERMINATE={by_call['INDETERMINATE']}, "
                    f"NEGATIVE={by_call['NEGATIVE']}\n")
            if not results and filtered_contigs:
                f.write(
                    f"NOTE: Step 2 produced {len(filtered_contigs)} contig(s), "
                    f"but BLASTN found no hit "
                    f"(identity>={args.blastn_identity}%, "
                    f"coverage>={args.blastn_coverage}%).\n"
                        f"      The report is therefore empty; see 'Filtered BLASTN "
                        f"results' and the diagnostic lines in the run log.\n")
            f.write("\n")
            for subject in sorted(decisions,
                                  key=lambda s: (decisions[s]['call'] != 'POSITIVE', s)):
                d = decisions[subject]
                gb = d['genome_breadth_pct']
                gb_txt = 'NA' if gb is None else f"{gb:.1f}%"
                mg = d['margin']
                mg_txt = 'NA' if mg is None else f"{mg:.2f}"
                f.write(f"- {subject}: {d['call']}\n")
                f.write(f"  contigs={d['n_contigs']} regions={d['n_genome_regions']} "
                        f"reads={d['reads_mapped']} depth={d['mean_depth']:.2f} "
                        f"contig_breadth={d['contig_breadth_max']:.1f}% "
                        f"genome_breadth={gb_txt} margin={mg_txt}\n")
                f.write(f"  flags={d['flags']}")
                if d['reason']:
                    f.write(f"  reason={d['reason']}")
                f.write("\n")

        checkpoints["completed_steps"].append("step6_report")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 6 (completed)")
    
    # Pipeline completed
    total_time = time.time() - start_time
    logger.info("=" * 60)
    logger.info(f"Pipeline completed in {total_time:.2f} seconds")
    logger.info(f"Detected {len(results)} virus contig(s)")
    logger.info(f"Results saved to: {args.output}")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
