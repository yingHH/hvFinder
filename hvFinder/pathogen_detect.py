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
import time
import argparse
import subprocess
import shutil
import json
from collections import defaultdict
from datetime import datetime

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
    """加载检查点"""
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as f:
            return json.load(f)
    return {"completed_steps": []}

def save_checkpoints(checkpoints, checkpoint_file):
    """保存检查点"""
    with open(checkpoint_file, 'w') as f:
        json.dump(checkpoints, f, indent=2)

def should_skip_step(step_name, checkpoints, output_files):
    """判断是否跳过步骤"""
    if step_name not in checkpoints.get("completed_steps", []):
        return False
    for f in output_files:
        if not os.path.exists(f):
            return False
    return True


def pre_check(args, logger):
    """全方位的环境与输入文件预检"""
    logger.info("--- [Phase: Pre-Check] Validating Environment ---")
    validation_failed = False
    
    input_files = {
        'diamond': args.diamond,
        'fq1': args.fq1,
        'fq2': args.fq2,
        'taxids': args.taxids,
        'megahit': args.megahit,
    }
    for name, path in input_files.items():
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
    
    with open(taxid_file, 'r') as f:
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


def filter_diamond(diamond_tsv, taxids, pident, qcov, evalue, logger=None):
    """Filter diamond results by taxid and parameters
    
    Args:
        diamond_tsv: Path to diamond TSV file
        taxids: Set of target taxids
        pident: Minimum identity threshold
        qcov: Minimum coverage threshold
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
    
    with open(diamond_tsv, 'r') as f:
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
            
            # Apply filters
            if pident_val < pident:
                continue
            if qcovhsp < qcov:
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


def extract_contigs_by_ids(input_fa, contig_ids, output_fa):
    """Extract specific contigs by ID from a FASTA file"""
    contig_ids_set = set(contig_ids)
    extracted_count = 0
    
    with open(output_fa, 'w') as out:
        write_seq = False
        with open(input_fa, 'r') as f:
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
    
    print(f"INFO: Extracted {extracted_count} contigs to {output_fa}")
    return extracted_count > 0


def extract_contigs(megahit_dir, contig_name, output_fa):
    """Extract specific contigs from megahit results"""
    # Try different possible contig filenames
    possible_files = [
        os.path.join(megahit_dir, f"{contig_name}.contigs.fa"),
        os.path.join(megahit_dir, "results.contigs.fa"),
        os.path.join(megahit_dir, "final.contigs.fa")
    ]
    
    input_fa = None
    for f in possible_files:
        if os.path.exists(f):
            input_fa = f
            break
    
    if not input_fa:
        print(f"ERROR: No contigs file found in {megahit_dir}")
        return False
    
    # Extract specific contigs
    target_contigs = set()
    with open(input_fa, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig_id = line[1:].strip().split()[0]
                target_contigs.add(contig_id)
    
    if not target_contigs:
        print(f"ERROR: No contigs found in {input_fa}")
        return False
    
    with open(output_fa, 'w') as out:
        write_seq = False
        with open(input_fa, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig_id = line[1:].strip().split()[0]
                    if contig_id in target_contigs:
                        out.write(line)
                        write_seq = True
                    else:
                        write_seq = False
                elif write_seq:
                    out.write(line)
    
    print(f"INFO: Extracted {len(target_contigs)} contigs to {output_fa}")
    return True


def run_bowtie2(bowtie2_bin, index_path, query_fa, sam_file, threads, logger):
    """Run bowtie2 alignment"""
    cmd = f"""
{bowtie2_bin} -x {index_path} -f {query_fa} \
    -S {sam_file} \
    --threads {threads} \
    --local \
    --no-unal \
    -f
"""
    logger.info(f"Running bowtie2 alignment...")
    logger.debug(f"Command: {cmd}")
    
    try:
        proc = subprocess.run(
            cmd.strip(),
            shell=True,
            capture_output=True,
            text=True,
            executable='/bin/bash'
        )
        
        if proc.returncode == 0:
            logger.info(f"Bowtie2 alignment completed successfully")
            return True
        else:
            logger.error(f"Bowtie2 failed: {proc.stderr}")
            return False
    except Exception as e:
        logger.error(f"Exception during bowtie2: {str(e)}")
        return False


def run_blastn(query_fa, db_path, out_file, threads, evalue, blastn_bin, logger):
    """Run BLASTN alignment (contigs vs virus reference)"""
    blastn_path = os.path.expanduser(blastn_bin)
    cmd = f"""
    {blastn_path} -query {query_fa} \
           -db {db_path} \
           -out {out_file} \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
           -evalue {evalue} \
           -num_threads {threads}
    """
    logger.info("Running BLASTN alignment...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"BLASTN failed: {result.stderr}")
        return False
    logger.info(f"BLASTN completed. Results saved to {out_file}")
    return True


def filter_blastn_results(input_file, output_file, identity, coverage, logger):
    """Filter BLASTN results by identity and coverage"""
    count = 0
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue
            pident = float(fields[2])
            qcov = float(fields[12])
            if pident >= identity and qcov >= coverage:
                f_out.write(line)
                count += 1
    logger.info(f"Filtered BLASTN results: {count} hits saved to {output_file}")


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
    result = subprocess.run(cmd_index, shell=True, capture_output=True, text=True)
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
    result = subprocess.run(cmd_map, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Failed to map reads: {result.stderr}")
        return False
    
    logger.info("Converting SAM to BAM and sorting...")
    cmd_bam = f"{samtools_path} view -bS {sam_file} | {samtools_path} sort -o {output_bam}"
    result = subprocess.run(cmd_bam, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Failed to create BAM: {result.stderr}")
        return False
    
    cmd_index = f"{samtools_path} index {output_bam}"
    subprocess.run(cmd_index, shell=True, capture_output=True, text=True)
    
    logger.info("Calculating coverage summary...")
    cmd_summary = f"""
    {samtools_path} depth {output_bam} | awk '{{ 
        sum[$1]+=$3; count[$1]++; 
        if($3>max[$1]) max[$1]=$3; 
        if(min[$1]=="" || $3<min[$1]) min[$1]=$3
    }} END {{
        for(c in sum) print c"\\t"sum[c]/count[c]"\\t"min[c]"\\t"max[c]"\\t"count[c]
    }}' > {output_summary}
    """
    result = subprocess.run(cmd_summary, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Failed to calculate coverage: {result.stderr}")
        return False
    
    logger.info("Generating BigWig...")
    if shutil.which(bedgraph_to_bigwig_path) or os.path.exists(bedgraph_to_bigwig_path):
        cmd_bedgraph = f"{bedtools_path} genomecov -ibam {output_bam} -bg > {output_bam}.coverage.bedgraph"
        result = subprocess.run(cmd_bedgraph, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to create bedGraph: {result.stderr}")
            return False
        
        cmd_sizes = f"{samtools_path} faidx {contigs_fa} && cut -f1,2 {contigs_fa}.fai > {contigs_fa}.sizes"
        result = subprocess.run(cmd_sizes, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to create sizes file: {result.stderr}")
            return False
        
        cmd_bw = f"{bedgraph_to_bigwig_path} {output_bam}.coverage.bedgraph {contigs_fa}.sizes {output_bw}"
        result = subprocess.run(cmd_bw, shell=True, capture_output=True, text=True)
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






def main():
    parser = argparse.ArgumentParser(description="Pathogen Detection Pipeline")
    parser.add_argument("--diamond", required=True, help="Diamond results TSV file")
    parser.add_argument("--fq1", required=True, help="Original forward reads (FASTQ)")
    parser.add_argument("--fq2", required=True, help="Original reverse reads (FASTQ)")
    parser.add_argument("--megahit", required=True, help="Megahit results directory")
    parser.add_argument("--taxids", required=True, help="Target taxid list file")
    parser.add_argument("--virus-index", required=True, help="BLASTN virus database path")
    parser.add_argument("--output", required=True, help="Output directory")
    
    # Diamond 过滤参数
    parser.add_argument("--pident", type=float, default=60, help="Diamond identity threshold (%%)")
    parser.add_argument("--qcov", type=float, default=50, help="Diamond coverage threshold (%%)")
    parser.add_argument("--evalue", type=float, default=1e-10, help="Diamond e-value threshold")
    
    # BLASTN 参数
    parser.add_argument("--blastn-identity", type=float, default=70, help="BLASTN identity threshold (%%)")
    parser.add_argument("--blastn-evalue", type=float, default=1e-4, help="BLASTN e-value threshold")
    parser.add_argument("--blastn-coverage", type=float, default=70, help="BLASTN coverage threshold (%%)")
    
    # 其他参数
    parser.add_argument("--threads", type=int, default=64, help="Number of threads")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose/debug logging")
    parser.add_argument("--no-resume", action="store_true", help="Ignore checkpoint and start fresh")
    
    parser.add_argument("--bowtie2-bin", default="bowtie2", help="bowtie2 executable path")
    parser.add_argument("--samtools-bin", default="samtools", help="samtools executable path")
    parser.add_argument("--bedtools-bin", default="bedtools", help="bedtools executable path")
    parser.add_argument("--bedgraph-to-bigwig-bin", default="bedGraphToBigWig", help="bedGraphToBigWig executable path")
    parser.add_argument("--blastn-bin", default="blastn", help="blastn executable path")
    
    args = parser.parse_args()
    
    # Create output directory FIRST (before logger tries to create log file)
    os.makedirs(args.output, exist_ok=True)
    
    # Setup logger (log file can now be created in existing directory)
    log_file = os.path.join(args.output, f"pathogen_detect_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    # 设置日志级别
    logger = setup_logger("PathogenDetect", log_file)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.info("Verbose logging enabled")
    
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
        with open(step1_taxids_cache, 'w') as f:
            json.dump({"taxids": list(taxids), "taxid_names": taxid_names}, f, indent=2)
        
        checkpoints["completed_steps"].append("step1_load_taxids")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 1 (completed)")
        with open(step1_taxids_cache, 'r') as f:
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
    
    # Step 2: Filter diamond results
    step2_filtered_cache = os.path.join(args.output, "step2_filtered_contigs.json")
    
    if not should_skip_step("step2_filter", checkpoints, [step2_filtered_cache]):
        logger.info("Step 2: Filtering diamond results...")
        logger.info(f"Filters: pident>={args.pident}%, qcov>={args.qcov}%, evalue<={args.evalue}")
        filtered_contigs = filter_diamond(args.diamond, taxids, args.pident, args.qcov, args.evalue, logger)
        logger.info(f"Found {len(filtered_contigs)} virus-related contigs after filtering")
        
        # Save filtered_contigs for checkpoint resume
        with open(step2_filtered_cache, 'w') as f:
            json.dump(filtered_contigs, f, indent=2)
        
        checkpoints["completed_steps"].append("step2_filter")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 2 (completed)")
        with open(step2_filtered_cache, 'r') as f:
            filtered_contigs = json.load(f)
        logger.info(f"Loaded {len(filtered_contigs)} filtered contigs from cache")
    
    # 详细报告
    if filtered_contigs:
        logger.info("-" * 50)
        logger.info("Filtered contigs summary:")
        logger.info(f"  Total contigs: {len(filtered_contigs)}")
        
        # 统计覆盖度和相似度分布
        pidents = [info['pident'] for info in filtered_contigs.values()]
        qcovs = [info['qcov'] for info in filtered_contigs.values()]
        
        logger.info(f"  Identity range: {min(pidents):.1f}% - {max(pidents):.1f}% (avg: {sum(pidents)/len(pidents):.1f}%)")
        logger.info(f"  Coverage range: {min(qcovs):.1f}% - {max(qcovs):.1f}% (avg: {sum(qcovs)/len(qcovs):.1f}%)")
        
        # 按病毒种类分组统计
        virus_groups = defaultdict(list)
        for contig_id, info in filtered_contigs.items():
            virus_groups[info['sscinames']].append(contig_id)
        
        logger.info(f"  Virus groups: {len(virus_groups)}")
        for virus_name, contig_list in sorted(virus_groups.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
            logger.info(f"    - {virus_name}: {len(contig_list)} contigs")
        
        logger.info("-" * 50)
    
    if not filtered_contigs:
        logger.warning("No virus-related contigs found!")
        # Write empty results
        output_tsv = os.path.join(args.output, "detection_results.tsv")
        output_summary = os.path.join(args.output, "detection_summary.txt")
        with open(output_tsv, 'w') as f:
            f.write("Virus_Seq\tBest_Hit_Seq\tCoverage(%)\tDepth\tContig_ID\n")
        with open(output_summary, 'w') as f:
            f.write("No pathogen detected.\n")
        return
    
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
        from collections import defaultdict as dd
        contig_ids_set = set(filtered_contigs.keys())
        extracted_count = 0
        with open(virus_contigs_fa, 'w') as out:
            write_seq = False
            with open(megahit_contigs_fa, 'r') as f:
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
            filter_blastn_results(blastn_raw, blastn_filtered, 
                                 args.blastn_identity, args.blastn_coverage, logger)
            checkpoints["completed_steps"].append("step4_blastn")
            save_checkpoints(checkpoints, checkpoint_file)
        else:
            logger.error("BLASTN failed!")
            sys.exit(1)
    else:
        logger.info("Skipping Step 4 (completed)")
    
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
    
    if not should_skip_step("step6_report", checkpoints, [output_tsv, output_summary]):
        logger.info("Step 6: Generating detection report...")
        
        # Parse BLASTN results and select best hit for each contig
        blastn_hits = defaultdict(list)
        with open(blastn_filtered, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 13:
                    continue
                contig_id = fields[0]
                virus_seq = fields[1]
                bitscore = float(fields[11])
                blastn_hits[contig_id].append({
                    'virus_seq': virus_seq,
                    'bitscore': bitscore,
                    'sseqid': virus_seq
                })
        
        # Select best hit for each contig (highest bitscore)
        best_hits = {}
        for contig_id, hits in blastn_hits.items():
            best_hit = max(hits, key=lambda x: x['bitscore'])
            best_hits[contig_id] = best_hit
        
        # Parse coverage summary
        coverage_stats = {}
        with open(coverage_summary, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    coverage_stats[fields[0]] = {
                        'avg_depth': float(fields[1]),
                        'min_depth': float(fields[2]),
                        'max_depth': float(fields[3]),
                        'length': int(fields[4])
                    }
        
        # Generate report (sorted by virus sequence)
        results = []
        for contig_id in sorted(best_hits.keys(), key=lambda x: best_hits[x]['virus_seq']):
            hit = best_hits[contig_id]
            cov = coverage_stats.get(contig_id, {})
            results.append({
                'Virus_Seq': hit['virus_seq'],
                'Best_Hit_Seq': hit['sseqid'],
                'Coverage': f"{cov.get('avg_depth', 0):.2f}",
                'Depth': f"{cov.get('avg_depth', 0):.2f}",
                'Contig_ID': contig_id
            })
        
        # Write TSV report
        with open(output_tsv, 'w', encoding='utf-8') as f:
            f.write('Virus_Seq\tBest_Hit_Seq\tCoverage(%)\tDepth\tContig_ID\n')
            for r in results:
                f.write(f"{r['Virus_Seq']}\t{r['Best_Hit_Seq']}\t{r['Coverage']}\t{r['Depth']}\t{r['Contig_ID']}\n")
        
        # Write summary
        with open(output_summary, 'w', encoding='utf-8') as f:
            f.write("Pathogen Detection Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Total pathogens detected: {len(results)}\n\n")
            
            # Group by virus sequence
            virus_groups = defaultdict(list)
            for r in results:
                virus_groups[r['Virus_Seq']].append(r['Contig_ID'])
            
            for virus_seq, contig_ids in sorted(virus_groups.items()):
                f.write(f"- {virus_seq}\n")
                f.write(f"  Contigs: {len(contig_ids)}\n")
                if len(contig_ids) <= 10:
                    f.write(f"  IDs: {', '.join(contig_ids)}\n\n")
                else:
                    f.write(f"  IDs: {', '.join(contig_ids[:10])}... (showing 10 of {len(contig_ids)})\n\n")
        
        checkpoints["completed_steps"].append("step6_report")
        save_checkpoints(checkpoints, checkpoint_file)
    else:
        logger.info("Skipping Step 6 (completed)")
    
    # Pipeline completed
    total_time = time.time() - start_time
    logger.info("=" * 60)
    logger.info(f"Pipeline completed in {total_time:.2f} seconds")
    logger.info(f"Detected {len(results) if 'results' in dir() else 0} pathogens")
    logger.info(f"Results saved to: {args.output}")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
