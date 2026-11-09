# -*- coding: utf-8 -*-
"""
Author: Ying Huang
Date: 2026-01-21
Description: Enhanced mNGS Pipeline with Host Removal (Bowtie2), structural logging, 
             performance tracking, and robust fault-tolerant checkpointing.
"""

import os
import sys
import time
import subprocess
import argparse
import yaml
import json
import shutil
import glob
import traceback
from datetime import datetime
from collections import namedtuple

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

StepResult = namedtuple('StepResult', ['success', 'stdout', 'stderr', 'duration'])

FASTP_MAX_THREADS = 16



def resume_frontier(step_order, completed_steps):
    """
    Last completed step - the only step whose artifacts resume validates.
    """
    done = [s for s in step_order if s in (completed_steps or {})]
    return done[-1] if done else None


def first_empty_paths(paths):
    """
    Return the entries of ``paths`` that are missing or zero-sized.
    """
    bad = []
    for p in paths or []:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            bad.append(p)
    return bad


class MNGSPipeline:
    def __init__(self, config_path, output_prefix):
        self.start_time = time.time()
        self.output_prefix = os.path.abspath(output_prefix)
        _out_dir = os.path.dirname(self.output_prefix)
        if _out_dir:
            os.makedirs(_out_dir, exist_ok=True)
        
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                self.cfg = yaml.safe_load(f)
        except Exception as e:
            print(f"FATAL: Could not load config file {config_path}: {e}")
            sys.exit(1)

        log_file = f"{output_prefix}_pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        self.logger = setup_logger("mNGS_Core", log_file)

        ribo_cfg = self.cfg.get('ribodetector', {}) or {}
        self.skip_ribo = (ribo_cfg.get('enabled', True) is False)

        self.logger.info("="*60)
        self.logger.info(f"Pipeline Session Started: {datetime.now().isoformat()}")
        self.logger.info(f"Using Config: {config_path}")
        self.logger.info(f"Output Prefix: {self.output_prefix}")
        self.logger.info(f"Process ID: {os.getpid()}")
        self.logger.info("="*60)

    def _execute_command(self, step_name, cmd):
        """
        Core runner: echo the command, time it, and capture failures.
        """
        cmd_clean = cmd.strip()
        self.logger.info(f"Executing Step [{step_name}]...")
        self.logger.debug(f"Command String:\n{cmd_clean}")
        
        start_step = time.time()
        try:
            proc = subprocess.run(
                cmd_clean, 
                shell=True, 
                capture_output=True, 
                encoding='utf-8', errors='replace', 
                executable='/bin/bash' 
            )
            duration = time.time() - start_step
            
            if proc.returncode == 0:
                self.logger.info(f"Step [{step_name}] finished successfully in {duration:.2f}s")
                return StepResult(True, proc.stdout, proc.stderr, duration)
            else:
                self.logger.error(f"Step [{step_name}] FAILED (Exit Code: {proc.returncode})")
                self.logger.error(f"Error Output:\n{proc.stderr}")
                return StepResult(False, proc.stdout, proc.stderr, duration)
                
        except Exception as e:
            self.logger.critical(f"Exception during Step [{step_name}]: {str(e)}")
            self.logger.error(traceback.format_exc())
            return StepResult(False, "", str(e), 0)

    def pre_check(self, fq1, fq2):
        """
        Pre-flight check of the environment and of the input files.
        """
        self.logger.info("--- [Phase: Pre-Check] Validating Environment ---")
        validation_failed = False

        for f in [fq1, fq2]:
            if not os.path.exists(f):
                self.logger.error(f"Input file missing: {f}")
                validation_failed = True

        required_tools = ['fastp', 'bowtie2', 'samtools', 'megahit', 'cdhit', 'diamond']
        if not getattr(self, 'skip_ribo', False):
            required_tools.insert(3, 'ribodetector')
        else:
            self.logger.info("rRNA removal disabled (skip_ribo=True); ribodetector not required")
        for tool in required_tools:
            path = os.path.expanduser(self.cfg[tool]['bin'])
            if not (shutil.which(path) or os.path.exists(path)):
                self.logger.error(f"Software path invalid: [{tool}] -> {path}")
                validation_failed = True
            else:
                self.logger.info(f"Verified tool: {tool:15} at {path}")

        bt2_index = os.path.expanduser(self.cfg['bowtie2']['params']['index_path'])
        if not (os.path.exists(f"{bt2_index}.1.bt2") or os.path.exists(f"{bt2_index}.1.bt2l")):
            self.logger.error(f"Bowtie2 Index missing or prefix invalid: {bt2_index}")
            validation_failed = True

        dm_db = os.path.expanduser(self.cfg['diamond']['db_path'])
        if not os.path.exists(dm_db):
            self.logger.error(f"Diamond Database missing: {dm_db}")
            validation_failed = True

        if validation_failed:
            self.logger.critical("Pre-check failed. Please rectify config.yaml or input paths.")
            raise RuntimeError("Environment validation failed.")
        
        self.logger.info("--- [Phase: Pre-Check] All systems GO ---\n")

    def _run_fastp(self, fq1, fq2, out_f, t):
        c = self.cfg['fastp']
        fastp_threads = max(1, min(int(t), FASTP_MAX_THREADS))
        if int(t) > FASTP_MAX_THREADS:
            self.logger.info(f"fastp threads capped at {FASTP_MAX_THREADS} "
                             f"(requested {t})")
        cmd = f"""
        {c['bin']} \\
            --thread {fastp_threads} \\
            --in1 {fq1} --in2 {fq2} \\
            --out1 {out_f[0]} --out2 {out_f[1]} \\
            --json {self.output_prefix}.fastp.json \\
            --html {self.output_prefix}.fastp.html
        """
        res = self._execute_command("fastp_qc", cmd)
        if res.success:
            res = self._check_fastp_reads(res)
        return res

    def _check_fastp_reads(self, res):
        """
        Zero surviving reads must fail at QC time.

        fastp exits 0 even when it filters everything out; continuing would open a
        multi-hundred-GB DIAMOND database with an empty FASTA and then fail with an
        opaque `Error detecting input file format`.
        """
        json_path = f"{self.output_prefix}.fastp.json"
        try:
            with open(json_path, encoding="utf-8") as f:
                summary = (json.load(f) or {}).get("summary", {})
            before = int(summary.get("before_filtering", {})
                         .get("total_reads", 0) or 0)
            after = int(summary.get("after_filtering", {})
                        .get("total_reads", 0) or 0)
        except (OSError, ValueError, TypeError, AttributeError):
            self.logger.warning(
                f"fastp JSON not readable ({json_path}); empty-read check skipped")
            return res
        if after > 0:
            return res
        msg = (f"No reads passed QC: input_total_reads={before}, "
               f"after_filtering=0. Empty or too-short FASTQ; "
               f"aborting before host removal / assembly.")
        self.logger.error(msg)
        return StepResult(False, res.stdout, msg, res.duration)

    def _run_bowtie2(self, in_f, out_f, t):
        c = self.cfg['bowtie2']
        s = self.cfg['samtools']
        index_path = c['params']['index_path']
        un_prefix = out_f[0].replace("_R1.fq.gz", "")
        
        cmd = f"""
        {c['bin']} \\
            -x {index_path} \\
            -1 {in_f[0]} -2 {in_f[1]} \\
            --threads {t} \\
            --very-sensitive \\
            --un-conc-gz {un_prefix}_R%.fq.gz \\
            2> {self.output_prefix}.bowtie2.log | \\
        {s['bin']} view -@ {t} -bS - > /dev/null
        """
        res = self._execute_command("host_removal_bowtie2", cmd)
        if res.success:
            for f in in_f:
                if os.path.exists(f): os.remove(f)
            self.logger.info("Host removal completed. Non-host reads preserved.")
        return res

    def _run_ribodetector(self, in_f, out_f, t):
        c = self.cfg['ribodetector']
        cmd = f"""
        {c['bin']} \\
            --threads {t} \\
            --input {in_f[0]} {in_f[1]} \\
            --output {out_f[0]} {out_f[1]} \\
            --len {c['params']['len']} \\
            -e {c['params']['enzyme']} \\
            --chunk_size {c['params']['chunk_size']}
        """
        res = self._execute_command("ribo_removal", cmd)
        if res.success:
            for f in in_f: 
                if os.path.exists(f): os.remove(f)
        return res

    def _run_megahit(self, in_f, out_dir, t):
        c = self.cfg['megahit']
        if os.path.exists(out_dir) and c['params'].get('force'):
            shutil.rmtree(out_dir)

        out_prefix = os.path.basename(self.output_prefix)

        cmd = f"""
        {c['bin']} \\
            --num-cpu-threads {t} \\
            -1 {in_f[0]} -2 {in_f[1]} \\
            --out-dir {out_dir} \\
            --out-prefix {out_prefix}
        """
        res = self._execute_command("assembly_megahit", cmd)
        if res.success:
            for f in in_f: 
                if os.path.exists(f): os.remove(f)
        return res

    def _run_cdhit(self, megahit_dir, out_fa, t):
        c = self.cfg['cdhit']
        in_fa = self._locate_contigs(megahit_dir)
        if in_fa is None:
            self.logger.error(f"No contigs found in {megahit_dir}; cannot run CD-HIT")
            return StepResult(False, "", f"no contigs.fa in {megahit_dir}", 0.0)

        cmd = f"""
        {c['bin']} \\
            -T {t} \\
            -M {c['params']['memory_limit']} \\
            -i {in_fa} \\
            -o {out_fa}
        """
        return self._execute_command("cdhit_dedup", cmd)

    def _locate_contigs(self, megahit_dir):
        """
        Locate the MEGAHIT assembly, tolerating different --out-prefix values and
        older naming.
        """
        base = os.path.basename(self.output_prefix)
        candidates = [
            os.path.join(megahit_dir, f"{base}.contigs.fa"),
            os.path.join(megahit_dir, "final.contigs.fa"),
        ]
        for p in candidates:
            if os.path.exists(p):
                return p
        found = sorted(glob.glob(os.path.join(megahit_dir, "*.contigs.fa")))
        if found:
            return found[0]
        for p in (f"{self.output_prefix}.contigs.fa",):
            if os.path.exists(p):
                return p
        return None

    def _prescreen_cfg(self):
        """
        Return (db_path, bin, params) for the pre-screen, which is always run.
        """
        c = self.cfg.get('prescreen') or {}
        params = c.get('params') or {}
        db = getattr(self, 'prescreen_db', None) or c.get('db_path') or ''
        bin_ = c.get('bin') or self.cfg['diamond']['bin']
        return os.path.expanduser(str(db)), bin_, params

    def _run_prescreen(self, in_fa, cand_fa, t):
        """
        Stage 1: one blastx pass against the high pathogenic virus protein database,
        keeping only the contigs that hit.

        Outputs under preprocess/: `<prefix>.prescreen.tsv` (hits),
        `<prefix>.prescreen.json` (counts and timing) and
        `prescreen_candidates.fasta` (the BLASTN input of stage 2).
        """
        db, bin_, params = self._prescreen_cfg()
        if not db or not os.path.exists(db):
            return StepResult(False, "", f"prescreen DB not found: {db}", 0.0)

        evalue = params.get('evalue', '1e-3')
        sens = params.get('sensitivity', '--sensitive')
        k = int(params.get('max_target_seqs', 1) or 1)
        hits = f"{self.output_prefix}.prescreen.tsv"

        cmd = f"""
        {bin_} blastx \\
            --threads {t} \\
            --query {in_fa} \\
            --db {db} \\
            --out {hits} \\
            {sens} \\
            --evalue {evalue} \\
            --max-target-seqs {k} \\
            --outfmt 6 qseqid sseqid pident evalue
        """
        self.logger.info(
            "PRESCREEN: db=%s evalue=%s k=%d", db, evalue, k)
        res = self._execute_command("prescreen_diamond", cmd)
        if not res.success:
            return res

        keep = set()
        if os.path.exists(hits):
            with open(hits, 'r', encoding='utf-8') as fh:
                for line in fh:
                    qid = line.split('\t', 1)[0].strip()
                    if qid:
                        keep.add(qid)

        n_in = n_out = 0
        with open(in_fa, 'r', encoding='utf-8') as fh, \
                open(cand_fa, 'w', encoding='utf-8') as out:
            write = False
            for line in fh:
                if line.startswith('>'):
                    n_in += 1
                    write = line[1:].split()[0] in keep
                    if write:
                        n_out += 1
                if write:
                    out.write(line)

        stats = {
            'db': db,
            'contigs_in': n_in,
            'contigs_kept': n_out,
            'hits': len(keep),
            'evalue': str(evalue),
            'max_target_seqs': k,
            'wall_s': round(res.duration, 2),
        }
        with open(f"{self.output_prefix}.prescreen.json", 'w', encoding='utf-8') as fh:
            json.dump(stats, fh, indent=2, ensure_ascii=False)
        self.logger.info(
            "PRESCREEN: kept %d/%d contig(s) in %.1fs; candidates written for the "
            "detection stage (BLASTN -> nr DIAMOND, tab3)",
            n_out, n_in, res.duration)
        return res

    def start_pipeline(self, fq1, fq2, threads=None, no_resume=False):
        """
        Main loop, with checkpoint-based resume.
        """
        t = threads or self.cfg['global']['threads']
        tmp_dir = f"{self.output_prefix}{self.cfg['global']['tmp_dir_suffix']}"
        status_path = os.path.join(os.path.dirname(self.output_prefix), "preprocess.checkpoints.json")
        
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        checkpoint = {'completed_steps': {}, 'performance': {}}
        if no_resume and os.path.exists(status_path):
            os.remove(status_path)
            self.logger.info("Removed old checkpoint (--no-resume specified)")
        
        if os.path.exists(status_path) and not no_resume:
            try:
                with open(status_path, 'r', encoding='utf-8') as f:
                    loaded = json.load(f)
                if not isinstance(loaded, dict) or \
                        not isinstance(loaded.get('completed_steps'), dict):
                    raise ValueError("unexpected checkpoint structure")
                checkpoint = loaded
                checkpoint.setdefault('performance', {})
                self.logger.info(
                    "Checkpoint detected. Resuming. Previously finished: "
                    f"{list(checkpoint['completed_steps'].keys())}")
            except (OSError, ValueError) as exc:
                self.logger.warning(
                    f"Checkpoint file is unreadable or corrupt ({exc}); "
                    "ignoring it and starting from scratch.")
                try:
                    os.replace(status_path, status_path + ".corrupt")
                except OSError:
                    pass
                checkpoint = {'completed_steps': {}, 'performance': {}}
        
        f_qc = [os.path.join(tmp_dir, "step1_qc_R1.fq.gz"), os.path.join(tmp_dir, "step1_qc_R2.fq.gz")]
        f_non_host = [os.path.join(tmp_dir, "step2_nonhost_R1.fq.gz"), os.path.join(tmp_dir, "step2_nonhost_R2.fq.gz")]
        if getattr(self, 'skip_ribo', False):
            f_non_ribo = f_non_host
            self.logger.info("rRNA removal skipped: MEGAHIT will consume host-depleted reads directly")
        else:
            f_non_ribo = [os.path.join(tmp_dir, "step3_noribo_R1.fq.gz"), os.path.join(tmp_dir, "step3_noribo_R2.fq.gz")]
        f_cdhit = os.path.join(tmp_dir, "step5_clustered.fasta")
        mega_dir = f"{self.output_prefix}_megahit_result"
        f_cand = os.path.join(os.path.dirname(self.output_prefix),
                              "prescreen_candidates.fasta")

        pipeline_map = [
            ("FASTP", self._run_fastp, (fq1, fq2, f_qc, t)),
            ("BOWTIE2", self._run_bowtie2, (f_qc, f_non_host, t)),
        ]
        if not getattr(self, 'skip_ribo', False):
            pipeline_map.append(("RIBODETECTOR", self._run_ribodetector, (f_non_host, f_non_ribo, t)))
        pipeline_map.extend([
            ("MEGAHIT", self._run_megahit, (f_non_ribo, mega_dir, t)),
            ("CDHIT", self._run_cdhit, (mega_dir, f_cdhit, t)),
            ("PRESCREEN", self._run_prescreen, (f_cdhit, f_cand, t)),
        ])

        def _step_outputs(step):
            if step == "FASTP":
                return list(f_qc)
            if step == "BOWTIE2":
                return list(f_non_host)
            if step == "RIBODETECTOR":
                return list(f_non_ribo)
            if step == "MEGAHIT":
                if not os.path.isdir(mega_dir):
                    return [mega_dir]
                found = [os.path.join(mega_dir, f) for f in os.listdir(mega_dir)
                         if f.endswith(".contigs.fa")]
                return found or [os.path.join(mega_dir, "*.contigs.fa")]
            if step == "CDHIT":
                return [f_cdhit]
            if step == "PRESCREEN":
                return [f_cand]
            return []

        def _clear_step_artifacts(step):
            """
            Drop the previous artifacts of this step (MEGAHIT refuses an existing
            output directory).
            """
            if step == "MEGAHIT" and os.path.isdir(mega_dir):
                shutil.rmtree(mega_dir, ignore_errors=True)
            elif step == "CDHIT" and os.path.exists(f_cdhit):
                os.remove(f_cdhit)
            elif step == "PRESCREEN":
                for p in (f_cand,
                          f"{self.output_prefix}.prescreen.tsv",
                          f"{self.output_prefix}.prescreen.json"):
                    if os.path.exists(p):
                        os.remove(p)

        step_order = [n for n, _f, _a in pipeline_map]

        frontier = resume_frontier(step_order, checkpoint['completed_steps'])

        for name, func, args in pipeline_map:
            if name in checkpoint['completed_steps']:
                stale = (first_empty_paths(_step_outputs(name))
                         if name == frontier else [])
                if stale:
                    self.logger.warning(
                        f"Checkpoint lists Step [{name}] as done, but output is "
                        f"missing/empty ({', '.join(stale)}). Invalidating it and "
                        f"re-running this and all downstream steps.")
                    _clear_step_artifacts(name)
                    for downstream in step_order[step_order.index(name):]:
                        checkpoint['completed_steps'].pop(downstream, None)
                else:
                    self.logger.info(f"Skipping Step [{name}] (Already done).")

                    if name == "MEGAHIT":
                        clean_reads_dir = os.path.join(os.path.dirname(self.output_prefix), "clean_reads")
                        clean_r1 = os.path.join(clean_reads_dir, "clean_R1.fq.gz")
                        clean_r2 = os.path.join(clean_reads_dir, "clean_R2.fq.gz")
                        if not os.path.exists(clean_r1) or not os.path.exists(clean_r2):
                            os.makedirs(clean_reads_dir, exist_ok=True)
                            if os.path.exists(f_non_ribo[0]) and os.path.exists(f_non_ribo[1]):
                                shutil.copy(f_non_ribo[0], clean_r1)
                                shutil.copy(f_non_ribo[1], clean_r2)
                                self.logger.info(f"Clean reads saved to: {clean_reads_dir} (from checkpoint resume)")
                            else:
                                self.logger.warning(f"Clean reads files not found in tmp_dir, cannot save clean_reads")
                    continue
            
            if name == "MEGAHIT":
                clean_reads_dir = os.path.join(os.path.dirname(self.output_prefix), "clean_reads")
                os.makedirs(clean_reads_dir, exist_ok=True)
                if os.path.exists(f_non_ribo[0]) and os.path.exists(f_non_ribo[1]):
                    shutil.copy(f_non_ribo[0], os.path.join(clean_reads_dir, "clean_R1.fq.gz"))
                    shutil.copy(f_non_ribo[1], os.path.join(clean_reads_dir, "clean_R2.fq.gz"))
                    self.logger.info(f"Clean reads saved to: {clean_reads_dir}")
            
            res = func(*args)
            if res.success:
                checkpoint['completed_steps'][name] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                checkpoint['performance'][name] = f"{res.duration:.2f}s"
                with open(status_path, 'w', encoding='utf-8') as f:
                    json.dump(checkpoint, f, indent=2)
            else:
                self.logger.critical(f"Pipeline halted at Step [{name}]. Manual intervention required.")
                sys.exit(1)

        total_time = time.time() - self.start_time
        self.logger.info("="*60)
        self.logger.info("PIPELINE SUMMARY")
        for step, p_time in checkpoint['performance'].items():
            self.logger.info(f"- {step:15}: {p_time}")
        self.logger.info(f"Total Execution Time: {total_time/60:.2f} minutes")
        
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        self.logger.info(f"Checkpoint saved to: {status_path}")
        self.logger.info("Cleanup temporary files. Pipeline finished successfully.")
        self.logger.info("="*60)

def main():
    parser = argparse.ArgumentParser(description="Professional mNGS Workflow (Fastp -> Bowtie2 -> RiboDetector -> Megahit -> CD-HIT -> Diamond)")
    parser.add_argument("-1", "--fq1", required=True, help="Forward reads")
    parser.add_argument("-2", "--fq2", required=True, help="Reverse reads")
    parser.add_argument("-o", "--output", required=True, help="Output prefix/directory for results")
    parser.add_argument("-c", "--config", help="Path to config.yaml (optional)")
    parser.add_argument("-t", "--threads", type=int, help="Override CPU threads")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--no-resume", action="store_true", help="Disable checkpoint resume and delete old checkpoint")
    parser.add_argument("--skip-ribo", action="store_true",
                        help="Skip rRNA removal (ribodetector); useful for DNA-only libraries")
    parser.add_argument("--node-mem-gb", type=float, default=None,
                        help="Deprecated: nr DIAMOND moved to the detection "
                             "stage (pathogen_detect.py --node-mem-gb). Accepted and "
                             "ignored for backward compatibility")
    parser.add_argument("--prescreen-db", default=None,
                        help="Protein DIAMOND DB for the prescreen stage "
                             "(overrides config prescreen.db_path)")

    args = parser.parse_args()
    
    config_file = args.config
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.dirname(os.path.abspath(args.output))
    
    if not config_file:
        local_output_config = os.path.join(output_dir, "config.yaml")
        if os.path.exists(local_output_config):
            config_file = local_output_config
        else:
            template_config = os.path.join(script_dir, "configs", "config2.yaml")
            if os.path.exists(template_config):
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                config_file = local_output_config
                shutil.copy2(template_config, config_file)
                print(f"INFO: No config found in output dir. Copied template from {template_config} to {config_file}")
            else:
                print(f"FATAL: Config file not provided and default template '{template_config}' not found.")
                sys.exit(1)
    
    pipe = MNGSPipeline(config_file, args.output)
    if args.skip_ribo:
        pipe.skip_ribo = True
        pipe.logger.info("rRNA removal will be skipped (--skip-ribo)")
    if args.node_mem_gb is not None:
        pipe.logger.warning(
            "--node-mem-gb is deprecated here; nr DIAMOND now runs in "
            "pathogen_detect.py (pass --node-mem-gb there): value %s GB ignored",
            args.node_mem_gb)
    if getattr(args, 'prescreen_db', None):
        pipe.prescreen_db = args.prescreen_db
    try:
        pipe.pre_check(args.fq1, args.fq2)
        pipe.start_pipeline(args.fq1, args.fq2, args.threads, no_resume=args.no_resume)
    except KeyboardInterrupt:
        pipe.logger.warning("Pipeline interrupted by user.")
    except Exception as e:
        if hasattr(pipe, 'logger'):
            pipe.logger.error(f"Unexpected Crash: {str(e)}")
            pipe.logger.error(traceback.format_exc())
        else:
            print(f"Unexpected Crash before logger init: {str(e)}")
            traceback.print_exc()

if __name__ == '__main__':
    main()