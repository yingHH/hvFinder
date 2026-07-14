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
import traceback
from datetime import datetime
from collections import namedtuple

# 模拟 logger 导入，如果没有外部库请确保 logger.py 存在
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

# 扩展命名元组，包含时间消耗
StepResult = namedtuple('StepResult', ['success', 'stdout', 'stderr', 'duration'])

class MNGSPipeline:
    def __init__(self, config_path, output_prefix):
        self.start_time = time.time()
        self.output_prefix = output_prefix
        
        # 1. 加载配置
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                self.cfg = yaml.safe_load(f)
        except Exception as e:
            print(f"FATAL: Could not load config file {config_path}: {e}")
            sys.exit(1)

        # 2. 初始化日志
        log_file = f"{output_prefix}_pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        self.logger = setup_logger("mNGS_Core", log_file)
        self.logger.info("="*60)
        self.logger.info(f"Pipeline Session Started: {datetime.now().isoformat()}")
        self.logger.info(f"Using Config: {config_path}")
        self.logger.info(f"Output Prefix: {self.output_prefix}")
        self.logger.info(f"Process ID: {os.getpid()}")
        self.logger.info("="*60)

    def _execute_command(self, step_name, cmd):
        """核心执行器：记录命令、监控耗时、捕获异常"""
        cmd_clean = cmd.strip()
        self.logger.info(f"Executing Step [{step_name}]...")
        self.logger.debug(f"Command String:\n{cmd_clean}")
        
        start_step = time.time()
        try:
            proc = subprocess.run(
                cmd_clean, 
                shell=True, 
                capture_output=True, 
                text=True, 
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
        """全方位的环境与输入文件预检"""
        self.logger.info("--- [Phase: Pre-Check] Validating Environment ---")
        validation_failed = False

        # 1. 检查输入文件
        for f in [fq1, fq2]:
            if not os.path.exists(f):
                self.logger.error(f"Input file missing: {f}")
                validation_failed = True

        # 2. 检查软件路径
        required_tools = ['fastp', 'bowtie2', 'samtools', 'ribodetector', 'megahit', 'cdhit', 'diamond']
        for tool in required_tools:
            path = os.path.expanduser(self.cfg[tool]['bin'])
            if not (shutil.which(path) or os.path.exists(path)):
                self.logger.error(f"Software path invalid: [{tool}] -> {path}")
                validation_failed = True
            else:
                self.logger.info(f"Verified tool: {tool:15} at {path}")

        # 3. 检查数据库
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
        cmd = f"""
        {c['bin']} \\
            --thread {t} \\
            --in1 {fq1} --in2 {fq2} \\
            --out1 {out_f[0]} --out2 {out_f[1]} \\
            --json {self.output_prefix}.fastp.json \\
            --html {self.output_prefix}.fastp.html
        """
        return self._execute_command("fastp_qc", cmd)

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
            
        cmd = f"""
        {c['bin']} \\
            --num-cpu-threads {t} \\
            -1 {in_f[0]} -2 {in_f[1]} \\
            --out-dir {out_dir} \\
            --out-prefix {self.output_prefix}
        """
        res = self._execute_command("assembly_megahit", cmd)
        if res.success:
            for f in in_f: 
                if os.path.exists(f): os.remove(f)
        return res

    def _run_cdhit(self, megahit_dir, out_fa, t):
        c = self.cfg['cdhit']
        in_fa = os.path.join(megahit_dir, f"{self.output_prefix}.contigs.fa")
        if not os.path.exists(in_fa):
            in_fa = os.path.join(megahit_dir, "final.contigs.fa")
            
        cmd = f"""
        {c['bin']} \\
            -T {t} \\
            -M {c['params']['memory_limit']} \\
            -i {in_fa} \\
            -o {out_fa}
        """
        return self._execute_command("cdhit_dedup", cmd)

    def _run_diamond(self, in_fa, out_txt, t):
        c = self.cfg['diamond']
        cmd = f"""
        {c['bin']} blastx \\
            --threads {t} \\
            --query {in_fa} \\
            --db {c['db_path']} \\
            --out {out_txt} \\
            {c['params']['sensitivity']} \\
            --evalue {c['params']['evalue']} \\
            --max-target-seqs {c['params']['max_target_seqs']} \\
            --outfmt {c['params']['outfmt']}
        """
        res = self._execute_command("taxonomy_diamond", cmd)
        if res.success and os.path.exists(in_fa):
            os.remove(in_fa)
        return res

    def start_pipeline(self, fq1, fq2, threads=None, no_resume=False):
        """断点续算主循环"""
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
            with open(status_path, 'r') as f:
                checkpoint = json.load(f)
            self.logger.info(f"Checkpoint detected. Resuming. Previously finished: {list(checkpoint['completed_steps'].keys())}")
        
        f_qc = [os.path.join(tmp_dir, "step1_qc_R1.fq.gz"), os.path.join(tmp_dir, "step1_qc_R2.fq.gz")]
        f_non_host = [os.path.join(tmp_dir, "step2_nonhost_R1.fq.gz"), os.path.join(tmp_dir, "step2_nonhost_R2.fq.gz")]
        f_non_ribo = [os.path.join(tmp_dir, "step3_noribo_R1.fq.gz"), os.path.join(tmp_dir, "step3_noribo_R2.fq.gz")]
        f_cdhit = os.path.join(tmp_dir, "step5_clustered.fasta")
        mega_dir = f"{self.output_prefix}_megahit_result"
        final_tsv = f"{self.output_prefix}.diamond.tsv"

        pipeline_map = [
            ("FASTP", self._run_fastp, (fq1, fq2, f_qc, t)),
            ("BOWTIE2", self._run_bowtie2, (f_qc, f_non_host, t)),
            ("RIBODETECTOR", self._run_ribodetector, (f_non_host, f_non_ribo, t)),
            ("MEGAHIT", self._run_megahit, (f_non_ribo, mega_dir, t)),
            ("CDHIT", self._run_cdhit, (mega_dir, f_cdhit, t)),
            ("DIAMOND", self._run_diamond, (f_cdhit, final_tsv, t))
        ]

        for name, func, args in pipeline_map:
            if name in checkpoint['completed_steps']:
                self.logger.info(f"Skipping Step [{name}] (Already done).")
                
                # 断点续算时，如果跳过 MEGAHIT 但 clean_reads 不存在，尝试保存
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
            
            # 在 MEGAHIT 执行之前保存 clean_reads（因为 MEGAHIT 会删除输入文件）
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
                with open(status_path, 'w') as f:
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

    args = parser.parse_args()
    
    # 配置文件定位逻辑
    config_file = args.config
    
    # 获取脚本所在目录和输出文件所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.dirname(os.path.abspath(args.output))
    
    # 如果没有指定 -c 参数，则自动搜寻
    if not config_file:
        # 1. 尝试使用输出目录中的 config.yaml
        local_output_config = os.path.join(output_dir, "config.yaml")
        if os.path.exists(local_output_config):
            config_file = local_output_config
        else:
            # 2. 尝试从脚本所在目录的 configs/config2.yaml 拷贝
            template_config = os.path.join(script_dir, "configs", "config2.yaml")
            if os.path.exists(template_config):
                # 确保输出目录存在
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                # 拷贝并更名为 config.yaml
                config_file = local_output_config
                shutil.copy2(template_config, config_file)
                print(f"INFO: No config found in output dir. Copied template from {template_config} to {config_file}")
            else:
                print(f"FATAL: Config file not provided and default template '{template_config}' not found.")
                sys.exit(1)
    
    # 执行流水线
    pipe = MNGSPipeline(config_file, args.output)
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