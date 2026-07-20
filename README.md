# hvFinder - mNGS Pipeline and GUI Tool

[![Version](https://img.shields.io/badge/version-0.0.1-blue.svg)](https://github.com/your-repo/hvFinder)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux-lightgrey.svg)]()

**hvFinder** is a comprehensive mNGS (metagenomic Next-Generation Sequencing) analysis pipeline with a user-friendly GUI interface. It provides automated pathogen detection from raw sequencing data to final reporting. This project is licensed under the [MIT License](https://choosealicense.com/licenses/mit/).

---

## ✨ Features

### Local GUI (Windows)
- 🖥️ Intuitive PySide6-based GUI
- 🔐 SSH connection management
- 📁 Remote project file management
- 📊 Real-time log display

### Remote mNGS Analysis Pipeline
- 🧬 Quality control (fastp)
- 🎯 Host read removal (bowtie2)
- 🧹 rRNA removal (ribodetector)
- 🔗 De novo assembly (megahit)
- 📌 Clustering (cd-hit)
- 🏷️ Taxonomic annotation (diamond)

### Pathogen Detection
- 🦠 Target virus identification
- 🔬 BLASTN alignment
- 📈 Coverage calculation (samtools)
- 📊 Detailed reporting with BigWig visualization

### Advanced Features
- ⏸️ Checkpoint resume support
- 🚀 Slurm job scheduler support
- ⚙️ Configurable parameters
- 📝 Verbose logging for debugging

---

## 🚀 Local GUI
We recommend using the compiled hvFinder.exe.

**Download**
1. Download `hvFinder.exe` from the [releases](https://github.com/yingHH/hvFinder/releases) page
2. Place it in a directory of your choice

**First Run**
1. Double-click `hvFinder.exe` to launch the GUI
2. Configure SSH connection in **Tab 1**
3. Create or select a project directory
4. Configure analysis parameters in **Tab 2**
5. Click **🚀 RUN PIPELINE** to start analysis

**File Structure**
```
your-directory/
├── hvFinder.exe      # Main executable
└── settings.ini      # Configuration file (auto-created)
```

If you want to manually compile the GUI from the source code.
**Requirements**
- Python 3.9+
- Windows 10/11 or Linux

**Step 1: Clone Repository**
```bash
git clone https://github.com/yingHH/hvFinder.git
cd hvFinder
```

**Step 2: Create Conda Environment**
```bash
conda create -n hvfinder python=3.9
conda activate hvfinder
```

**Step 3: Install Dependencies**
```bash
pip install PySide6 paramiko PyYAML cryptography
```

**Step 4: Launch GUI**
```bash
# Windows
python hvFinder/gui.py

# Or use the entry point
hvFinder-gui
```

**Step 5: (Optional) Build EXE**
```bash
pip install pyinstaller pyinstaller-hooks-contrib
python -m PyInstaller --clean --name hvFinder \
    --windowed --onefile \
    --add-data="settings.ini:." \
    --add-data="configs:configs" \
    --hidden-import paramiko \
    --hidden-import PySide6 \
    --hidden-import yaml \
    --hidden-import cryptography \
    --icon=hvFinder.ico \
    hvFinder/gui.py
```

---

## 🖥️ Remote Server Setup

The remote Linux server must have the following tools installed by conda environments (essential):

**Core Tools**
```bash
# Sequence processing
fastp
bowtie2
samtools

# rRNA removal
ribodetector

# Assembly
megahit
cd-hit

# Annotation
diamond
blastn

# Visualization
bedtools
bedGraphToBigWig
```

**Installation (Conda)**
```bash
conda install -c bioconda fastp bowtie2 samtools ribodetector megahit cd-hit diamond blastn bedtools ucsc-bedgraphtobigwig
```

**Python Package**
Install hvFinder on remote server:
```bash
cd /path/to/hvFinder
```

**Databases**
The constructed databases are listed in [figshare](https://doi.org/10.6084/m9.figshare.32984498). Download the **virus2taxid** file and the **Highly pathetic virus database folder**. For host genome removal with Bowtie2, you may build a corresponding bowtie2 index by `bowtie2-build` command.

---

## 📖 Usage

#### Tab 1: SSH & Files

Configure SSH connection to remote Linux server.

**Settings:**
- **Host IP**: Remote server address
- **Port**: SSH port (default: 22)
- **Username**: SSH username
- **Password**: SSH password
- **Base Directory**: Remote base directory for projects

**Actions:**
- **List Files**: Browse remote directory
- **Create Project**: Create new project directory

#### Tab 2: Pre-processing

Configure and run the mNGS analysis pipeline.

**Input Data:**
- **Reads R1**: Forward reads path (remote)
- **Reads R2**: Reverse reads path (remote)

**Execution Mode:**
- **Direct Python**: Run directly with Python
- **Slurm (srun)**: Submit to Slurm cluster

**Advanced Parameters:**
- fastp parameters
- bowtie2 parameters
- megahit parameters
- cd-hit parameters
- diamond parameters

**Settings:**
- ☑️ Save Settings
- ☐ Verbose Logging
- ☑️ Resume from checkpoint

#### Tab 3: Pathogen Detection

Detect specific pathogens from assembled contigs.

**Input Settings:**
- **TaxID List File**: Target virus taxids (remote path)
- **Virus Index Path**: BLASTN database path (remote)

**Execution Mode:**
- **Direct Python** / **Slurm (srun)**

**Filter Parameters (Diamond):**
- **Identity (%)**: Minimum identity (default: 60)
- **Coverage (%)**: Minimum coverage (default: 50)
- **E-value**: Maximum E-value (default: 1e-10)

**BLASTN Parameters:**
- **Identity (%)**: BLASTN identity threshold (default: 70)
- **E-value**: BLASTN E-value threshold (default: 1e-4)
- **Coverage (%)**: BLASTN coverage threshold (default: 70)

**Output Settings:**
- **Folder Name**: Output directory name (default: pathogen_results)

**Actions:**
- **💾 View Results**: View detection results in table
- **💾 Export CSV**: Export results to CSV file

---

## ⚙️ Configuration

#### settings.ini

Configuration file is stored in the same directory as the executable (or in the project directory for development).
There is an example setting.ini in the root directory.

**Example:**
```ini
[SSH]
host = ssh.example.com
port = 22
user = username
pwd = password
base = /home/data/mNGS_projects
script_base = /path/to/hvFinder

[PIPELINE]
fq1 = /path/to/reads_1.fq
fq2 = /path/to/reads_2.fq
run_mode = Direct Python
verbose = False
resume = True

[PATHOGEN]
taxid_file = /path/to/taxids.txt
virus_index = /path/to/virus_db
run_mode = Direct Python
output_folder = pathogen_results
diamond_pident = 60
diamond_qcov = 50
diamond_evalue = 1e-10
blastn_identity = 70
blastn_evalue = 1e-4
blastn_coverage = 70
verbose = False
resume = True
```

---

## ⏸️ Checkpoint Resume

Both Tab 2 and Tab 3 support checkpoint resume functionality.

### How It Works
1. After each major step, a checkpoint file is created
2. If the pipeline is interrupted, restart will skip completed steps
3. Checkpoint files:
   - Tab 2: `checkpoint_status.pkl`
   - Tab 3: `pathogen_detect.checkpoints.json`

### Usage
- ✅ **Resume from checkpoint (if exists)**: Automatically detect and resume
- ❌ Uncheck to force restart from beginning

### Manual Reset
Delete checkpoint files to force fresh start:
```bash
# Tab 2
rm /path/to/project/preprocess/checkpoint_status.pkl

# Tab 3
rm /path/to/project/pathogen_results/pathogen_detect.checkpoints.json
```

---

## 📊 Output Files

#### Tab 2 (Pre-processing)
```
project_name/
└── preprocess/
    ├── results.fastp.html      # Quality control report
    ├── results.fastp.json      # QC statistics
    ├── results.bowtie2.log     # Host removal log
    ├── results.diamond.tsv     # Taxonomic annotation
    └── results_megahit_result/ # Assembly results
        ├── results.contigs.fa
        └── ...
```

#### Tab 3 (Pathogen Detection)
```
project_name/
└── pathogen_results/
    ├── blastn_results.tsv           # BLASTN alignments
    ├── contigs_mapping.bam          # Reads to contigs alignment
    ├── contigs_mapping.bam.bai      # BAM index
    ├── contig_coverage_summary.tsv  # Coverage statistics
    ├── contig_coverage.bw           # BigWig coverage track
    ├── detection_results.tsv        # Final report
    └── detection_summary.txt        # Summary report
```

#### Detection Report Format
| Column | Description |
|--------|-------------|
| Virus_Name | Virus name |
| TaxID | NCBI Taxonomy ID |
| Best_Hit_Seq | Best hit virus sequence name |
| Coverage(%) | Average coverage percentage |
| Depth | Average sequencing depth |
| Contig_ID | Contig identifier |

---

## ❓ FAQ

**Q: hvFinder.exe fails to start**  
**A:** Make sure you have downloaded the complete executable. Antivirus software may block it - try adding an exception.

**Q: SSH connection fails**  
**A:** Check network connection, SSH credentials, and firewall settings. Test connection with:
```bash
ssh username@host -p port
```

**Q: Pipeline fails at specific step**  
**A:** 
1. Enable **Verbose Logging** in Settings
2. Check remote server logs
3. Verify all required tools are installed on remote server
4. Check disk space on remote server

**Q: Checkpoint resume not working**  
**A:** 
1. Verify checkpoint file exists
2. Ensure **Resume from checkpoint** is checked
3. Delete checkpoint file and restart if needed

**Q: No pathogens detected**  
**A:** 
1. Verify taxid list file is correct
2. Check Diamond results contain target viruses
3. Adjust filtering thresholds (lower identity/coverage)
4. Verify BLASTN database is correctly built

**Q: BigWig file not generated**  
**A:** Ensure `bedGraphToBigWig` is installed on remote server:
```bash
conda install -c bioconda ucsc-bedgraphtobigwig
```

---

## 🙏 Acknowledgments

- **fastp**: Ultra-fast all-in-one FASTQ preprocessor
- **bowtie2**: Fast gapped-read alignment
- **megahit**: Efficient metagenomics assembler
- **diamond**: Fast protein alignment
- **PySide6**: Python Qt6 bindings
- **PyInstaller**: Python to executable converter

---
