# -*- coding: utf-8 -*-

"""
Author: Ying Huang
Date: 2026-01-27
Description: hvFinder Pro 
"""

import os
os.environ["QT_LOGGING_RULES"] = "*.debug=false;qt.qpa.xcb=false"

import sys
import shlex
import shutil
import time
import yaml
import configparser
from datetime import datetime
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QTextEdit, QTabWidget, QGroupBox, 
    QFormLayout, QScrollArea, QCheckBox, QMessageBox, QInputDialog,
    QComboBox, QTableWidget, QTableWidgetItem, QDialog, QFileDialog,
    QSpinBox
)
from PySide6.QtCore import Signal, QThread, QTimer
from PySide6.QtGui import QTextCursor


def _paramiko():
    """
    Import paramiko lazily.

    paramiko drags in cryptography (with a native extension on Windows), which
    costs 0.5-1.0 s at import time. Importing it here keeps GUI start-up fast;
    the module is still imported only once per process.
    """
    import paramiko
    return paramiko


def parse_script_base(text):
    """
    Extract the remote hvFinder installation directory from
    `get_script_base` output (ANSI codes and chatter are ignored).
    """
    import re
    clean = re.sub(r"\x1b\[[0-9;]*m", "", text or "")
    lines = [line.strip() for line in clean.strip().splitlines() if line.strip()]
    for line in reversed(lines):
        if "hvFinder" in line and line.startswith("/"):
            return line
    return ""


# ---------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------
_PARAM_KEY_LABELS = {
    'bin': 'Bin',
    'db_path': 'DB Path',
    'index_path': 'Index Path',
    'sensitivity': 'Sensitivity',
    'evalue': 'E-value',
    'max_target_seqs': 'Max Target Seqs',
    'outfmt': 'Outfmt',
    'memory_limit': 'Memory Limit',
    'chunk_size': 'Chunk Size',
    'gb_per_block': 'GB per Block',
    'mem_reserve_gb': 'Mem Reserve (GB)',
    'index_chunks': 'Index Chunks',
    'node_mem_gb': 'Node Mem (GB)',
    'tmp_dir_suffix': 'Tmp Dir Suffix',
    'threads': 'Threads',
}

_PARAM_ROW_LABELS = {
    ('prescreen', 'bin'): 'Diamond Bin',
    ('diamond', 'bin'): 'Diamond Bin',
    ('prescreen', 'db_path'): 'DB Path',
    ('diamond', 'db_path'): 'DB Path',
    ('diamond', 'params', 'blocks'): 'Blocks (-b, 0 = Auto)',
}

_PARAM_HINTS = {
    ('prescreen', 'bin'):
        'Diamond bin. Leave empty = reuse the nr DIAMOND bin from Tab 3.',
    ('prescreen', 'db_path'):
        'Screening database: high-consequence virus protein subset.',
    ('diamond', 'bin'):
        'Diamond bin. Leave empty = use "diamond" from PATH.',
    ('diamond', 'db_path'):
        'Comprehensive nr protein database used by the nr DIAMOND step.',
}


def param_label(map_key):
    """
    Render a config2.yaml key path as a form label (Title Case, acronyms kept).
    """
    map_key = tuple(map_key)
    if map_key in _PARAM_ROW_LABELS:
        return _PARAM_ROW_LABELS[map_key]
    key = map_key[-1] if map_key else ''
    if key in _PARAM_KEY_LABELS:
        return _PARAM_KEY_LABELS[key]
    return ' '.join(word.capitalize() for word in key.replace('_', ' ').split())


def parse_detection_results(content):
    """
    Parse `detection_results.tsv` into (rows, has_taxid).

    Two headers are supported: the legacy layout with a TaxID column and the
    current one without it. The split is decided by the header only, never by
    the column count; empty input or a header-only file yields an empty list.
    """
    results = []
    lines = (content or "").strip().splitlines()
    if len(lines) <= 1:
        return results, False
    has_taxid = "TaxID" in lines[0]
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.strip().split('\t')
        if has_taxid:
            if len(fields) >= 6:
                results.append({
                    'Virus_Seq': fields[0],
                    'TaxID': fields[1],
                    'Best_Hit_Seq': fields[2],
                    'Coverage': fields[3],
                    'Depth': fields[4],
                    'Contig_ID': fields[5],
                })
        elif len(fields) >= 5:
            results.append({
                'Virus_Seq': fields[0],
                'Best_Hit_Seq': fields[1],
                'Coverage': fields[2],
                'Depth': fields[3],
                'Contig_ID': fields[4],
            })
    return results, has_taxid


def app_settings_path():
    """
    Absolute path of settings.ini.

    Running from source: next to this file. Frozen build: next to the exe
    (a bundled template is copied out of _MEIPASS on first start).
    """
    if getattr(sys, "frozen", False):
        base = os.path.dirname(sys.executable)
    else:
        base = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base, "settings.ini")
    if getattr(sys, "frozen", False) and not os.path.exists(path):
        bundled = os.path.join(getattr(sys, "_MEIPASS", ""), "settings.ini")
        if os.path.exists(bundled):
            try:
                shutil.copy2(bundled, path)
            except OSError:
                pass
    return path

def read_ini_file(path):
    """
    Read an INI file as UTF-8 (a BOM is tolerated).

    configparser defaults to the locale encoding on Windows, which would raise
    UnicodeDecodeError on UTF-8 content, so the encoding is always explicit.
    """
    parser = configparser.ConfigParser(interpolation=None)
    try:
        with open(path, "r", encoding="utf-8-sig") as f:
            parser.read_file(f)
    except (OSError, configparser.Error, UnicodeDecodeError):
        return configparser.ConfigParser(interpolation=None)
    return parser


def pipeline_template_path():
    """
    Path of the config2.yaml template (frozen: _MEIPASS/configs).
    """
    if getattr(sys, "frozen", False):
        return os.path.join(getattr(sys, "_MEIPASS", ""), "configs", "config2.yaml")
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(here, "..", "server", "configs", "config2.yaml"))


class SSHWorker(QThread):
    log_signal = Signal(str)
    status_signal = Signal(int)
    finished_signal = Signal(bool, str) 

    def __init__(self, host, port, user, pwd, command=None, upload_cfg=None, target_dir=None, task_type="cmd", detached=False, upload_files=None):
        super().__init__()
        self.host = host
        self.port = port
        self.user = user
        self.pwd = pwd
        self.command = command
        self.upload_cfg = upload_cfg
        self.upload_files = upload_files or {}
        self.target_dir = target_dir
        self.task_type = task_type
        self.detached = detached

    def run(self):
        ssh = None
        try:
            pk = _paramiko()
            ssh = pk.SSHClient()
            ssh.set_missing_host_key_policy(pk.AutoAddPolicy())
            ssh.connect(self.host, port=self.port, username=self.user, password=self.pwd, timeout=10)
            
            self.status_signal.emit(1)
            self.emit_log(f"✅ Connected to {self.host}")

            if self.task_type == "mkdir":
                ssh.exec_command(f"mkdir -p {shlex.quote(self.target_dir)}")
                self.emit_log(f"📁 Created workspace: {self.target_dir}")
                self.finished_signal.emit(True, "CREATED")
                return

            if self.target_dir:
                stdin, stdout, stderr = ssh.exec_command(
                    f"mkdir -p {shlex.quote(self.target_dir)}")
                stdout.channel.recv_exit_status()
                
                if self.upload_cfg:
                    sftp = ssh.open_sftp()
                    remote_cfg_path = f"{self.target_dir}/config.yaml"
                    with sftp.file(remote_cfg_path, "w") as f:
                        f.write(self.upload_cfg)
                    sftp.close()
                    self.emit_log(f"📄 Config uploaded to: {remote_cfg_path}")

                if self.upload_files:
                    sftp = ssh.open_sftp()
                    for name, content in self.upload_files.items():
                        remote_path = f"{self.target_dir}/{name}"
                        with sftp.file(remote_path, "w") as f:
                            f.write(content)
                        self.emit_log(f"📄 Uploaded: {remote_path}")
                    sftp.close()

            if self.command:
                self.status_signal.emit(2)
                if self.detached:
                    final_cmd = self._detached_wrapper()
                else:
                    final_cmd = (f"cd {shlex.quote(self.target_dir)} && {self.command}"
                                if self.target_dir else self.command)
                self.emit_log(f"🚀 Running command: {self.command}")
                
                stdin, stdout, stderr = ssh.exec_command(
                    f"bash -l -c {shlex.quote(final_cmd)} 2>&1", get_pty=True)
                for line in stdout:
                    self.log_signal.emit(line.strip())
                
                stdout.channel.recv_exit_status()
                self.status_signal.emit(1)
            
            self.finished_signal.emit(True, "SUCCESS")

        except Exception as e:
            self.emit_log(f"❌ Error: {str(e)}")
            self.status_signal.emit(0)
            self.finished_signal.emit(False, str(e))
        finally:
            if ssh:
                ssh.close()

    def emit_log(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_signal.emit(f"[{timestamp}] {message}")

    RUN_PID_NAME = ".hvfinder.run.pid"
    RUN_LOG_NAME = ".hvfinder.run.log"
    RUN_JOB_NAME = ".hvfinder.slurm.jobid"
    RUN_RC_NAME = ".hvfinder.run.rc"
    JOB_SCRIPT_NAME = ".hvfinder.job.sh"
    SLURM_OUT_NAME = ".hvfinder.slurm.out"

    def _detached_wrapper(self):
        """
        Run the command in a new setsid session and report its exit code.

        The remote job therefore survives a dropped SSH connection; the log is
        followed with `tail --pid`, so it ends when the job does. The exit code is
        written to `.hvfinder.run.rc` and echoed as `[hvfinder] exit=N`, which is
        what lets the GUI tell success from failure.
        """
        workdir = self.target_dir or "."
        inner = (f"cd {shlex.quote(self.target_dir)} && {self.command}"
                 if self.target_dir else self.command)
        pid = f"{workdir}/{self.RUN_PID_NAME}"
        log = f"{workdir}/{self.RUN_LOG_NAME}"
        rc = f"{workdir}/{self.RUN_RC_NAME}"
        job = f"{workdir}/{self.RUN_JOB_NAME}"
        launcher = (f"echo $$ > {shlex.quote(pid)}; "
                    f"bash -l -c {shlex.quote(inner)} "
                    f"> {shlex.quote(log)} 2>&1 < /dev/null; "
                    f"rc=$?; echo $rc > {shlex.quote(rc)}; exit $rc")
        return (
            f"mkdir -p {shlex.quote(workdir)} && "
            f"rm -f {shlex.quote(pid)} {shlex.quote(job)} {shlex.quote(rc)} && "
            f"setsid bash -c {shlex.quote(launcher)} & "
            f"sleep 1; P=$(cat {shlex.quote(pid)} 2>/dev/null); "
            f"echo '[hvfinder] detached, pid='\"$P\"', log={log}'; "
            f"if [ -n \"$P\" ]; then tail -n +1 -s 1 --pid=\"$P\" -f {shlex.quote(log)}; "
            f"else tail -n +1 -f {shlex.quote(log)}; fi; "
            f"wait; RC=$(cat {shlex.quote(rc)} 2>/dev/null || echo '?'); "
            f"echo '[hvfinder] exit='\"$RC\""
        )


class SSHProbe(QThread):
    """
    Background probe of the saved SSH credentials (connect and close).
    """
    result_signal = Signal(bool, str)

    def __init__(self, host, port, user, pwd, timeout=8):
        super().__init__()
        self.host = host
        self.port = port
        self.user = user
        self.pwd = pwd
        self.timeout = timeout

    def run(self):
        try:
            pk = _paramiko()
            ssh = pk.SSHClient()
            ssh.set_missing_host_key_policy(pk.AutoAddPolicy())
            ssh.connect(
                self.host, port=self.port, username=self.user, password=self.pwd,
                timeout=self.timeout, banner_timeout=self.timeout,
                auth_timeout=self.timeout, look_for_keys=False, allow_agent=False,
            )
            ssh.close()
            self.result_signal.emit(True, "connected")
        except Exception as exc:  # noqa: BLE001
            self.result_signal.emit(False, f"{type(exc).__name__}: {exc}")


class RemoteQuery(QThread):
    """
    Background execution of one remote command; the result comes back by signal.
    """
    result_signal = Signal(bool, str)

    def __init__(self, host, port, user, pwd, command, timeout=10):
        super().__init__()
        self.host = host
        self.port = port
        self.user = user
        self.pwd = pwd
        self.command = command
        self.timeout = timeout
        self.context = None

    def run(self):
        try:
            pk = _paramiko()
            ssh = pk.SSHClient()
            ssh.set_missing_host_key_policy(pk.AutoAddPolicy())
            ssh.connect(
                self.host, port=self.port, username=self.user, password=self.pwd,
                timeout=self.timeout, banner_timeout=self.timeout,
                auth_timeout=self.timeout, look_for_keys=False, allow_agent=False,
            )
            try:
                _, stdout, _ = ssh.exec_command(self.command, timeout=self.timeout)
                out = stdout.read().decode("utf-8", "replace")
            finally:
                ssh.close()
            self.result_signal.emit(True, out)
        except Exception as exc:  # noqa: BLE001
            self.result_signal.emit(False, f"{type(exc).__name__}: {exc}")


# --- Results Dialog ---
class ResultsDialog(QDialog):
    def __init__(self, results, output_dir, has_taxid=False, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Pathogen Detection Results")
        self.resize(900, 500)
        self.results = results
        self.output_dir = output_dir
        self.has_taxid = has_taxid
        
        layout = QVBoxLayout(self)
        
        # Results table
        self.table = QTableWidget()
        if has_taxid:
            # Old format with TaxID
            self.table.setColumnCount(6)
            self.table.setHorizontalHeaderLabels([
                "Virus Seq", "TaxID", "Best Hit Seq", "Coverage(%)", "Depth", "Contig ID"
            ])
        else:
            # New format without TaxID
            self.table.setColumnCount(5)
            self.table.setHorizontalHeaderLabels([
                "Virus Seq", "Best Hit Seq", "Coverage(%)", "Depth", "Contig ID"
            ])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setAlternatingRowColors(True)
        
        # Populate table
        for row, r in enumerate(results):
            self.table.insertRow(row)
            self.table.setItem(row, 0, QTableWidgetItem(str(r.get('Virus_Seq', ''))))
            if has_taxid:
                self.table.setItem(row, 1, QTableWidgetItem(str(r.get('TaxID', ''))))
                self.table.setItem(row, 2, QTableWidgetItem(str(r.get('Best_Hit_Seq', ''))))
                self.table.setItem(row, 3, QTableWidgetItem(str(r.get('Coverage', ''))))
                self.table.setItem(row, 4, QTableWidgetItem(str(r.get('Depth', ''))))
                self.table.setItem(row, 5, QTableWidgetItem(str(r.get('Contig_ID', ''))))
            else:
                self.table.setItem(row, 1, QTableWidgetItem(str(r.get('Best_Hit_Seq', ''))))
                self.table.setItem(row, 2, QTableWidgetItem(str(r.get('Coverage', ''))))
                self.table.setItem(row, 3, QTableWidgetItem(str(r.get('Depth', ''))))
                self.table.setItem(row, 4, QTableWidgetItem(str(r.get('Contig_ID', ''))))
        
        layout.addWidget(self.table)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.export_btn = QPushButton("💾 Export CSV")
        self.export_btn.clicked.connect(self.export_csv)
        
        self.ok_btn = QPushButton("OK")
        self.ok_btn.clicked.connect(self.accept)
        
        btn_layout.addWidget(self.export_btn)
        btn_layout.addStretch()
        btn_layout.addWidget(self.ok_btn)
        layout.addLayout(btn_layout)
    
    def export_csv(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Results", "pathogen_detection_results.csv", "CSV Files (*.csv)"
        )
        if file_path:
            import csv
            with open(file_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                if self.has_taxid:
                    writer.writerow(["Virus Seq", "TaxID", "Best Hit Seq", "Coverage(%)", "Depth", "Contig ID"])
                    for r in self.results:
                        writer.writerow([
                            r.get('Virus_Seq', ''), r.get('TaxID', ''), r.get('Best_Hit_Seq', ''),
                            r.get('Coverage', ''), r.get('Depth', ''), r.get('Contig_ID', '')
                        ])
                else:
                    writer.writerow(["Virus Seq", "Best Hit Seq", "Coverage(%)", "Depth", "Contig ID"])
                    for r in self.results:
                        writer.writerow([
                            r.get('Virus_Seq', ''), r.get('Best_Hit_Seq', ''),
                            r.get('Coverage', ''), r.get('Depth', ''), r.get('Contig_ID', '')
                        ])
            QMessageBox.information(self, "Success", f"Results exported to:\n{file_path}")

PRESET_VIRUS_INDEXES = {
    "Test Database": "data/test_idx/test",
    "hvn (high pathogenic virus nucleotide DB)": "data/hvn/hvn.fa",
    "Custom": ""
}

PRESET_TAXID_FILES = {
    "Default Virus TaxIDs": "data/default_virus_taxids.txt",
    "High Pathogenic Virus TaxIDs": "data/virus_taxids_high_pathogenic.txt",
    "Custom": ""
}

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("hvFinder Pro - Management Console")
        self.resize(1000, 850)
        
        self.settings_file = app_settings_path()
        
        config_path = pipeline_template_path()
        self.pipeline_cfg = self.load_local_template(config_path)
        
        self.current_project_dir = ""
        self._script_base = ""
        self._slurm_jobs = {}
        self._poll_timers = {}
        self._poll_state = {}         # target_dir -> {offset, file, task_type, ...}
        self._task_rc = {}
        self._task_started = {}       # id(log_widget) -> time.monotonic()
        
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        self.init_ssh_tab()
        self.init_pipeline_tab()
        self.init_pathogen_tab()

        self._loading = False
        self.load_saved_settings()
        self.link_execution_controls()
        self.setup_autosave()
        self.apply_button_styles()
        self._sync_execution_rows()

        self.tabs.currentChanged.connect(self.on_tab_changed)

        QTimer.singleShot(0, self.try_auto_connect)

    def apply_button_styles(self):
        """
        Only recolour buttons; keep the native shape.
        """
        blue_style = """
            QPushButton {
                background-color: #2E86C1;
                color: white;
                border: 1px solid #2E86C1;
            }
            QPushButton:hover {
                background-color: #3498DB;
                border: 1px solid #3498DB;
            }
            QPushButton:pressed {
                background-color: #1A5276;
                border: 1px solid #1A5276;
            }
            QPushButton:disabled {
                background-color: #D5D8DC;
                color: #909497;
                border: 1px solid #BDC3C7;
            }
        """
        self.list_dir_btn.setStyleSheet(blue_style)
        self.select_project_btn.setStyleSheet(blue_style)
        self.toggle_btn.setStyleSheet(blue_style)
        self.run_btn.setStyleSheet(blue_style)
        
        # Pathogen tab buttons
        self.pathogen_toggle_btn.setStyleSheet(blue_style)
        self.pathogen_run_btn.setStyleSheet(blue_style)
        self.view_results_btn.setStyleSheet(blue_style)
        self.export_csv_btn.setStyleSheet(blue_style)

    def load_local_template(self, pth):
        """
        Load a YAML template from the given path.
        """
        if os.path.exists(pth):
            try:
                with open(pth, "r", encoding="utf-8") as f:
                    return yaml.safe_load(f) or {}
            except Exception as e:
                print(f"Error loading YAML: {e}")
        return {}

    def init_ssh_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setSpacing(10)
        layout.setContentsMargins(10, 10, 10, 18)

        auth_group = QGroupBox("Server Authentication")
        auth_form = QFormLayout(auth_group)
        self.host_edit = QLineEdit()
        self.port_edit = QLineEdit("22")
        self.user_edit = QLineEdit()
        self.pwd_edit = QLineEdit(); self.pwd_edit.setEchoMode(QLineEdit.Password)
        auth_form.addRow("Host IP:", self.host_edit)
        auth_form.addRow("Port:", self.port_edit)
        auth_form.addRow("Username:", self.user_edit)
        auth_form.addRow("Password:", self.pwd_edit)

        self.remember_cb = QCheckBox("Save Settings")
        self.remember_cb.setChecked(True)
        auth_form.addRow("", self.remember_cb)

        layout.addWidget(auth_group)

        dir_group = QGroupBox("Project Root & Project Selection")
        dir_form = QFormLayout(dir_group)
        
        self.search_dir_edit = QLineEdit("")
        self.search_dir_edit.setPlaceholderText(
            "Root directory holding the projects, e.g. /data/mNGS_projects")
        self.search_dir_edit.setToolTip(
            "Root directory that holds all projects (used as the output path prefix).\n"
            "Actual output: <root>/<project>/<step>/results/")
        
        btn_layout = QHBoxLayout()
        self.select_project_btn = QPushButton("📁 Select / New Project…")
        self.select_project_btn.setStyleSheet("background-color: #D4E6F1; font-weight: bold;")
        self.select_project_btn.setToolTip(
            "List existing projects under this root (one level only), or create a new one")
        self.select_project_btn.clicked.connect(self.on_select_project)
        self.list_dir_btn = QPushButton("🔍 List Files")
        self.list_dir_btn.setToolTip(
            "List files and sub-folders one level below the Active Project")
        self.list_dir_btn.clicked.connect(self.on_list_remote)
        self.list_dir_btn.setEnabled(False)
        
        btn_layout.addWidget(self.select_project_btn)
        btn_layout.addWidget(self.list_dir_btn)

        self.sync_info_label = QLabel("Active Project: None")
        self.sync_info_label.setStyleSheet("color: #2980B9; font-style: bold;")

        dir_form.addRow("Project Root (prefix):", self.search_dir_edit)
        dir_form.addRow(btn_layout)
        dir_form.addRow(self.sync_info_label)
        layout.addWidget(dir_group)

        status_layout = QHBoxLayout()
        self.status_label = QLabel("● Disconnected")
        self.status_label.setStyleSheet("color: #7F8C8D; font-weight: bold;")
        status_layout.addStretch()
        status_layout.addWidget(self.status_label)
        layout.addLayout(status_layout)

        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setStyleSheet(
            "background-color: #FFFFFF; color: #000000; font-family: 'Consolas';")
        layout.addWidget(self.log_text, stretch=1)
        
        self.tabs.addTab(tab, "1. SSH & Files")

    def init_pipeline_tab(self):
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setSpacing(10)
        layout.setContentsMargins(10, 10, 10, 18)

        # Project info display (read-only) + run status
        top_row = QHBoxLayout()
        self.pipeline_project_label = QLabel("📁 Project: None")
        self.pipeline_project_label.setStyleSheet("color: #2980B9; font-weight: bold;")
        self.pipe_status_label = QLabel("● Idle")
        self.pipe_status_label.setStyleSheet("color: #7F8C8D; font-weight: bold;")
        top_row.addWidget(self.pipeline_project_label)
        top_row.addStretch()
        top_row.addWidget(self.pipe_status_label)
        layout.addLayout(top_row)

        io_group = QGroupBox("Input Data (Remote Paths)")
        io_form = QFormLayout(io_group)
        self.fq1_edit = QLineEdit()
        self.fq2_edit = QLineEdit()
        io_form.addRow("Reads R1:", self.fq1_edit)
        io_form.addRow("Reads R2:", self.fq2_edit)
        layout.addWidget(io_group)

        exec_group = QGroupBox("Execution Mode")
        exec_form = QFormLayout(exec_group)
        self.run_mode_combo = QComboBox()
        self.run_mode_combo.addItems(["Direct Python", "Slurm (srun)"])
        exec_form.addRow("Run Mode:", self.run_mode_combo)

        self.cores_spin = QSpinBox()
        self.cores_spin.setRange(1, 256)
        self.cores_spin.setValue(int(self.pipeline_cfg.get('global', {}).get('threads', 64) or 64))
        self.cores_spin.setToolTip(
            "CPU cores: passed to mNGS_pip -t in Direct mode, and also used as "
            "srun --cpus-per-task in Slurm mode.\n"
            "Note: fastp uses at most 16 threads; extra cores are ignored.")
        exec_form.addRow("CPU Cores:", self.cores_spin)

        self.slurm_extra_edit = QLineEdit()
        self.slurm_extra_edit.setPlaceholderText("Example: --exclusive")
        self.slurm_extra_edit.setToolTip("Appended to the srun command in Slurm mode only")
        exec_form.addRow("Extra Slurm Args:", self.slurm_extra_edit)
        self._slurm_rows = [
            (exec_form.labelForField(self.slurm_extra_edit), self.slurm_extra_edit)]
        layout.addWidget(exec_group)

        self.toggle_btn = QPushButton("Advanced Parameters ▼")
        self.toggle_btn.setCheckable(True)
        layout.addWidget(self.toggle_btn)

        self.adv_widget = QWidget()
        adv_vbox = QVBoxLayout(self.adv_widget)
        adv_vbox.setContentsMargins(0, 0, 0, 0)
        self.ui_param_map = {}

        skip_sections = {'diamond'}
        for section, data in self.pipeline_cfg.items():
            if section == 'global' and 'threads' not in data: continue
            if section in skip_sections: continue
            
            sec_group = QGroupBox(f"Setting: {section}")
            sec_form = QFormLayout(sec_group)
            
            def add_rows(item_data, prefix_keys=[]):
                for k, v in item_data.items():
                    if isinstance(v, dict):
                        add_rows(v, prefix_keys + [k])
                    else:
                        map_key = tuple([section] + prefix_keys + [k])
                        edit = QLineEdit(str(v))
                        hint = _PARAM_HINTS.get(map_key)
                        if hint:
                            edit.setPlaceholderText(hint)
                        sec_form.addRow(f"{param_label(map_key)}:", edit)
                        self.ui_param_map[map_key] = edit

            if isinstance(data, dict):
                add_rows(data)
                adv_vbox.addWidget(sec_group)

        self.adv_scroll = QScrollArea()
        self.adv_scroll.setWidget(self.adv_widget)
        self.adv_scroll.setWidgetResizable(True)
        self.adv_scroll.setMinimumHeight(200)
        self.adv_scroll.setMaximumHeight(300)
        self.adv_scroll.setVisible(False)
        layout.addWidget(self.adv_scroll)
        self.toggle_btn.clicked.connect(
            lambda: self.adv_scroll.setVisible(self.toggle_btn.isChecked()))

        settings_group = QGroupBox("Settings")
        settings_row = QHBoxLayout(settings_group)
        self.save_pipeline_cb = QCheckBox("Save Settings")
        self.save_pipeline_cb.setChecked(True)
        self.pipeline_verbose_cb = QCheckBox("Verbose Logging")
        self.pipeline_verbose_cb.setChecked(False)
        self.pipeline_resume_cb = QCheckBox("Resume from checkpoint (if exists)")
        self.pipeline_resume_cb.setChecked(True)
        settings_row.addWidget(self.save_pipeline_cb)
        settings_row.addWidget(self.pipeline_verbose_cb)
        settings_row.addWidget(self.pipeline_resume_cb)
        self.sync_db_btn = QPushButton("⇪ Sync DB to Server")
        self.sync_db_btn.setToolTip(
            "Write the database paths (diamond db / bowtie2 index / BLASTN "
            "virus index / taxid list) from Settings into the server-side "
            "configs/config2.yaml.")
        self.sync_db_btn.clicked.connect(
            lambda: self.sync_db_to_server(self.pipe_log_text))
        settings_row.addWidget(self.sync_db_btn)
        settings_row.addStretch()
        layout.addWidget(settings_group)

        self.run_btn = QPushButton("🚀 RUN PIPELINE")
        self.run_btn.clicked.connect(lambda: self.on_launch(step_name="preprocess"))
        layout.addWidget(self.run_btn)

        self.stop_pipe_btn = QPushButton("⏸ STOP PIPELINE")
        self.stop_pipe_btn.setToolTip("Are you sure to terminate the running task?")
        self.stop_pipe_btn.clicked.connect(self.on_stop_pipeline)
        layout.addWidget(self.stop_pipe_btn)

        self.pipe_log_text = QTextEdit()
        self.pipe_log_text.setReadOnly(True)
        self.pipe_log_text.setPlaceholderText(
            "Not started. Select a project and R1/R2, then click RUN PIPELINE.")
        self.pipe_log_text.setStyleSheet(
            "background-color: #FFFFFF; color: #000000; font-family: 'Consolas';")
        self.pipe_log_text.setMinimumHeight(140)
        layout.addWidget(self.pipe_log_text, stretch=1)

        self.tabs.addTab(content, "2. Pre-processing")

    def init_pathogen_tab(self):
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setSpacing(10)
        layout.setContentsMargins(10, 10, 10, 18)
        
        # Project info display (read-only) + run status
        top_row = QHBoxLayout()
        self.pathogen_project_label = QLabel("📁 Project: None")
        self.pathogen_project_label.setStyleSheet("color: #2980B9; font-weight: bold;")
        self.pathogen_status_label = QLabel("● Idle")
        self.pathogen_status_label.setStyleSheet("color: #7F8C8D; font-weight: bold;")
        top_row.addWidget(self.pathogen_project_label)
        top_row.addStretch()
        top_row.addWidget(self.pathogen_status_label)
        layout.addLayout(top_row)
        
        # Input Settings
        input_group = QGroupBox("Input Settings")
        input_form = QFormLayout(input_group)
        
        # TaxID file dropdown (editable for Custom)
        self.taxid_file_combo = QComboBox()
        self.taxid_file_combo.setEditable(True)
        self.taxid_file_combo.addItems(list(PRESET_TAXID_FILES.keys()))
        self.taxid_file_combo.lineEdit().setPlaceholderText(
            "Enter custom TaxID file path.")
        self.taxid_file_combo.currentTextChanged.connect(self.on_taxid_file_changed)
        input_form.addRow("TaxID List File:", self.taxid_file_combo)
        
        # Virus Index dropdown (editable for Custom)
        self.virus_index_combo = QComboBox()
        self.virus_index_combo.setEditable(True)
        self.virus_index_combo.addItems(list(PRESET_VIRUS_INDEXES.keys()))
        self.virus_index_combo.lineEdit().setPlaceholderText(
            "Enter custom virus index path.")
        self.virus_index_combo.currentTextChanged.connect(self.on_virus_index_changed)
        input_form.addRow("Virus Index Path:", self.virus_index_combo)
        layout.addWidget(input_group)
        
        # Execution Mode
        exec_group = QGroupBox("Execution Mode")
        exec_form = QFormLayout(exec_group)
        self.pathogen_run_mode_combo = QComboBox()
        self.pathogen_run_mode_combo.addItems(["Direct Python", "Slurm (srun)"])
        exec_form.addRow("Run Mode:", self.pathogen_run_mode_combo)

        self.pathogen_cores_spin = QSpinBox()
        self.pathogen_cores_spin.setRange(1, 256)
        self.pathogen_cores_spin.setValue(
            int(self.pipeline_cfg.get('global', {}).get('threads', 64) or 64))
        self.pathogen_cores_spin.setToolTip(
            "CPU cores: passed to pathogen_detect --threads (shared by blastn / "
            "bowtie2 / samtools).\n"
            "Also used as srun --cpus-per-task in Slurm mode.")
        exec_form.addRow("CPU Cores:", self.pathogen_cores_spin)

        self.pathogen_slurm_extra_edit = QLineEdit()
        self.pathogen_slurm_extra_edit.setPlaceholderText(
            "Example: --exclusive")
        exec_form.addRow("Extra Slurm Args:", self.pathogen_slurm_extra_edit)
        self._slurm_rows.append(
            (exec_form.labelForField(self.pathogen_slurm_extra_edit),
             self.pathogen_slurm_extra_edit))
        layout.addWidget(exec_group)
        
        # Advanced Parameters Toggle
        self.pathogen_toggle_btn = QPushButton("▼ Advanced Parameters")
        self.pathogen_toggle_btn.setCheckable(True)
        layout.addWidget(self.pathogen_toggle_btn)
        
        # Advanced Parameters Area
        self.pathogen_adv_widget = QWidget()
        adv_vbox = QVBoxLayout(self.pathogen_adv_widget)
        adv_vbox.setContentsMargins(0, 0, 0, 0)

        self.pathogen_param_map = {}
        nr_cfg = self.pipeline_cfg.get('diamond', {}) or {}
        nr_params = nr_cfg.get('params', {}) or {}

        nr_group = QGroupBox("Setting: diamond (nr database)")
        nr_form = QFormLayout(nr_group)

        self.pathogen_diamond_db_edit = QLineEdit(str(nr_cfg.get('db_path', '') or ''))
        self.pathogen_diamond_db_edit.setPlaceholderText(
            _PARAM_HINTS[('diamond', 'db_path')])
        self.pathogen_diamond_db_edit.setToolTip(
            "Comprehensive nr protein database used by the nr DIAMOND step "
            "(runs after BLASTN, on the contigs already confirmed).")
        nr_form.addRow(f"{param_label(('diamond', 'db_path'))}:",
                       self.pathogen_diamond_db_edit)

        self.pathogen_diamond_bin_edit = QLineEdit(
            str(nr_cfg.get('bin', 'diamond') or 'diamond'))
        self.pathogen_diamond_bin_edit.setPlaceholderText(
            _PARAM_HINTS[('diamond', 'bin')])
        nr_form.addRow(f"{param_label(('diamond', 'bin'))}:",
                       self.pathogen_diamond_bin_edit)

        self.pathogen_nr_sensitivity_edit = QLineEdit(
            str(nr_params.get('sensitivity', '--sensitive')))
        self.pathogen_nr_evalue_edit = QLineEdit(str(nr_params.get('evalue', '1e-5')))
        self.pathogen_nr_max_target_edit = QLineEdit(
            str(nr_params.get('max_target_seqs', 10)))
        self.pathogen_nr_outfmt_edit = QLineEdit(str(nr_params.get('outfmt', '')))
        self.pathogen_nr_blocks_edit = QLineEdit(str(nr_params.get('blocks', 0)))
        self.pathogen_nr_gb_per_block_edit = QLineEdit(
            str(nr_params.get('gb_per_block', 6)))
        self.pathogen_nr_mem_reserve_edit = QLineEdit(
            str(nr_params.get('mem_reserve_gb', 30)))
        self.pathogen_nr_index_chunks_edit = QLineEdit(
            str(nr_params.get('index_chunks', 1)))
        nr_form.addRow(f"{param_label(('diamond', 'params', 'sensitivity'))}:",
                       self.pathogen_nr_sensitivity_edit)
        nr_form.addRow(f"{param_label(('diamond', 'params', 'evalue'))}:",
                       self.pathogen_nr_evalue_edit)
        nr_form.addRow(f"{param_label(('diamond', 'params', 'max_target_seqs'))}:",
                       self.pathogen_nr_max_target_edit)
        nr_form.addRow(f"{param_label(('diamond', 'params', 'outfmt'))}:",
                       self.pathogen_nr_outfmt_edit)
        nr_form.addRow(f"{param_label(('diamond', 'params', 'blocks'))}:",
                       self.pathogen_nr_blocks_edit)
        nr_form.addRow(f"{param_label(('diamond', 'params', 'gb_per_block'))}:",
                       self.pathogen_nr_gb_per_block_edit)
        nr_form.addRow(f"{param_label(('diamond', 'params', 'mem_reserve_gb'))}:",
                       self.pathogen_nr_mem_reserve_edit)
        nr_form.addRow(f"{param_label(('diamond', 'params', 'index_chunks'))}:",
                       self.pathogen_nr_index_chunks_edit)

        self.pathogen_node_mem_spin = QSpinBox()
        self.pathogen_node_mem_spin.setRange(0, 8192)
        self.pathogen_node_mem_spin.setSingleStep(16)
        self.pathogen_node_mem_spin.setSuffix(" GB")
        self.pathogen_node_mem_spin.setSpecialValueText("not set (DIAMOND defaults)")
        self.pathogen_node_mem_spin.setValue(int(nr_params.get('node_mem_gb', 0) or 0))
        self.pathogen_node_mem_spin.setToolTip(
            "Total memory of the compute node in GB. Used to derive the nr "
            "DIAMOND -b (memory ~ 6 GB x (b / c), c fixed to 1):\n"
            "231 GB node - 30 GB reserve = 201 GB -> b = 33.\n"
            "0 = not set: DIAMOND keeps its own defaults (-b 2 -c 4).\n"
            "This is NOT a Slurm memory limit; no --mem is requested.")
        nr_form.addRow(f"{param_label(('diamond', 'params', 'node_mem_gb'))}:",
                       self.pathogen_node_mem_spin)
        adv_vbox.addWidget(nr_group)

        for path_tuple, widget in (
            (('diamond', 'db_path'), self.pathogen_diamond_db_edit),
            (('diamond', 'bin'), self.pathogen_diamond_bin_edit),
            (('diamond', 'params', 'sensitivity'), self.pathogen_nr_sensitivity_edit),
            (('diamond', 'params', 'evalue'), self.pathogen_nr_evalue_edit),
            (('diamond', 'params', 'max_target_seqs'), self.pathogen_nr_max_target_edit),
            (('diamond', 'params', 'outfmt'), self.pathogen_nr_outfmt_edit),
            (('diamond', 'params', 'blocks'), self.pathogen_nr_blocks_edit),
            (('diamond', 'params', 'gb_per_block'), self.pathogen_nr_gb_per_block_edit),
            (('diamond', 'params', 'mem_reserve_gb'), self.pathogen_nr_mem_reserve_edit),
            (('diamond', 'params', 'index_chunks'), self.pathogen_nr_index_chunks_edit),
        ):
            self.pathogen_param_map[path_tuple] = widget

        # Filter Parameters (Diamond)
        filter_group = QGroupBox("Filter Parameters (Diamond)")
        filter_form = QFormLayout(filter_group)
        self.diamond_pident_edit = QLineEdit("50")
        self.diamond_evalue_edit = QLineEdit("1e-2")
        filter_form.addRow("Identity (%):", self.diamond_pident_edit)
        filter_form.addRow("E-value:", self.diamond_evalue_edit)
        adv_vbox.addWidget(filter_group)
        
        # BLASTN Parameters
        blastn_group = QGroupBox("BLASTN Parameters")
        blastn_form = QFormLayout(blastn_group)
        self.blastn_identity_edit = QLineEdit("70")
        self.blastn_evalue_edit = QLineEdit("1e-4")
        self.blastn_coverage_edit = QLineEdit("70")
        blastn_form.addRow("Identity (%):", self.blastn_identity_edit)
        blastn_form.addRow("E-value:", self.blastn_evalue_edit)
        blastn_form.addRow("Coverage (%):", self.blastn_coverage_edit)
        adv_vbox.addWidget(blastn_group)
        
        # Output Settings
        output_group = QGroupBox("Output Settings")
        output_form = QFormLayout(output_group)
        self.pathogen_output_folder_edit = QLineEdit("pathogen_results")
        output_form.addRow("Folder Name:", self.pathogen_output_folder_edit)
        adv_vbox.addWidget(output_group)
        
        #  bedtools.bin / blastn.bin / bedgraph_to_bigwig.bin）
        tools_hint = QLabel(
            "Tool paths are configured once in the Advanced panel "
            "(config2.yaml sections: bowtie2 / samtools / bedtools / "
            "blastn / bedgraph_to_bigwig; the nr DIAMOND bin sits above.)")
        tools_hint.setWordWrap(True)
        adv_vbox.addWidget(tools_hint)
        
        self.pathogen_adv_scroll = QScrollArea()
        self.pathogen_adv_scroll.setWidget(self.pathogen_adv_widget)
        self.pathogen_adv_scroll.setWidgetResizable(True)
        self.pathogen_adv_scroll.setMinimumHeight(200)
        self.pathogen_adv_scroll.setMaximumHeight(300)
        self.pathogen_adv_scroll.setVisible(False)
        layout.addWidget(self.pathogen_adv_scroll)
        self.pathogen_toggle_btn.clicked.connect(
            lambda: self.pathogen_adv_scroll.setVisible(
                self.pathogen_toggle_btn.isChecked())
        )
        
        settings_group = QGroupBox("Settings")
        settings_row = QHBoxLayout(settings_group)
        self.save_pathogen_cb = QCheckBox("Save Settings")
        self.save_pathogen_cb.setChecked(True)
        self.verbose_cb = QCheckBox("Verbose Logging")
        self.verbose_cb.setChecked(False)
        self.resume_cb = QCheckBox("Resume from checkpoint (if exists)")
        self.resume_cb.setChecked(True)
        settings_row.addWidget(self.save_pathogen_cb)
        settings_row.addWidget(self.verbose_cb)
        settings_row.addWidget(self.resume_cb)
        self.pathogen_sync_db_btn = QPushButton("⇪ Sync DB to Server")
        self.pathogen_sync_db_btn.setToolTip(
            "Write the database paths (diamond db / bowtie2 index / BLASTN "
            "virus index / taxid list) from Settings into the server-side "
            "configs/config2.yaml.")
        self.pathogen_sync_db_btn.clicked.connect(
            lambda: self.sync_db_to_server(self.pathogen_log_text))
        settings_row.addWidget(self.pathogen_sync_db_btn)
        settings_row.addStretch()
        layout.addWidget(settings_group)
        
        # Run Button
        self.pathogen_run_btn = QPushButton("🚀 RUN DETECTION")
        self.pathogen_run_btn.clicked.connect(self.on_pathogen_run)
        layout.addWidget(self.pathogen_run_btn)

        self.stop_pathogen_btn = QPushButton("⏸ STOP DETECTION")
        self.stop_pathogen_btn.setToolTip("Are you sure to terminate the running task?")
        self.stop_pathogen_btn.clicked.connect(self.on_stop_pathogen)
        layout.addWidget(self.stop_pathogen_btn)
        
        self.pathogen_log_text = QTextEdit()
        self.pathogen_log_text.setReadOnly(True)
        self.pathogen_log_text.setPlaceholderText(
            "Not started. Run the pre-processing pipeline first, then click "
            "RUN DETECTION.")
        self.pathogen_log_text.setStyleSheet(
            "background-color: #FFFFFF; color: #000000; font-family: 'Consolas';")
        self.pathogen_log_text.setMinimumHeight(140)
        layout.addWidget(self.pathogen_log_text, stretch=1)
        
        result_btn_layout = QHBoxLayout()
        self.view_results_btn = QPushButton("💾 View Results")
        self.view_results_btn.clicked.connect(self.on_view_pathogen_results)
        
        self.export_csv_btn = QPushButton("💾 Export CSV")
        self.export_csv_btn.clicked.connect(self.on_export_pathogen_csv)
        
        result_btn_layout.addWidget(self.view_results_btn)
        result_btn_layout.addWidget(self.export_csv_btn)
        result_btn_layout.addStretch()
        layout.addLayout(result_btn_layout)

        self.tabs.addTab(content, "3. Pathogen Detection")

    NEW_PROJECT_LABEL = "➕ New Project…"

    def on_list_remote(self):
        """
        List one level below the active project (`ls -lah`).
        """
        path = (self.current_project_dir or "").strip()
        if not path:
            QMessageBox.warning(
                self, "No Active Project",
                "Please click “📁 Select / New Project…” first.\n"
                "“🔍 List Files” only lists files inside the Active Project.")
            return
        quoted = shlex.quote(path)
        command = f"ls -lah {quoted}"
        self.log_text.append(f"📂 Listing files in: {path}")
        self.run_ssh_task(command=command, log_widget=self.log_text)

    def _ssh_capture(self, command, timeout=20):
        """
        Run a remote command synchronously and return (stdout, stderr).
        """
        pk = _paramiko()
        ssh = pk.SSHClient()
        ssh.set_missing_host_key_policy(pk.AutoAddPolicy())
        ssh.connect(self.host_edit.text().strip(), port=int(self.port_edit.text()),
                    username=self.user_edit.text().strip(),
                    password=self.pwd_edit.text(), timeout=timeout)
        try:
            _, stdout, stderr = ssh.exec_command(command, timeout=timeout)
            return (stdout.read().decode("utf-8", "replace"),
                    stderr.read().decode("utf-8", "replace"))
        finally:
            ssh.close()

    def _set_active_project(self, path):
        changed = (path or "") != (self.current_project_dir or "")
        old_project = self.current_project_dir or ""
        self.current_project_dir = path
        self.sync_info_label.setText(f"Active Project: {path}")
        if hasattr(self, "list_dir_btn"):
            self.list_dir_btn.setEnabled(bool(path))
        if hasattr(self, "pipeline_project_label"):
            self.pipeline_project_label.setText(f"📁 Project: {path}")
        if hasattr(self, "pathogen_project_label"):
            self.pathogen_project_label.setText(f"📁 Project: {path}")
        if changed and hasattr(self, "pipe_log_text"):
            self._reset_tab_outputs(old_project)
        self.schedule_autosave()
        if changed and path:
            self.refresh_tab_statuses()

    def _reset_tab_outputs(self, old_project=""):
        """
        Clear log boxes and status when the project changes, so that leftover
        messages from the previous project cannot bleed into the new one.
        """
        if old_project:
            prefix = old_project.rstrip("/") + "/"
            for target in list(self._poll_timers):
                if target == old_project or target.startswith(prefix):
                    self._stop_log_poll(target)
        for box in (getattr(self, "pipe_log_text", None),
                    getattr(self, "pathogen_log_text", None)):
            if box is None:
                continue
            box.clear()
            self._task_rc.pop(id(box), None)
            self._task_started.pop(id(box), None)
        self._once_keys = set()
        self.set_tab_status(1)

    def on_select_project(self):
        """
        List first-level directories under the project root; the first entry creates
        a new project.
        """
        base = self.search_dir_edit.text().strip().rstrip("/")
        if not base:
            QMessageBox.warning(self, "Path Error", "Please fill in the project root first.")
            return
        if not self.host_edit.text().strip():
            QMessageBox.warning(self, "SSH Error", "Please fill in the server settings first.")
            return

        try:
            out, _ = self._ssh_capture(
                f"find {shlex.quote(base)} -mindepth 1 -maxdepth 1 -type d -printf '%f\\n' "
                f"2>/dev/null | sort")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "SSH Error", f"Cannot list projects under {base}:\n{exc}")
            return

        projects = [line.strip() for line in out.splitlines() if line.strip()]
        items = [self.NEW_PROJECT_LABEL] + projects
        choice, ok = QInputDialog.getItem(
            self, "Select Project",
            f"Project root:\n{base}\n\nSelect an existing project or create a new one:",
            items, 0, False)
        if not ok:
            return

        if choice == self.NEW_PROJECT_LABEL:
            name, ok2 = QInputDialog.getText(self, "New Project", "Project folder name:")
            name = (name or "").strip().strip("/")
            if not (ok2 and name):
                return
            if "/" in name or name in (".", ".."):
                QMessageBox.warning(self, "Invalid Name",
                                    "Folder name must not contain '/' / be '.' or '..'")
                return
            self.run_ssh_task(target_dir=f"{base}/{name}", task_type="mkdir",
                              log_widget=self.log_text)
        else:
            self._set_active_project(f"{base}/{choice}")
            self.log_text.append(f"📁 Using existing project: {self.current_project_dir}")

    def on_task_finished(self, success, msg):
        task_type = getattr(self.worker, 'task_type_ext', 'default')
        target_dir = getattr(self.worker, 'target_dir', None) or ""

        if task_type == "pathogen":
            if target_dir and not self._once(("finished", target_dir)):
                return
            self._report_task_end(
                self.pathogen_log_text, target_dir, "Pathogen detection",
                success, msg, self.pathogen_output_folder)
            return

        if task_type == "preprocess":
            if target_dir and not self._once(("finished", target_dir)):
                return
            step_dir = target_dir or f"{self.current_project_dir}/preprocess"
            self._report_task_end(
                self.pipe_log_text, target_dir, "Pre-processing pipeline",
                success, msg, step_dir,
                log_path=f"{step_dir}/{SSHWorker.RUN_LOG_NAME}")
            return
        
        if success and msg == "CREATED":
            self._set_active_project(self.worker.target_dir)
            QMessageBox.information(self, "Success", f"Workspace ready:\n{self.current_project_dir}")
        elif msg == "EXISTS":
            QMessageBox.warning(self, "Conflict", "Directory already exists!")

    def _elapsed_text(self, log_widget, target_dir):
        """
        Return `elapsed 12m 05s`, or an empty string when the start time is unknown.
        """
        started = self._task_started.pop(id(log_widget), None)
        if started is None:
            started = (self._poll_state.get(target_dir) or {}).get("started")
        if not started:
            return ""
        secs = max(0, int(time.monotonic() - started))
        if secs < 60:
            return f"elapsed {secs}s"
        if secs < 3600:
            return f"elapsed {secs // 60}m {secs % 60:02d}s"
        return f"elapsed {secs // 3600}h {(secs % 3600) // 60:02d}m"

    def _report_task_end(self, log_widget, target_dir, label, success, msg,
                         out_dir=None, log_path=None):
        """
        Single end-of-task report: success/failure decided by the remote exit code,
        with the elapsed time attached.
        """
        rc = self._task_rc.pop(id(log_widget), None)
        ok = bool(success) and rc in (None, 0)
        elapsed = self._elapsed_text(log_widget, target_dir)
        log_widget.append("=" * 50)
        if ok:
            log_widget.append(f"✅ {label} finished!"
                              + (f"   ({elapsed})" if elapsed else ""))
        else:
            if rc is not None:
                reason = ("exit code unknown" if rc < 0
                          else f"exit code {rc}")
                if msg and msg != "SUCCESS":
                    reason += f"; {msg}"
            else:
                reason = msg or "unknown error"
            log_widget.append(f"❌ {label} failed ({reason})"
                              + (f"   ({elapsed})" if elapsed else ""))
            log_widget.append("   Check the log above for the first failing step.")
        if out_dir:
            log_widget.append(f"Output dir: {out_dir}")
        if log_path:
            log_widget.append(f"Full log : {log_path}")
        log_widget.append("=" * 50)
        self._set_tab_status(log_widget, 3 if ok else 1)

    @staticmethod
    def _user_sets_cpus(extra):
        """
        Whether the user pinned cpus-per-task in the extra Slurm arguments.
        """
        return any(tok == "-c" or tok.startswith("--cpus-per-task")
                   for tok in (extra or "").split())

    def _slurm_script(self, base_cmd, cores=None, job_name="",
                      target_dir=""):
        """
        Build the sbatch script.

        Resources are expressed as `#SBATCH` directives so the script is
        self-describing: `sbatch .hvfinder.job.sh` by hand gives the same resources,
        and the GUI only adds `--parsable` plus the user's extra arguments. Memory is
        deliberately not limited - the requirement depends on the input, and a fixed
        ceiling would kill legitimate jobs; set `--mem=` in the extra arguments if a
        limit is really wanted.
        """
        lines = ["#!/bin/bash",
                 "# Generated by hvFinder. Change parameters in the GUI, not here."]
        if job_name:
            lines.append(f"#SBATCH --job-name={job_name}")
        if cores:
            lines.append(f"#SBATCH --cpus-per-task={cores}")
        if target_dir:
            lines.append(f"#SBATCH --chdir={shlex.quote(target_dir)}")
            out = shlex.quote(f"{target_dir}/{SSHWorker.SLURM_OUT_NAME}")
            lines.append(f"#SBATCH --output={out}")
            lines.append(f"#SBATCH --error={out}")
        lines.append("set -uo pipefail")
        lines.append(f"bash -l -c {shlex.quote(base_cmd)}")
        lines.append("rc=$?")
        lines.append('echo "[hvfinder] exit=$rc"')
        lines.append("exit $rc")
        return "\n".join(lines) + "\n"

    def _compose_command(self, base_cmd, mode_text, extra, job_name,
                         jobid_file=None, target_dir=None):
        """
        Compose the final command for the selected run mode.

        Slurm: `sbatch --parsable [extra] <target_dir>/.hvfinder.job.sh`; cores and
        output live in the script. Direct: run it right away, no scheduler, no
        ulimit.
        """
        if mode_text.startswith("Slurm"):
            extra = (extra or "").strip()
            opts = ["--parsable"]
            if extra:
                opts.append(extra)
            if target_dir:
                out = shlex.quote(f"{target_dir}/{SSHWorker.SLURM_OUT_NAME}")
                script = shlex.quote(
                    f"{target_dir}/{SSHWorker.JOB_SCRIPT_NAME}")
                clean = f"rm -f {out}"
                if jobid_file:
                    clean += f" {shlex.quote(jobid_file)}"
                clean += "; "
            else:
                script = SSHWorker.JOB_SCRIPT_NAME
                clean = ""
            cmd = f"{clean}sbatch {' '.join(opts)} {script}"
            if jobid_file:
                cmd += f" | head -1 | tee {shlex.quote(jobid_file)}"
            return cmd
        return base_cmd

    def on_launch(self, step_name="preprocess"):
        if not self.current_project_dir:
            QMessageBox.warning(self, "Path Error", "No project folder selected!")
            return

        r1 = self.fq1_edit.text().strip()
        r2 = self.fq2_edit.text().strip()
        if not r1 or not r2:
            QMessageBox.warning(self, "Input Error",
                                "Please select both R1 and R2 FASTQ files.")
            return
        if r1 == r2:
            QMessageBox.warning(self, "Input Error",
                                "R1 and R2 must be two different files.")
            return

        self._guard_slurm_submission(
            f"mNGS_{step_name}", self.pipe_log_text,
            lambda: self._launch_pipeline(step_name))

    def _launch_pipeline(self, step_name="preprocess"):
        """
        Submit the pre-processing run (duplicate Slurm jobs already ruled out).
        """
        target_sub_dir = f"{self.current_project_dir}/{step_name}"
        self._reset_once(("finished", target_sub_dir), ("status-error", target_sub_dir))
        sync_yaml = self.get_ui_as_yaml()
        
        if step_name == "preprocess":
            extra_args = ""
            if hasattr(self, 'pipeline_verbose_cb') and self.pipeline_verbose_cb.isChecked():
                extra_args += " --verbose"
            if hasattr(self, 'pipeline_resume_cb') and not self.pipeline_resume_cb.isChecked():
                extra_args += " --no-resume"
            
            cores = str(self.cores_spin.value())
            base_cmd = (f"mNGS_pip -1 {shlex.quote(self.fq1_edit.text().strip())} "
                        f"-2 {shlex.quote(self.fq2_edit.text().strip())} "
                        f"-o results -t {cores}{extra_args}")
        else:
            QMessageBox.warning(self, "Unsupported", f"Unknown step: {step_name}")
            return
        
        mode = self.run_mode_combo.currentText()
        is_slurm = mode.startswith("Slurm")
        slurm_extra = self.slurm_extra_edit.text().strip()
        upload_files = {}
        if is_slurm:
            upload_files[SSHWorker.JOB_SCRIPT_NAME] = self._slurm_script(
                base_cmd,
                cores=None if self._user_sets_cpus(slurm_extra) else cores,
                job_name=f"mNGS_{step_name}", target_dir=target_sub_dir)
        cmd = self._compose_command(base_cmd, mode, slurm_extra,
                                    f"mNGS_{step_name}",
                                    jobid_file=f"{target_sub_dir}/{SSHWorker.RUN_JOB_NAME}",
                                    target_dir=target_sub_dir)

        self.pipe_log_text.clear()
        self._task_started[id(self.pipe_log_text)] = time.monotonic()
        self._task_rc.pop(id(self.pipe_log_text), None)
        self.pipe_log_text.append("Starting pre-processing pipeline...")
        self.pipe_log_text.append(f"Project: {self.current_project_dir}")
        self.pipe_log_text.append(f"Output: {target_sub_dir}")
        self.pipe_log_text.append(f"CPU cores: {cores} (no --mem requested)")
        self.pipe_log_text.append(f"Run mode: {mode}")
        self.pipe_log_text.append("-" * 50)
        self.sync_db_to_server(self.pipe_log_text)
        self.run_ssh_task(command=cmd, upload_cfg=sync_yaml, target_dir=target_sub_dir,
                          log_widget=self.pipe_log_text,
                          task_type="submit" if is_slurm else "preprocess",
                          detached=not is_slurm, upload_files=upload_files)
        if is_slurm:
            self.pipe_log_text.append(
                f"📄 log: {target_sub_dir}/{self._log_glob('preprocess')}")
            self._poll_state.setdefault(target_sub_dir, {})["expect_rc"] = True
            self._start_log_poll(target_sub_dir, self.pipe_log_text,
                                 task_type="preprocess", reset=True)
        self._announce_job_id(target_sub_dir, self.pipe_log_text)

    def adv_value(self, section, key, default=""):
        """
        Read one value from the Advanced panel; `default` when the key is absent.
        """
        edit = self.ui_param_map.get((section, key))
        value = edit.text().strip() if edit is not None else ""
        return value or default

    def on_pathogen_run(self):
        if not self.current_project_dir:
            QMessageBox.warning(self, "Path Error", "No project folder selected!")
            return

        self.ensure_script_base()

        if not self.get_taxid_file_path():
            QMessageBox.warning(self, "Input Error",
                                "Please select or enter TaxID List File path!")
            return
        if not self.get_virus_index_path():
            QMessageBox.warning(self, "Input Error",
                                "Please select or enter Virus Index Path!")
            return

        self._check_candidates_then_run()

    def _check_candidates_then_run(self):
        """
        Confirm asynchronously that `preprocess/prescreen_candidates.fasta` exists
        and is not empty.
        """
        path = (f"{self.current_project_dir}/preprocess/"
                f"prescreen_candidates.fasta")
        if not self.host_edit.text().strip():
            self._guard_slurm_submission("pathogen", self.pathogen_log_text,
                                        self._run_pathogen_now)
            return
        cmd = f"if [ -s {shlex.quote(path)} ]; then echo OK; else echo MISSING; fi"
        self._start_remote_query(cmd, self._on_candidates_check,
                                 context=(path,), attr="_cand_query")

    def _on_candidates_check(self, ok, out, context, query):
        if not self._is_current_query(query, getattr(self, "_cand_query", None)):
            return
        path = context[0]
        lines = [ln.strip() for ln in (out or "").splitlines() if ln.strip()]
        verdict = lines[-1] if lines else ""
        if ok and verdict == "OK":
            self._guard_slurm_submission("pathogen", self.pathogen_log_text,
                                        self._run_pathogen_now)
            return
        if not ok:
            self.pathogen_log_text.append(
                f"⚠️ Cannot verify {path} ({(out or '').strip()}); "
                f"submitting anyway.")
            self._guard_slurm_submission("pathogen", self.pathogen_log_text,
                                        self._run_pathogen_now)
            return
        QMessageBox.warning(
            self, "Missing input",
            "Pre-processing output not found or empty:\n\n"
            f"{path}\n\n"
            "Run the Pre-processing pipeline (Tab 2) first, and make sure "
            "PRESCREEN kept at least one candidate contig.")
        self.pathogen_log_text.append(
            f"⏹ Submission cancelled: missing or empty {path}")

    def _run_pathogen_now(self):
        """
        Submit the detection run (inputs validated, duplicate jobs ruled out).
        """
        taxid_file = self.get_taxid_file_path()
        virus_index = self.get_virus_index_path()
        
        # Get output folder
        self.pathogen_output_folder = f"{self.current_project_dir}/{self.pathogen_output_folder_edit.text().strip()}"
        self._reset_once(("finished", self.pathogen_output_folder),
                         ("status-error", self.pathogen_output_folder))
        
        # Build command
        pident = self.diamond_pident_edit.text().strip()
        evalue = self.diamond_evalue_edit.text().strip()
        threads = str(self.pathogen_cores_spin.value())
        blastn_id = self.blastn_identity_edit.text().strip()
        blastn_evalue = self.blastn_evalue_edit.text().strip()
        blastn_cov = self.blastn_coverage_edit.text().strip()
        
        tab2_step = "preprocess"
        candidates_fa = f"{self.current_project_dir}/{tab2_step}/prescreen_candidates.fasta"
        megahit_results = f"{self.current_project_dir}/{tab2_step}/results_megahit_result"
        diamond_db = self.pathogen_diamond_db_edit.text().strip()
        diamond_bin = self.pathogen_diamond_bin_edit.text().strip() or 'diamond'
        node_mem_gb = int(self.pathogen_node_mem_spin.value())
        
        fq1 = self.fq1_edit.text()
        fq2 = self.fq2_edit.text()
        
        extra_args = ""
        if self.verbose_cb.isChecked():
            extra_args += " --verbose"
        if not self.resume_cb.isChecked():
            extra_args += " --no-resume"
        
        bowtie2_bin = self.adv_value('bowtie2', 'bin', 'bowtie2')
        samtools_bin = self.adv_value('samtools', 'bin', 'samtools')
        bedtools_bin = self.adv_value('bedtools', 'bin', 'bedtools')
        blastn_bin = self.adv_value('blastn', 'bin', 'blastn')
        bigwig_bin = self.adv_value('bedgraph_to_bigwig', 'bin', 'bedGraphToBigWig')

        base_cmd = f"pathogen_detect \
            --candidates {candidates_fa} \
            --diamond-db {diamond_db} \
            --diamond-bin {diamond_bin} \
            --node-mem-gb {node_mem_gb} \
            --fq1 {fq1} --fq2 {fq2} \
            --megahit {megahit_results} \
            --taxids {taxid_file} \
            --virus-index {virus_index} \
            --output {self.pathogen_output_folder} \
            --bowtie2-bin {bowtie2_bin} \
            --samtools-bin {samtools_bin} \
            --bedtools-bin {bedtools_bin} \
            --bedgraph-to-bigwig-bin {bigwig_bin} \
            --blastn-bin {blastn_bin} \
            --pident {pident} --evalue {evalue} \
            --blastn-identity {blastn_id} --blastn-evalue {blastn_evalue} --blastn-coverage {blastn_cov} \
            --threads {threads}{extra_args}"
        
        # Handle run mode
        mode = self.pathogen_run_mode_combo.currentText()
        is_slurm = mode.startswith("Slurm")
        slurm_extra = self.pathogen_slurm_extra_edit.text().strip()
        upload_files = {}
        if is_slurm:
            upload_files[SSHWorker.JOB_SCRIPT_NAME] = self._slurm_script(
                base_cmd,
                cores=None if self._user_sets_cpus(slurm_extra) else threads,
                job_name="pathogen",
                target_dir=self.pathogen_output_folder)
        cmd = self._compose_command(
            base_cmd, mode, slurm_extra, "pathogen",
            jobid_file=f"{self.pathogen_output_folder}/{SSHWorker.RUN_JOB_NAME}",
            target_dir=self.pathogen_output_folder)
        
        # Run
        self.pathogen_log_text.clear()
        self._task_started[id(self.pathogen_log_text)] = time.monotonic()
        self._task_rc.pop(id(self.pathogen_log_text), None)
        self.pathogen_log_text.append(f"Starting pathogen detection...")
        self.pathogen_log_text.append(f"Project: {self.current_project_dir}")
        self.pathogen_log_text.append(f"Output: {self.pathogen_output_folder}")
        self.pathogen_log_text.append(f"CPU cores: {threads} (no --mem requested)")
        self.pathogen_log_text.append(f"Run mode: {mode}")
        self.pathogen_log_text.append("-" * 50)
        
        self.sync_db_to_server(self.pathogen_log_text)
        self.run_ssh_task(
            command=cmd,
            target_dir=self.pathogen_output_folder,
            log_widget=self.pathogen_log_text,
            task_type="submit" if is_slurm else "pathogen",
            finished_callback=self.on_task_finished,
            detached=not is_slurm,
            upload_files=upload_files,
        )
        if is_slurm:
            self.pathogen_log_text.append(
                f"📄 log: {self.pathogen_output_folder}/{self._log_glob('pathogen')}")
            self._poll_state.setdefault(
                self.pathogen_output_folder, {})["expect_rc"] = True
            self._start_log_poll(self.pathogen_output_folder,
                                 self.pathogen_log_text,
                                 task_type="pathogen", reset=True)
        self._announce_job_id(self.pathogen_output_folder, self.pathogen_log_text)

    def on_view_pathogen_results(self):
        if not self.current_project_dir:
            QMessageBox.warning(self, "Path Error", "No project folder selected!")
            return
        
        results_file = f"{self.current_project_dir}/{self.pathogen_output_folder_edit.text().strip()}/detection_results.tsv"
        
        # Try to read results from remote
        try:
            pk = _paramiko()
            ssh = pk.SSHClient()
            ssh.set_missing_host_key_policy(pk.AutoAddPolicy())
            ssh.connect(self.host_edit.text().strip(), port=int(self.port_edit.text()),
                       username=self.user_edit.text().strip(), password=self.pwd_edit.text(), timeout=10)
            
            stdin, stdout, stderr = ssh.exec_command(f"cat {results_file}")
            content = stdout.read().decode('utf-8')
            ssh.close()
            
            old_header_v1 = "Virus_Name\tTaxID\tCoverage\tDepth\tReads_Count\tContigs_Count\tContig_IDs"
            old_header_v2 = "Virus_Name\tTaxID\tBest_Hit_Seq\tCoverage(%)\tDepth\tContig_ID"
            new_header = "Virus_Seq\tBest_Hit_Seq\tCoverage(%)\tDepth\tContig_ID"
            if not content or content.strip() == old_header_v1 or content.strip() == old_header_v2 or content.strip() == new_header:
                QMessageBox.information(self, "No Results", "No pathogen detection results found.")
                return

            results, has_taxid = parse_detection_results(content)
            if not results:
                QMessageBox.information(self, "No Results", "No pathogen detected.")
                return
            dialog = ResultsDialog(results, self.current_project_dir, has_taxid, self)
            dialog.exec_()

        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to read results:\n{str(e)}")

    def get_remote_script_base(self, timeout=15):
        """
        Resolve the remote hvFinder directory synchronously (blocking).

        Only used when the cached value is still missing at Run time.
        """
        try:
            out, _ = self._ssh_capture("bash -l -c 'get_script_base' 2>&1",
                                       timeout=timeout)
        except Exception:  # noqa: BLE001
            return ""
        self._script_base = parse_script_base(out)
        return self._script_base

    def request_script_base(self):
        """
        Resolve the remote hvFinder directory in a background thread (cached).

        Doing this synchronously from the combo-box handlers used to add 6-12 s to
        GUI start-up.
        """
        if self._script_base or not self.host_edit.text().strip():
            return
        running = getattr(self, "_script_base_query", None)
        if self._thread_running(running):
            return
        self._start_remote_query("bash -l -c 'get_script_base' 2>&1",
                                 self._on_script_base, attr="_script_base_query")

    def _on_script_base(self, ok, out, context, query):
        if not self._is_current_query(query, getattr(self, "_script_base_query", None)):
            return
        if ok:
            self._script_base = parse_script_base(out)

    def ensure_script_base(self):
        """
        Make sure the remote script root is known (cache first, then one blocking
        call - only when Run is pressed).
        """
        if self._script_base:
            return self._script_base
        return self.get_remote_script_base()

    def on_virus_index_changed(self, text):
        """
        Handle a change of the Virus Index combo box.
        """
        self.request_script_base()

        if text == "Custom":
            self.virus_index_combo.lineEdit().clear()

    def on_taxid_file_changed(self, text):
        """
        Handle a change of the TaxID File combo box.
        """
        self.request_script_base()

        if text == "Custom":
            self.taxid_file_combo.lineEdit().clear()

    def resolve_preset_path(self, relative_path):
        """
        Turn a preset relative path into an absolute one (falls back to the relative
        path while the remote root is still unknown).
        """
        if not relative_path:
            return ""
        if not self._script_base:
            self.request_script_base()
            return relative_path
        return f"{self._script_base}/{relative_path}"

    def get_virus_index_path(self):
        """
        Currently selected Virus Index path.
        """
        current_text = self.virus_index_combo.currentText()
        if current_text in PRESET_VIRUS_INDEXES and current_text != "Custom":
            return self.resolve_preset_path(PRESET_VIRUS_INDEXES[current_text])
        else:
            return current_text

    def get_taxid_file_path(self):
        """
        Currently selected TaxID File path.
        """
        current_text = self.taxid_file_combo.currentText()
        if current_text in PRESET_TAXID_FILES and current_text != "Custom":
            return self.resolve_preset_path(PRESET_TAXID_FILES[current_text])
        else:
            return current_text

    def on_export_pathogen_csv(self):
        self.on_view_pathogen_results()

    def _worker_busy_for(self, log_widget):
        """
        Whether a worker is already streaming into this log box.
        """
        worker = getattr(self, "worker", None)
        return (self._thread_running(worker)
                and getattr(worker, "log_widget", None) is log_widget)

    def _once(self, key):
        """
        Let a key through exactly once (True on the first call).
        """
        seen = getattr(self, "_once_keys", None)
        if seen is None:
            seen = self._once_keys = set()
        if key in seen:
            return False
        seen.add(key)
        return True

    def _reset_once(self, *keys):
        """
        Clear the one-shot markers (new task, or a task was seen running again).
        """
        seen = getattr(self, "_once_keys", None)
        if not seen:
            return
        for key in keys:
            seen.discard(key)

    def _start_remote_query(self, command, handler, context=None,
                            attr="_reattach_query", timeout=10):
        """
        Run a remote command on a background thread and call handler(ok, out, ctx).
        """
        host = self.host_edit.text().strip()
        user = self.user_edit.text().strip()
        if not host or not user:
            return None
        try:
            port = int(self.port_edit.text().strip() or 22)
        except ValueError:
            port = 22
        query = RemoteQuery(host, port, user, self.pwd_edit.text(), command,
                            timeout=timeout)
        query.context = context
        query.result_signal.connect(
            lambda ok, out, q=query: handler(ok, out, q.context, q))
        query.finished.connect(query.deleteLater)
        setattr(self, attr, query)
        query.start()
        return query

    @staticmethod
    def _is_current_query(query, attr_value):
        return query is not None and query is attr_value

    @staticmethod
    def _thread_running(worker):
        """
        Whether a QThread is still alive.

        Background queries call deleteLater(), so the C++ object may already be
        gone; calling into the Python wrapper then raises RuntimeError. Such a
        thread simply counts as not running.
        """
        if worker is None:
            return False
        try:
            return bool(worker.isRunning())
        except RuntimeError:
            return False
        except Exception:
            return False

    def _status_probe_command(self, target_dir, task_type="preprocess",
                              fallback_dir=None):
        """
        Build the remote status probe.

        Priority: (1) `.hvfinder.slurm.jobid` -> squeue while queued, sacct once it
        leaves the queue; (2) pid file + kill -0 (Direct mode); (3) NO_PID / NONE;
        (4) no GUI record but pipeline logs present (batch or CLI runs) ->
        PIPE_DONE / PIPE_PARTIAL / PIPE_FAILED, so a finished project is not shown as
        `Not started`.

        One line is printed: SLURM_RUNNING <jobid> <STATE> <pid|->,
        SLURM_DONE <jobid> <STATE|GONE>, PY_RUNNING/PY_DEAD <pid>, NO_PID, NONE,
        PIPE_DONE/PIPE_PARTIAL/PIPE_FAILED <log>.
        """
        pid_q = shlex.quote(f"{target_dir}/{SSHWorker.RUN_PID_NAME}")
        job_q = shlex.quote(f"{target_dir}/{SSHWorker.RUN_JOB_NAME}")
        log_q = shlex.quote(f"{target_dir}/{SSHWorker.RUN_LOG_NAME}")
        glob = self._log_glob(task_type)
        if task_type == "pathogen":
            done_test = '[ -f "$D/detection_results.tsv" ]'
        else:
            done_test = ('tail -n 60 "$LOG" | grep -q '
                         '"Pipeline finished successfully"')
        fallback = ""
        if fallback_dir:
            fallback = (f"D2={shlex.quote(fallback_dir)}; "
                        'if [ ! -d "$D" ] && [ -d "$D2" ]; then D="$D2"; fi; ')
        script = (
            f"D={shlex.quote(target_dir)}; " + fallback +
            f"P={pid_q}; J={job_q}; L={log_q}; "
            'PIDV=""; if [ -f "$P" ]; then Q=$(cat "$P" 2>/dev/null); '
            'if [ -n "$Q" ] && kill -0 "$Q" 2>/dev/null; then PIDV="$Q"; fi; fi; '
            'if [ -f "$J" ] && [ -n "$(cat "$J" 2>/dev/null)" ]; then '
            'JID=$(cat "$J"); JID=${JID%%;*}; '
            'ST=""; if command -v squeue >/dev/null 2>&1; then '
            'ST=$(squeue -h -j "$JID" -o "%T" 2>/dev/null | head -1); fi; '
            'if [ -n "$ST" ]; then '
            'echo "SLURM_RUNNING $JID $ST ${PIDV:--}"; '
            'elif ! command -v squeue >/dev/null 2>&1; then '
            'echo "SLURM_UNKNOWN $JID"; '
            "else "
            'FIN=$(sacct -n -X -j "$JID" -o State 2>/dev/null | head -1 '
            '| tr -d " "); '
            'echo "SLURM_DONE $JID ${FIN:-GONE}"; '
            "fi; "
            'elif [ -n "$PIDV" ]; then echo "PY_RUNNING $PIDV"; '
            'elif [ -f "$P" ]; then echo "PY_DEAD $(cat "$P" 2>/dev/null)"; '
            'elif [ -f "$L" ]; then echo "NO_PID"; '
            "else "
            f'LOG=$(ls -t "$D"/{glob} 2>/dev/null | head -1); '
            'if [ -n "$LOG" ]; then '
            'if tail -n 60 "$LOG" | grep -q "Traceback (most recent call last)" '
            '|| tail -n 60 "$LOG" | grep -q "Pipeline halted at Step"; then '
            'echo "PIPE_FAILED $LOG"; '
            f'elif {done_test}; then echo "PIPE_DONE $LOG"; '
            'else echo "PIPE_PARTIAL $LOG"; fi; '
            'else echo "NONE"; fi; '
            "fi"
        )
        return f"bash -l -c {shlex.quote(script)}"

    def reattach_running_task(self, target_dir, log_widget, label, task_type="cmd",
                              fallback_dir=None):
        """
        Re-attach to a still-running remote task (after a reconnect or GUI restart).
        Slurm is followed by job id, Direct mode by pid. The probe runs on a
        background thread.
        """
        if not target_dir or not self.host_edit.text().strip():
            return
        if self._worker_busy_for(log_widget):
            return
        if self._poll_active(target_dir):
            return

        log_file = f"{target_dir}/{SSHWorker.RUN_LOG_NAME}"
        log_q = shlex.quote(log_file)
        self._start_remote_query(
            self._status_probe_command(target_dir, task_type, fallback_dir),
            self._on_reattach_state,
            context=(target_dir, log_widget, label, task_type, log_file, log_q))

    def _on_reattach_state(self, ok, out, context, query):
        if not self._is_current_query(query, getattr(self, "_reattach_query", None)):
            return
        target_dir, log_widget, label, task_type, log_file, log_q = context
        if not ok:
            if self._once(("status-error", target_dir)):
                log_widget.append(f"⚠️ Cannot query remote task status: {out}")
            return

        state = out.strip().split()
        head = state[0] if state else "NONE"
        if head == "SLURM_UNKNOWN":
            job_id = state[1] if len(state) > 1 else "?"
            self._slurm_jobs[target_dir] = job_id
            if self._once(("state-unknown", target_dir, "reattach")):
                log_widget.append(
                    f"⚠️ Cannot query Slurm state for job {job_id} "
                    f"(squeue/sacct unavailable); following the saved log instead.")
            self._start_log_poll(target_dir, log_widget, task_type=task_type,
                                 reset=True)
            return
        if head == "SLURM_RUNNING":
            job_id = state[1] if len(state) > 1 else "?"
            slurm_state = state[2] if len(state) > 2 else "?"
            self._slurm_jobs[target_dir] = job_id
            self._reset_once(("finished", target_dir))
            self._set_tab_status(log_widget, 2)
            log_widget.append(
                f"🔄 Re-attached to running {label} "
                f"(slurm job {job_id}, state {slurm_state})")
            log_widget.append("   following the saved log (polling)")
            self._start_log_poll(target_dir, log_widget, task_type=task_type,
                                 reset=True)
            return
        elif head in ("PY_RUNNING", "RUNNING"):
            pid = state[1] if len(state) > 1 else "?"
            self._reset_once(("finished", target_dir))
            if self._worker_busy_for(log_widget):
                return
            self._set_tab_status(log_widget, 2)
            log_widget.append(f"🔄 Re-attached to running {label} (pid {pid})")
            log_widget.append(f"   log: {log_file}")
            self.run_ssh_task(command=f"tail -n 200 -s 1 --pid={pid} -f {log_q}",
                              target_dir=target_dir,
                              log_widget=log_widget, task_type=task_type,
                              finished_callback=self.on_task_finished)
        elif head in ("SLURM_DONE", "PY_DEAD", "DEAD", "NO_PID"):
            self._set_tab_status(log_widget, 3)
            if not self._once(("finished", target_dir)):
                return
            if head == "SLURM_DONE":
                job_id = state[1] if len(state) > 1 else "?"
                job_state = state[2] if len(state) > 2 else "GONE"
                self._slurm_jobs[target_dir] = job_id
                log_widget.append(
                    f"ℹ️ Previous {label} has finished "
                    f"(slurm job {job_id}, final state {job_state})")
                log_widget.append(f"   log: {log_file}")
                self._start_log_poll(target_dir, log_widget,
                                     task_type=task_type, once=True, reset=True)
                return
            elif head in ("PY_DEAD", "DEAD"):
                pid_txt = state[1] if len(state) > 1 else "?"
                log_widget.append(
                    f"ℹ️ Previous {label} has finished (pid {pid_txt})")
            else:
                log_widget.append(
                    f"ℹ️ Previous {label} has finished (no pid file found)")
            log_widget.append(f"   log: {log_file}")
            self._dump_remote_log(log_q, log_file, log_widget)
        elif head in ("PIPE_DONE", "PIPE_PARTIAL", "PIPE_FAILED"):
            done = head == "PIPE_DONE"
            self._set_tab_status(log_widget, 3 if done else 1)
            if not self._once(("finished", target_dir)):
                return
            prev_log = state[1] if len(state) > 1 else ""
            prev_log_q = shlex.quote(prev_log) if prev_log else ""
            if done:
                log_widget.append(
                    f"✅ Previous {label} finished (results found; "
                    f"no GUI session record)")
            elif head == "PIPE_FAILED":
                log_widget.append(
                    f"❌ Previous {label} failed (no GUI session record)")
            else:
                log_widget.append(
                    f"⏸ Previous {label} did not finish "
                    f"(no GUI session record)")
            if prev_log:
                log_widget.append(f"   log: {prev_log}")
                self._dump_remote_log(prev_log_q, prev_log, log_widget)

    def _dump_remote_log(self, log_q, log_file, log_widget, lines=200):
        """
        Replay the last N lines of a remote log into the output box (background).
        """
        cmd = f"tail -n {int(lines)} {log_q} 2>/dev/null"
        self._start_remote_query(
            cmd, self._on_log_dump,
            context=(log_widget, log_file, int(lines)),
            attr=self._dump_attr(log_file))

    @staticmethod
    def _dump_attr(log_file):
        """
        Thread reference per log file, so the two tabs cannot overwrite each other.
        """
        return "_dump_query" + (log_file or "").replace("/", "_")

    def _on_log_dump(self, ok, out, context, query):
        log_widget, log_file, lines = context
        if not self._is_current_query(query, getattr(self, self._dump_attr(log_file), None)):
            return
        log_widget.append(f"--- last {lines} lines of {log_file} ---")
        if not ok:
            log_widget.append(f"   ⚠️ cannot read log: {out}")
            return
        body = out.rstrip("\n")
        if body:
            for line in body.splitlines():
                log_widget.append(line)
        else:
            log_widget.append("   (log is empty)")
        log_widget.append("--- end of log ---")

    def _announce_job_id(self, target_dir, log_widget, attempt=1):
        """
        Read `.hvfinder.slurm.jobid` shortly after submission and log the job id.

        Non-blocking: the file is written by the job itself and may appear a little
        later, hence up to 3 retries. Direct mode never creates it.
        """
        if not self.run_mode_combo.currentText().startswith("Slurm"):
            return
        job_file = shlex.quote(f"{target_dir}/{SSHWorker.RUN_JOB_NAME}")
        delay = 3000 if attempt == 1 else 7000
        QTimer.singleShot(
            delay,
            lambda: self._start_remote_query(
                f"cat {job_file} 2>/dev/null", self._on_job_id,
                context=(target_dir, log_widget, attempt), attr="_jobid_query"))

    def _on_job_id(self, ok, out, context, query):
        if not self._is_current_query(query, getattr(self, "_jobid_query", None)):
            return
        target_dir, log_widget, attempt = context
        raw = out.strip().splitlines()[0].strip() if out.strip() else ""
        job_id = raw.split(";")[0].strip()
        if ok and job_id:
            self._slurm_jobs[target_dir] = job_id
            log_widget.append(f"🧾 Slurm job id: {job_id}")
            return
        if attempt < 3:
            self._announce_job_id(target_dir, log_widget, attempt + 1)

    def _slurm_jobs_command(self, job_name):
        """
        List same-named Slurm jobs of the current user (`jobid,STATE` per line).
        """
        script = (f"squeue -h --user=\"$USER\" --name={shlex.quote(job_name)} "
                  f"-o '%i,%T' 2>/dev/null")
        return f"bash -l -c {shlex.quote(script)}"

    def _guard_slurm_submission(self, job_name, log_widget, proceed):
        """
        Check for a duplicate Slurm job before submitting.

        `Nothing shown locally but the job is already queued` is the usual way to
        double-submit, so squeue is consulted first (background thread).
        """
        if not self.run_mode_combo.currentText().startswith("Slurm"):
            proceed()
            return
        self._start_remote_query(
            self._slurm_jobs_command(job_name), self._on_slurm_guard,
            context=(job_name, log_widget, proceed), attr="_guard_query")

    def _on_slurm_guard(self, ok, out, context, query):
        if not self._is_current_query(query, getattr(self, "_guard_query", None)):
            return
        job_name, log_widget, proceed = context
        if not ok:
            log_widget.append(f"⚠️ Cannot check Slurm queue ({out}); submitting anyway.")
            proceed()
            return

        jobs = []
        for line in (out or "").splitlines():
            if "," in line:
                jobid, state = line.strip().split(",", 1)
                jobs.append((jobid, state))
        if not jobs:
            proceed()
            return

        detail = "\n".join(f"  job {jid} ({state})" for jid, state in jobs)
        log_widget.append(f"⚠️ Slurm already has job(s) named {job_name}:")
        for line in detail.splitlines():
            log_widget.append(line)
        answer = QMessageBox.question(
            self, "Job already running",
            f'A Slurm job named "{job_name}" is already queued or running:\n\n'
            f"{detail}\n\nSubmit another one anyway?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if answer == QMessageBox.Yes:
            proceed()
        else:
            log_widget.append("⏹ Submission cancelled: existing job detected.")

    POLL_INTERVAL_MS = 3000
    POLL_IDLE_INTERVAL_MS = 15000
    POLL_IDLE_NOTICE_S = 600
    POLL_QUERY_TIMEOUT = 30
    POLL_EMPTY_CONFIRM = 3

    @staticmethod
    def _log_glob(task_type):
        """
        Log file name pattern produced by logger.py.
        """
        if task_type == "pathogen":
            return "pathogen_detect_*.log"
        return "results_pipeline_*.log"

    @staticmethod
    def _poll_attr(target_dir):
        return "_poll_query" + (target_dir or "").replace("/", "_")

    def _poll_command(self, target_dir, task_type, offset, out_offset=0):
        """
        Build one polling command: new log lines, Slurm's own stdout, and state.

        The output is machine-parsable and wraps the log body:
        @@FILE:<path>, lines, @@SIZE:<bytes>, @@OUT:<path>, lines, @@OUTSIZE:<bytes>,
        @@JOB:<jobid>, @@STATE:<RUNNING|PENDING|empty>.
        """
        dir_q = shlex.quote(target_dir)
        job_q = shlex.quote(f"{target_dir}/{SSHWorker.RUN_JOB_NAME}")
        out_q = shlex.quote(f"{target_dir}/{SSHWorker.SLURM_OUT_NAME}")
        off = int(offset)
        ooff = int(out_offset)
        head = ('if [ -n "$F" ]; then S=$(wc -c < "$F"); '
                f'if [ {off} -eq 0 ]; then tail -n 200 "$F"; '
                f'else tail -c +{off + 1} "$F" | head -c $((S-{off})); fi; '
                'echo; echo "@@SIZE:$S"; else echo "@@SIZE:0"; fi; ')
        out_head = ('if [ -f "$O" ]; then OS=$(wc -c < "$O"); '
                    f'if [ {ooff} -eq 0 ]; then tail -n 50 "$O"; '
                    f'else tail -c +{ooff + 1} "$O" | head -c $((OS-{ooff})); fi; '
                    'echo; echo "@@OUTSIZE:$OS"; else echo "@@OUTSIZE:0"; fi; ')
        script = (
            f"D={dir_q}; J={job_q}; O={out_q}; "
            f'F=$(ls -t "$D"/{self._log_glob(task_type)} 2>/dev/null | head -1); '
            'echo "@@FILE:$F"; '
            + head +
            'echo "@@OUT:$O"; '
            + out_head +
            'JID=""; if [ -f "$J" ]; then JID=$(head -1 "$J" 2>/dev/null); '
            'JID=${JID%%;*}; fi; '
            'echo "@@JOB:$JID"; '
            'ST=""; if [ -z "$JID" ]; then ST="NOJOB"; '
            'elif command -v squeue >/dev/null 2>&1; then '
            'ST=$(squeue -h -j "$JID" -o "%T" 2>/dev/null | head -1); '
            'else ST="UNKNOWN"; fi; '
            'echo "@@STATE:$ST"'
        )
        return f"bash -l -c {shlex.quote(script)}"

    def _start_log_poll(self, target_dir, log_widget, task_type="preprocess",
                        label=None, once=False, reset=False):
        """
        Start polling the remote log (Slurm mode).
        """
        if not target_dir or not self.host_edit.text().strip():
            return
        st = self._poll_state.setdefault(target_dir, {})
        if reset or "offset" not in st:
            st["offset"] = 0
            st["out_offset"] = 0
            st["file"] = ""
        st["task_type"] = task_type
        st["log_widget"] = log_widget
        st["label"] = label or ("pathogen detection" if task_type == "pathogen"
                                else "pre-processing pipeline")
        st["once"] = bool(once)
        st.setdefault("last_change", time.monotonic())
        st.setdefault("idle_bucket", 0)
        st.setdefault("empty_polls", 0)

        timer = self._poll_timers.get(target_dir)
        if timer is None:
            timer = QTimer(self)
            timer.setInterval(self.POLL_INTERVAL_MS)
            timer.timeout.connect(lambda d=target_dir: self._poll_tick(d))
            self._poll_timers[target_dir] = timer
        if not timer.isActive():
            timer.start()
        if not once:
            self._set_tab_status(log_widget, 2)
        self._poll_tick(target_dir)

    def _poll_active(self, target_dir):
        """
        Whether this target is already being polled.
        """
        timer = self._poll_timers.get(target_dir)
        if timer is None:
            return False
        try:
            return bool(timer.isActive())
        except RuntimeError:
            return False

    def _stop_log_poll(self, target_dir):
        timer = self._poll_timers.get(target_dir)
        if timer is not None:
            try:
                timer.stop()
            except RuntimeError:
                pass
        st = self._poll_state.get(target_dir)
        if st is not None:
            st["busy"] = False

    def _poll_tick(self, target_dir):
        st = self._poll_state.get(target_dir)
        timer = self._poll_timers.get(target_dir)
        if not st or timer is None or st.get("busy"):
            return
        st["busy"] = True
        self._start_remote_query(
            self._poll_command(target_dir, st["task_type"],
                               st.get("offset", 0), st.get("out_offset", 0)),
            self._on_poll_result, context=(target_dir,),
            attr=self._poll_attr(target_dir),
            timeout=self.POLL_QUERY_TIMEOUT)

    def _on_poll_result(self, ok, out, context, query):
        target_dir = context[0]
        st = self._poll_state.get(target_dir) or {}
        st["busy"] = False
        if not self._is_current_query(query, getattr(self, self._poll_attr(target_dir), None)):
            return

        log_widget = st.get("log_widget")
        if log_widget is None:
            return
        if not ok:
            if self._once(("poll-error", target_dir)):
                log_widget.append(
                    f"⚠️ Log polling hiccup (will keep retrying): {out}")
            return

        new_file = ""
        section = None
        app_lines, slurm_lines = [], []
        size = out_size = None
        job_id, state = "", None
        for line in out.splitlines():
            if line.startswith("@@"):
                name, _, value = line.partition(":")
                value = value.strip()
                section = None
                if name == "@@FILE":
                    new_file, section = value, "app"
                elif name == "@@OUT":
                    section = "out"
                elif name == "@@SIZE":
                    try:
                        size = int(value)
                    except ValueError:
                        size = None
                elif name == "@@OUTSIZE":
                    try:
                        out_size = int(value)
                    except ValueError:
                        out_size = None
                elif name == "@@JOB":
                    job_id = value
                elif name == "@@STATE":
                    state = value
                continue
            if section == "app":
                app_lines.append(line)
            elif section == "out":
                slurm_lines.append(line)

        if job_id:
            self._slurm_jobs[target_dir] = job_id

        if new_file and st.get("file") and new_file != st.get("file"):
            st["file"] = new_file
            st["offset"] = 0
            st["out_offset"] = 0
            QTimer.singleShot(0, lambda d=target_dir: self._poll_tick(d))
            return
        if new_file:
            st["file"] = new_file

        app_text = "\n".join(app_lines).strip("\n")
        if app_text:
            for line in app_text.splitlines():
                log_widget.append(line)
        if size is not None and size >= 0:
            st["offset"] = size

        slurm_text = "\n".join(slurm_lines).strip("\n")
        if slurm_text:
            for line in slurm_text.splitlines():
                log_widget.append(f"[slurm] {line}")
                marker = line.strip()
                if marker.startswith("[hvfinder] exit="):
                    try:
                        st["exit_code"] = int(marker.split("=", 1)[1])
                    except ValueError:
                        pass
        if out_size is not None and out_size >= 0:
            st["out_offset"] = out_size

        got_output = bool(app_text or slurm_text)
        now = time.monotonic()
        if got_output:
            st["last_change"] = now
            st["idle_bucket"] = 0
            timer = self._poll_timers.get(target_dir)
            if timer is not None and timer.interval() != self.POLL_INTERVAL_MS:
                timer.setInterval(self.POLL_INTERVAL_MS)
        else:
            idle_s = now - st.get("last_change", now)
            timer = self._poll_timers.get(target_dir)
            if timer is not None and idle_s >= 300 and \
                    timer.interval() != self.POLL_IDLE_INTERVAL_MS:
                timer.setInterval(self.POLL_IDLE_INTERVAL_MS)
            if state and state not in ("NOJOB", "UNKNOWN") and not st.get("file"):
                if self._once(("waiting", target_dir, state)):
                    jid = job_id or self._slurm_jobs.get(target_dir, "?")
                    log_widget.append(
                        f"⏳ Slurm job {jid} is {state}; waiting for the log file…")
            elif state in ("NOJOB", "UNKNOWN"):
                if self._once(("state-unknown", target_dir, state)):
                    jid = job_id or self._slurm_jobs.get(target_dir, "?")
                    log_widget.append(
                        f"⚠️ Cannot determine job state yet ({state}); "
                        f"job {jid}. Continuing to poll.")
            elif state and idle_s >= self.POLL_IDLE_NOTICE_S:
                bucket = int(idle_s // self.POLL_IDLE_NOTICE_S)
                if bucket > st.get("idle_bucket", 0):
                    st["idle_bucket"] = bucket
                    jid = job_id or self._slurm_jobs.get(target_dir, "?")
                    mins = int(idle_s // 60)
                    log_widget.append(
                        f"⏳ {st.get('label', 'task')} still RUNNING "
                        f"(job {jid}, no new log output for {mins} min)")

        if state is not None and state == "":
            st["empty_polls"] = st.get("empty_polls", 0) + 1
        else:
            st["empty_polls"] = 0
        finished = st.get("empty_polls", 0) >= self.POLL_EMPTY_CONFIRM
        if (finished and st.get("expect_rc") and st.get("exit_code") is None
                and not st.get("final_check")):
            st["final_check"] = True
            st["empty_polls"] = 0
            return
        if finished or st.get("once"):
            self._stop_log_poll(target_dir)
            job_running = bool(state) and state not in ("NOJOB", "UNKNOWN")
            self._set_tab_status(
                log_widget, 2 if (job_running and not finished) else 3)
            if finished and self._once(("finished", target_dir)):
                self._announce_task_finished(target_dir, log_widget, job_id,
                                            st.get("file") or new_file)

    def _announce_task_finished(self, target_dir, log_widget, job_id, log_file):
        st = self._poll_state.get(target_dir) or {}
        label = st.get("label", "task")
        rc = st.get("exit_code")
        elapsed = self._elapsed_text(log_widget, target_dir)
        log_widget.append("=" * 50)
        if rc in (None, 0):
            log_widget.append(f"✅ {label} finished"
                              + (f" (slurm job {job_id})" if job_id else "")
                              + (f"   ({elapsed})" if elapsed else ""))
        else:
            rc_txt = "unknown" if rc < 0 else str(rc)
            log_widget.append(f"❌ {label} failed (slurm job {job_id}, "
                              f"exit code {rc_txt})"
                              + (f"   ({elapsed})" if elapsed else ""))
            log_widget.append("   Check the log above for the first failing step.")
        if log_file:
            log_widget.append(f"Full log: {log_file}")
        log_widget.append("=" * 50)

    def on_tab_changed(self, index):
        """
        On returning to a tab, detect and re-attach a remote task that still runs.
        """
        if not self.current_project_dir:
            return
        if index == 1:
            self.reattach_running_task(f"{self.current_project_dir}/preprocess",
                                       self.pipe_log_text, "pre-processing pipeline",
                                       task_type="preprocess")
        elif index == 2:
            out_dir = (f"{self.current_project_dir}/"
                       f"{self.pathogen_output_folder_edit.text().strip()}")
            self.reattach_running_task(
                out_dir, self.pathogen_log_text, "pathogen detection",
                task_type="pathogen",
                fallback_dir=f"{self.current_project_dir}/pathogen")

    def _confirm_and_stop(self, target_dir, log_widget, label):
        """
        After confirmation, kill the remote process group recorded for target_dir
        (the setsid session).
        """
        if not target_dir:
            QMessageBox.warning(self, "No Target",
                                "No task directory found; cannot locate the remote process.")
            return
        pid_file = f"{target_dir}/{SSHWorker.RUN_PID_NAME}"
        job_file = f"{target_dir}/{SSHWorker.RUN_JOB_NAME}"
        job_id = self._slurm_jobs.get(target_dir, "")
        job_txt = f"{job_id} (from cache)" if job_id else "unknown"
        answer = QMessageBox.question(
            self, f"Stop {label}",
            f"Stop the remote {label} that is currently running?\n\n"
            f"PID file: {pid_file}\n"
            f"Slurm job: {job_txt}\n"
            f"Finished steps and intermediate files are kept, so you can resume "
            f"later (unless --no-resume was used).\n\n"
            f"This cannot be undone.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if answer != QMessageBox.Yes:
            log_widget.append("ℹ️ Stop cancelled.")
            return

        self._stop_log_poll(target_dir)
        self._set_tab_status(log_widget, 1)
        pid_q = shlex.quote(pid_file)
        job_q = shlex.quote(job_file)
        cmd = (
            f"if [ -f {job_q} ]; then "
            f"JID=$(cat {job_q} 2>/dev/null); "
            f"if [ -n \"$JID\" ]; then "
            f"echo \"[hvfinder] scancel $JID\"; "
            f"scancel \"$JID\" 2>/dev/null || echo '[hvfinder] scancel failed'; "
            f"fi; fi; "
            f"if [ -f {pid_q} ]; then "
            f"P=$(cat {pid_q}); "
            f"echo \"[hvfinder] stopping process group $P\"; "
            f"kill -TERM -- -$P 2>/dev/null || kill -TERM $P 2>/dev/null; "
            f"sleep 3; "
            f"kill -KILL -- -$P 2>/dev/null || true; "
            f"rm -f {pid_q}; "
            f"echo '[hvfinder] stopped'; "
            f"else echo '[hvfinder] no pid file (task not running or already finished)'; fi"
        )
        log_widget.append(f"⏹ Requested stop for: {label}")
        self.run_ssh_task(command=cmd, log_widget=log_widget)

    def on_stop_pipeline(self):
        if not self.current_project_dir:
            QMessageBox.warning(self, "Path Error", "No project folder selected!")
            return
        self._confirm_and_stop(f"{self.current_project_dir}/preprocess",
                               self.pipe_log_text, "pre-processing pipeline")

    def on_stop_pathogen(self):
        if not self.current_project_dir:
            QMessageBox.warning(self, "Path Error", "No project folder selected!")
            return
        out_dir = f"{self.current_project_dir}/{self.pathogen_output_folder_edit.text().strip()}"
        self._confirm_and_stop(out_dir, self.pathogen_log_text, "pathogen detection")

    def run_ssh_task(self, command=None, upload_cfg=None, target_dir=None, task_type="cmd", log_widget=None, finished_callback=None, detached=False, upload_files=None):
        self.save_settings()
        self.worker = SSHWorker(
            self.host_edit.text().strip(), int(self.port_edit.text()),
            self.user_edit.text().strip(), self.pwd_edit.text(),
            command=command, upload_cfg=upload_cfg, target_dir=target_dir,
            task_type=task_type, detached=detached, upload_files=upload_files
        )
        self.worker.log_widget = log_widget or self.log_text
        if log_widget:
            self.worker.log_signal.connect(lambda msg: self.safe_log(msg, log_widget))
        else:
            self.worker.log_signal.connect(lambda msg: self.safe_log(msg, self.log_text))
            
        self.worker.status_signal.connect(
            lambda code, w=self.worker.log_widget: self.update_status_ui(
                code, self._status_label_for(w)))
        
        if finished_callback:
            self.worker.finished_signal.connect(finished_callback)
        else:
            self.worker.finished_signal.connect(self.on_task_finished)
        
        self.worker.task_type_ext = task_type
        self.worker.start()

    def get_ui_as_yaml(self):
        """
        Write the UI values back into the nested dict and serialise it as YAML.
        """
        import copy
        new_cfg = copy.deepcopy(self.pipeline_cfg)

        def _apply(mapping):
            for path_tuple, edit in mapping.items():
                val = edit.text()
                if val.lower() == 'true': val = True
                elif val.lower() == 'false': val = False
                elif val.isdigit(): val = int(val)

                temp = new_cfg
                for key in path_tuple[:-1]:
                    temp = temp.setdefault(key, {})
                temp[path_tuple[-1]] = val

        _apply(self.ui_param_map)
        _apply(getattr(self, "pathogen_param_map", {}))   # Tab 3（nr DIAMOND）

        if hasattr(self, "cores_spin"):
            new_cfg.setdefault('global', {})['threads'] = int(self.cores_spin.value())
        if hasattr(self, "pathogen_node_mem_spin"):
            new_cfg.setdefault('diamond', {}).setdefault('params', {})['node_mem_gb'] = \
                int(self.pathogen_node_mem_spin.value())
        
        return yaml.dump(new_cfg, sort_keys=False, allow_unicode=True)

    def _adv_raw(self, *path, default=""):
        """
        Read an Advanced-panel value by its config2.yaml key path.
        """
        edit = self.ui_param_map.get(tuple(path))
        value = edit.text().strip() if edit is not None else ""
        return value or default

    def db_settings_block(self):
        """
        Collect the database paths shown in the UI (same source as [ADVANCED] /
        [PATHOGEN]).
        """
        return {
            "diamond_db": self.pathogen_diamond_db_edit.text().strip(),
            "bowtie2_index": self._adv_raw("bowtie2", "params", "index_path"),
            "blastn_index": self.get_virus_index_path(),
            "virus_taxids": self.get_taxid_file_path(),
        }

    def build_server_config_yaml(self):
        """
        Build the full server-side config2.yaml (including the `databases` block).
        """
        cfg = yaml.safe_load(self.get_ui_as_yaml()) or {}
        cfg.setdefault("databases", {}).update(self.db_settings_block())
        return yaml.dump(cfg, sort_keys=False, allow_unicode=True)

    def sync_db_to_server(self, log_widget=None):
        """
        Write the database settings into `<server>/configs/config2.yaml`.
        """
        log_widget = log_widget or getattr(self, "log_text", None)
        if not self.host_edit.text().strip():
            if log_widget is not None:
                log_widget.append("⚠️ Server config not synced: SSH host is empty.")
            return False
        base = self._script_base or self.ensure_script_base()
        if not base:
            if log_widget is not None:
                log_widget.append(
                    "⚠️ Server config not synced: cannot resolve the remote "
                    "hvFinder directory (get_script_base).")
            return False

        worker = getattr(self, "_config_sync_worker", None)
        if worker is not None and self._thread_running(worker):
            if log_widget is not None:
                log_widget.append("ℹ️ Server config sync is already running.")
            return False

        remote_dir = f"{base}/configs"
        content = self.build_server_config_yaml()
        w = SSHWorker(
            self.host_edit.text().strip(), int(self.port_edit.text() or 22),
            self.user_edit.text().strip(), self.pwd_edit.text(),
            command=None, target_dir=remote_dir, task_type="sync-config",
            upload_files={"config2.yaml": content})
        w.log_widget = log_widget
        if log_widget is not None:
            w.log_signal.connect(lambda msg: self.safe_log(msg, log_widget))
            w.finished_signal.connect(lambda ok, msg: log_widget.append(
                f"✅ Server config synced: {remote_dir}/config2.yaml "
                f"(databases block updated)"
                if ok else f"❌ Server config sync failed: {msg}"))
        self._config_sync_worker = w
        w.start()
        return True

    STATUS_TEXT = {
        0: ("● Disconnected", "#7F8C8D"),
        1: ("● Online", "#27AE60"),
        2: ("● Running Task", "#E67E22"),
    }
    TAB_STATUS_TEXT = {
        0: ("● Idle", "#7F8C8D"),
        1: ("● Ready", "#27AE60"),
        2: ("● Running Task", "#E67E22"),
        3: ("● Finished", "#2980B9"),
    }

    def _status_label_for(self, log_widget):
        """
        Map a log box to the status label of its tab.
        """
        if log_widget is getattr(self, "pipe_log_text", None):
            return getattr(self, "pipe_status_label", None)
        if log_widget is getattr(self, "pathogen_log_text", None):
            return getattr(self, "pathogen_status_label", None)
        return None

    def update_status_ui(self, code, label=None):
        """
        Update the status indicator (Tab 1: connection, Tabs 2/3: task state).
        """
        text, color = self.STATUS_TEXT.get(code, self.STATUS_TEXT[0])
        self.status_label.setText(text)
        self.status_label.setStyleSheet(f"color: {color}; font-weight: bold;")
        if label is not None:
            tab_text, tab_color = self.TAB_STATUS_TEXT.get(
                code, self.TAB_STATUS_TEXT[0])
            label.setText(tab_text)
            label.setStyleSheet(f"color: {tab_color}; font-weight: bold;")

    def _set_tab_status(self, log_widget, code):
        """
        Update the status indicator of the tab that owns this log box.
        """
        self.update_status_ui(code, self._status_label_for(log_widget))

    def set_tab_status(self, code):
        """
        Update only Tabs 2/3 (leave Tab 1's connection status alone).
        """
        text, color = self.TAB_STATUS_TEXT.get(code, self.TAB_STATUS_TEXT[0])
        for lbl in (getattr(self, "pipe_status_label", None),
                    getattr(self, "pathogen_status_label", None)):
            if lbl is not None:
                lbl.setText(text)
                lbl.setStyleSheet(f"color: {color}; font-weight: bold;")

    def safe_log(self, text, widget):
        stripped = (text or "").strip()
        if stripped.startswith("[hvfinder] exit="):
            try:
                self._task_rc[id(widget)] = int(stripped.split("=", 1)[1])
            except ValueError:
                self._task_rc[id(widget)] = -1
        widget.append(text)
        widget.moveCursor(QTextCursor.End)

    def _copy_exec_value(self, src, dst, kind):
        """
        Copy src into dst; kind is mode (QComboBox), spin (QSpinBox) or text
        (QLineEdit).
        """
        if kind == "mode":
            idx = dst.findText(src.currentText())
            if idx >= 0:
                dst.setCurrentIndex(idx)
        elif kind == "spin":
            dst.setValue(src.value())
        elif dst.text() != src.text():
            dst.setText(src.text())

    def _make_mirror(self, src, dst, kind):
        """
        Build a slot that mirrors src into dst, guarded against recursion.
        """
        def push(*_args):
            if getattr(self, "_syncing_exec", False):
                return
            self._syncing_exec = True
            try:
                self._copy_exec_value(src, dst, kind)
            finally:
                self._syncing_exec = False
        return push

    def link_execution_controls(self):
        """
        Keep Run Mode, CPU cores and extra Slurm arguments in sync across tabs.

        Both tabs drive the same server and the same resources, so separate copies
        would drift (64 cores on Tab 2, 8 on Tab 3). Tab 2 is the initial source of
        truth.
        """
        self._syncing_exec = False

        pairs = (
            (self.run_mode_combo, self.pathogen_run_mode_combo, "mode"),
            (self.pathogen_run_mode_combo, self.run_mode_combo, "mode"),
            (self.cores_spin, self.pathogen_cores_spin, "spin"),
            (self.pathogen_cores_spin, self.cores_spin, "spin"),
            (self.slurm_extra_edit, self.pathogen_slurm_extra_edit, "text"),
            (self.pathogen_slurm_extra_edit, self.slurm_extra_edit, "text"),
        )
        for src, dst, kind in pairs:
            if kind == "mode":
                signal = src.currentTextChanged
            elif kind == "spin":
                signal = src.valueChanged
            else:
                signal = src.textEdited
            signal.connect(self._make_mirror(src, dst, kind))

        self.run_mode_combo.currentTextChanged.connect(
            lambda *_: self._sync_execution_rows())

        self._syncing_exec = True
        try:
            self._copy_exec_value(self.run_mode_combo, self.pathogen_run_mode_combo, "mode")
            self._copy_exec_value(self.cores_spin, self.pathogen_cores_spin, "spin")
            self._copy_exec_value(self.slurm_extra_edit,
                                  self.pathogen_slurm_extra_edit, "text")
        finally:
            self._syncing_exec = False

    def _sync_execution_rows(self):
        """
        Direct mode does not use srun arguments, so hide that row.
        """
        visible = self.run_mode_combo.currentText().startswith("Slurm")
        for label, field in getattr(self, "_slurm_rows", []):
            if label is not None:
                label.setVisible(visible)
            field.setVisible(visible)

    def try_auto_connect(self):
        """
        Log in automatically when credentials are saved; report either way.
        """
        host = self.host_edit.text().strip()
        user = self.user_edit.text().strip()
        if not host or not user:
            self.log_text.append("ℹ️ No saved SSH credentials — auto-login skipped.")
            return
        try:
            port = int(self.port_edit.text().strip() or 22)
        except ValueError:
            port = 22

        self.log_text.append(f"🔌 Auto-login: connecting to {user}@{host}:{port} ...")
        self.update_status_ui(2)
        self._ssh_probe = SSHProbe(host, port, user, self.pwd_edit.text(), timeout=8)
        self._ssh_probe.result_signal.connect(self._on_autoconnect_result)
        self._ssh_probe.start()

    def refresh_root_folders(self):
        """
        List first-level folders under the project root (called after auto-login,
        asynchronously).
        """
        root = (self.search_dir_edit.text() or "").strip().rstrip("/")
        if not root:
            return
        self.log_text.append(f"📁 Project root: {root}")
        cmd = (f"find {shlex.quote(root)} -mindepth 1 -maxdepth 1 -type d "
               f"-printf '%f\\n' 2>/dev/null | sort")
        self._start_remote_query(cmd, self._on_root_folders,
                                 context=(self.log_text, root), attr="_root_query")

    def _on_root_folders(self, ok, out, context, query):
        if not self._is_current_query(query, getattr(self, "_root_query", None)):
            return
        log_widget, root = context
        if not ok:
            log_widget.append(f"   ⚠️ Cannot list folders: {out}")
            return
        folders = [line.strip() for line in out.splitlines() if line.strip()]
        if not folders:
            log_widget.append("   (empty) no sub-directories")
            return
        log_widget.append(f"   {len(folders)} folder(s):")
        for name in folders:
            log_widget.append(f"     - {name}")

    def refresh_tab_statuses(self):
        """
        Probe the real task state of Tabs 2/3 in the background (start polling if
        needed).
        """
        if not self.current_project_dir:
            return
        self.reattach_running_task(f"{self.current_project_dir}/preprocess",
                                   self.pipe_log_text, "pre-processing pipeline",
                                   task_type="preprocess")
        out_dir = (f"{self.current_project_dir}/"
                   f"{self.pathogen_output_folder_edit.text().strip()}")
        self.reattach_running_task(out_dir, self.pathogen_log_text,
                                   "pathogen detection", task_type="pathogen",
                                   fallback_dir=f"{self.current_project_dir}/pathogen")

    def _verify_active_project(self):
        """
        Confirm in the background that the last active project still exists.
        """
        path = (self.current_project_dir or "").strip()
        if not path:
            return
        cmd = f"[ -d {shlex.quote(path)} ] && echo EXISTS || echo MISSING"
        self._start_remote_query(cmd, self._on_verify_project,
                                 context=(path, self.log_text), attr="_verify_query")

    def _on_verify_project(self, ok, out, context, query):
        if not self._is_current_query(query, getattr(self, "_verify_query", None)):
            return
        path, log_widget = context
        if not ok:
            return
        if out.strip() == "EXISTS":
            log_widget.append(f"📁 Restored last project: {path}")
        elif self.current_project_dir == path:
            self._set_active_project("")
            log_widget.append(
                f"⚠️ Last project no longer exists on the server: {path}")

    def _on_autoconnect_result(self, ok, message):
        if ok:
            self.update_status_ui(1)
            self.set_tab_status(1)      # Tab 2/3：Ready
            self.log_text.append(
                f"✅ Auto-login succeeded: {self.user_edit.text().strip()}@"
                f"{self.host_edit.text().strip()}")
            self.refresh_root_folders()
            self._verify_active_project()
            self.refresh_tab_statuses()
        else:
            self.update_status_ui(0)
            self.set_tab_status(0)      # Tab 2/3：Idle
            self.log_text.append("❌ Auto-login failed — check host / port / username / password.")
            self.log_text.append(f"   reason: {message}")
        self.log_text.moveCursor(QTextCursor.End)

    def setup_autosave(self):
        """
        Save 500 ms after an edit, so parameter changes are not lost.
        """
        self._autosave_timer = QTimer(self)
        self._autosave_timer.setSingleShot(True)
        self._autosave_timer.setInterval(500)
        self._autosave_timer.timeout.connect(lambda: self.save_settings(quiet=True))

        edits = [self.host_edit, self.port_edit, self.user_edit, self.pwd_edit,
                 self.search_dir_edit, self.fq1_edit, self.fq2_edit,
                 self.slurm_extra_edit, self.pathogen_slurm_extra_edit]
        edits += [w for w in self.ui_param_map.values() if isinstance(w, QLineEdit)]
        edits += [w for w in getattr(self, "pathogen_param_map", {}).values()
                  if isinstance(w, QLineEdit)]
        for name in ("diamond_pident_edit", "diamond_evalue_edit",
                     "blastn_identity_edit", "blastn_evalue_edit", "blastn_coverage_edit",
                     "pathogen_output_folder_edit"):
            w = getattr(self, name, None)
            if isinstance(w, QLineEdit):
                edits.append(w)
        for w in edits:
            w.editingFinished.connect(self.schedule_autosave)

        for spin in (getattr(self, "cores_spin", None),
                     getattr(self, "pathogen_cores_spin", None),
                     getattr(self, "pathogen_node_mem_spin", None)):
            if spin is not None:
                spin.valueChanged.connect(self.schedule_autosave)

        for combo in (self.run_mode_combo, getattr(self, "pathogen_run_mode_combo", None)):
            if combo is not None:
                combo.currentIndexChanged.connect(self.schedule_autosave)
        for cb in (self.remember_cb, self.save_pipeline_cb,
                   self.pipeline_verbose_cb, self.pipeline_resume_cb):
            cb.toggled.connect(self.schedule_autosave)

    def schedule_autosave(self, *args):
        if getattr(self, "_loading", False):
            return
        if hasattr(self, "_autosave_timer"):
            self._autosave_timer.start()

    def stop_background_threads(self):
        """
        Wait for background queries, so no QThread is destroyed while running
        (which aborts the process).
        """
        for timer in list(getattr(self, "_poll_timers", {}).values()):
            try:
                timer.stop()
            except RuntimeError:
                pass
        attrs = ["_ssh_probe", "_reattach_query", "_root_query",
                 "_verify_query", "_script_base_query", "_jobid_query",
                 "_guard_query", "_config_sync_worker"]
        attrs += [k for k in list(self.__dict__) if k.startswith("_dump_query")]
        attrs += [k for k in list(self.__dict__) if k.startswith("_poll_query")]
        for attr in attrs:
            worker = getattr(self, attr, None)
            if not self._thread_running(worker):
                setattr(self, attr, None)
                continue
            try:
                worker.result_signal.disconnect()
            except (RuntimeError, TypeError):
                pass
            try:
                worker.wait(9000)
            except RuntimeError:
                pass

    def closeEvent(self, event):
        try:
            self.save_settings(quiet=True)
        except Exception:
            pass
        self.stop_background_threads()
        super().closeEvent(event)

    def save_settings(self, quiet=False):
        if hasattr(self, "remember_cb") and not self.remember_cb.isChecked():
            if not quiet:
                self.statusBar().showMessage(
                    "Settings not saved (Save Settings is unchecked)", 3000)
            return

        config = configparser.ConfigParser()
        config['SSH'] = {
            'host': self.host_edit.text(), 'port': self.port_edit.text(),
            'user': self.user_edit.text(), 'pwd': self.pwd_edit.text(),
            'base': self.search_dir_edit.text(),
            'active_project': self.current_project_dir or '',
        }
        if hasattr(self, 'save_pipeline_cb') and self.save_pipeline_cb.isChecked():
            config['PIPELINE'] = {
                'fq1': self.fq1_edit.text(),
                'fq2': self.fq2_edit.text(),
                'run_mode': self.run_mode_combo.currentText(),
                'threads': str(self.cores_spin.value()),
                'node_mem_gb': str(self.pathogen_node_mem_spin.value()),
                'slurm_extra': self.slurm_extra_edit.text().strip(),
                'verbose': str(self.pipeline_verbose_cb.isChecked()),
                'resume': str(self.pipeline_resume_cb.isChecked())
            }
            config['ADVANCED'] = {}
            for path_tuple, edit in self.ui_param_map.items():
                key = ".".join(path_tuple)
                config['ADVANCED'][key] = edit.text()
            for path_tuple, edit in getattr(
                    self, "pathogen_param_map", {}).items():
                key = ".".join(path_tuple)
                config['ADVANCED'][key] = edit.text()
        
        if hasattr(self, 'save_pathogen_cb') and self.save_pathogen_cb.isChecked():
            config['PATHOGEN'] = {
                'taxid_file_preset': self.taxid_file_combo.currentText(),
                'virus_index_preset': self.virus_index_combo.currentText(),
                'output_folder': self.pathogen_output_folder_edit.text(),
                'diamond_pident': self.diamond_pident_edit.text(),
                'diamond_evalue': self.diamond_evalue_edit.text(),
                'blastn_identity': self.blastn_identity_edit.text(),
                'blastn_evalue': self.blastn_evalue_edit.text(),
                'blastn_coverage': self.blastn_coverage_edit.text(),
                'verbose': str(self.verbose_cb.isChecked()),
                'resume': str(self.resume_cb.isChecked())
            }
        
        try:
            with open(self.settings_file, 'w', encoding='utf-8') as f:
                config.write(f)
        except OSError as exc:
            if not quiet:
                QMessageBox.warning(self, "Save Error",
                                    f"Cannot write {self.settings_file}:\n{exc}")
            return
        if not quiet:
            self.statusBar().showMessage(f"Settings saved: {self.settings_file}", 3000)

    def load_saved_settings(self):
        self._loading = True
        if os.path.exists(self.settings_file):
            config = read_ini_file(self.settings_file)
            if 'SSH' in config:
                self.host_edit.setText(config['SSH'].get('host', ''))
                self.port_edit.setText(config['SSH'].get('port', '22'))
                self.user_edit.setText(config['SSH'].get('user', ''))
                self.pwd_edit.setText(config['SSH'].get('pwd', ''))
                self.search_dir_edit.setText(config['SSH'].get('base', ''))
                saved_project = config['SSH'].get('active_project', '').strip()
                if saved_project:
                    self._set_active_project(saved_project)
            
            if 'PIPELINE' in config:
                self.fq1_edit.setText(config['PIPELINE'].get('fq1', ''))
                self.fq2_edit.setText(config['PIPELINE'].get('fq2', ''))
                run_mode = config['PIPELINE'].get('run_mode', 'Direct Python')
                idx = self.run_mode_combo.findText(run_mode)
                if idx >= 0:
                    self.run_mode_combo.setCurrentIndex(idx)
                try:
                    self.cores_spin.setValue(
                        int(config['PIPELINE'].get('threads', self.cores_spin.value())))
                except ValueError:
                    pass
                self.slurm_extra_edit.setText(config['PIPELINE'].get('slurm_extra', ''))
                try:
                    self.pathogen_node_mem_spin.setValue(
                        int(config['PIPELINE'].get('node_mem_gb',
                                                   self.pathogen_node_mem_spin.value()) or 0))
                except ValueError:
                    pass
                
                if hasattr(self, 'pipeline_verbose_cb') and 'verbose' in config['PIPELINE']:
                    self.pipeline_verbose_cb.setChecked(config['PIPELINE']['verbose'] == 'True')
                if hasattr(self, 'pipeline_resume_cb') and 'resume' in config['PIPELINE']:
                    self.pipeline_resume_cb.setChecked(config['PIPELINE']['resume'] == 'True')
            
            if 'ADVANCED' in config:
                for key, value in config['ADVANCED'].items():
                    path_tuple = tuple(key.split("."))
                    if path_tuple in self.ui_param_map:
                        self.ui_param_map[path_tuple].setText(value)
                    elif path_tuple in getattr(self, "pathogen_param_map", {}):
                        self.pathogen_param_map[path_tuple].setText(value)
            
            if 'PATHOGEN' in config:
                taxid_preset = config['PATHOGEN'].get('taxid_file_preset', '')
                if taxid_preset:
                    idx = self.taxid_file_combo.findText(taxid_preset)
                    if idx >= 0:
                        self.taxid_file_combo.setCurrentIndex(idx)
                    else:
                        self.taxid_file_combo.setCurrentText(taxid_preset)
                else:
                    old_taxid = config['PATHOGEN'].get('taxid_file', '')
                    if old_taxid:
                        self.taxid_file_combo.setCurrentText(old_taxid)
                
                virus_preset = config['PATHOGEN'].get('virus_index_preset', '')
                if virus_preset:
                    idx = self.virus_index_combo.findText(virus_preset)
                    if idx >= 0:
                        self.virus_index_combo.setCurrentIndex(idx)
                    else:
                        self.virus_index_combo.setCurrentText(virus_preset)
                else:
                    old_virus = config['PATHOGEN'].get('virus_index', '')
                    if old_virus:
                        self.virus_index_combo.setCurrentText(old_virus)
                
                self.pathogen_output_folder_edit.setText(config['PATHOGEN'].get('output_folder', 'pathogen_results'))
                self.diamond_pident_edit.setText(config['PATHOGEN'].get('diamond_pident', '50'))
                self.diamond_evalue_edit.setText(config['PATHOGEN'].get('diamond_evalue', '1e-10'))
                self.blastn_identity_edit.setText(config['PATHOGEN'].get('blastn_identity', '70'))
                self.blastn_evalue_edit.setText(config['PATHOGEN'].get('blastn_evalue', '1e-4'))
                self.blastn_coverage_edit.setText(config['PATHOGEN'].get('blastn_coverage', '70'))
                
                if 'verbose' in config['PATHOGEN']:
                    self.verbose_cb.setChecked(config['PATHOGEN']['verbose'] == 'True')
                if 'resume' in config['PATHOGEN']:
                    self.resume_cb.setChecked(config['PATHOGEN']['resume'] == 'True')

        self._loading = False


def main():
    """
    GUI entry point (also exposed as the `hvfinder-gui` console script).
    """
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    code = app.exec()
    window.stop_background_threads()
    return code


if __name__ == "__main__":
    sys.exit(main())