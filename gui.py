# -*- coding: utf-8 -*-

"""
Author: Ying Huang
Date: 2026-01-27
Description: hvFinder Pro 
"""

# 屏蔽 libpng iCCP 警告（必须放在 PySide6 导入之前）
import os
os.environ["QT_LOGGING_RULES"] = "*.debug=false;qt.qpa.xcb=false"

import sys
import yaml
import paramiko
import configparser
from datetime import datetime
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QTextEdit, QTabWidget, QGroupBox, 
    QFormLayout, QScrollArea, QCheckBox, QMessageBox, QInputDialog,
    QComboBox, QTableWidget, QTableWidgetItem, QHeaderView, QDialog,
    QDialogButtonBox, QFileDialog
)
from PySide6.QtCore import Qt, Signal, QThread
from PySide6.QtGui import QTextCursor

# --- 强化型 SSH/SFTP 执行线程 ---
class SSHWorker(QThread):
    log_signal = Signal(str)
    status_signal = Signal(int)
    # 扩展信号：成功标志, 消息内容
    finished_signal = Signal(bool, str) 

    def __init__(self, host, port, user, pwd, command=None, upload_cfg=None, target_dir=None, task_type="cmd"):
        super().__init__()
        self.host = host
        self.port = port
        self.user = user
        self.pwd = pwd
        self.command = command
        self.upload_cfg = upload_cfg
        self.target_dir = target_dir
        self.task_type = task_type # "cmd" 代表执行命令, "mkdir" 代表创建目录

    def run(self):
        ssh = None
        try:
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(self.host, port=self.port, username=self.user, password=self.pwd, timeout=10)
            
            self.status_signal.emit(1) # 已连接
            self.emit_log(f"✅ Connected to {self.host}")

            # 任务类型 1：创建目录
            if self.task_type == "mkdir":
                ssh.exec_command(f"mkdir -p {self.target_dir}")
                self.emit_log(f"📁 Created workspace: {self.target_dir}")
                self.finished_signal.emit(True, "CREATED")
                return

            # 任务类型 2：常规 Pipeline 执行 (严格顺序：创建目录 -> 上传配置 -> 执行命令)
            if self.target_dir:
                # 1. 确保目录存在
                stdin, stdout, stderr = ssh.exec_command(f"mkdir -p {self.target_dir}")
                stdout.channel.recv_exit_status() # 等待目录创建指令执行完毕
                
                # 2. 如果有配置需要上传，则执行 SFTP
                if self.upload_cfg:
                    sftp = ssh.open_sftp()
                    remote_cfg_path = f"{self.target_dir}/config.yaml"
                    with sftp.file(remote_cfg_path, "w") as f:
                        f.write(self.upload_cfg)
                    sftp.close()
                    self.emit_log(f"📄 Config uploaded to: {remote_cfg_path}")

            # 3. 最后运行命令，确保 config.yaml 已经就绪
            if self.command:
                self.status_signal.emit(2) # 忙碌
                # 切换到目标目录执行，确保脚本能找到本地的 config.yaml
                final_cmd = f"cd {self.target_dir} && {self.command}" if self.target_dir else self.command
                self.emit_log(f"🚀 Running command: {final_cmd}")
                
                stdin, stdout, stderr = ssh.exec_command(f"bash -l -c '{final_cmd}' 2>&1", get_pty=True)
                for line in stdout:
                    self.log_signal.emit(line.strip())
                
                # 阻塞直到命令完成
                stdout.channel.recv_exit_status()
                self.status_signal.emit(1) # 恢复在线
            
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

# --- 预设配置（相对路径） ---
PRESET_VIRUS_INDEXES = {
    "Test Database": "data/test_idx/test",
    "High Pathogenic Virus DB": "data/high_pathogenic_virus_db/virus_sequences_unique_ids.fa",
    "Custom": ""
}

PRESET_TAXID_FILES = {
    "Default Virus TaxIDs": "data/default_virus_taxids.txt",
    "High Pathogenic Virus TaxIDs": "data/virus_taxids_high_pathogenic.txt",
    "Custom": ""
}

# --- 主窗口界面 ---
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("hvFinder Pro - Management Console")
        self.resize(1000, 850)
        
        self.settings_file = "settings.ini"
        
        # 指定 configs/config2.yaml 的路径
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(script_dir, "configs", "config2.yaml")
        self.pipeline_cfg = self.load_local_template(config_path)
        
        self.current_project_dir = "" # 用于跨标签页同步的项目路径
        
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        self.init_ssh_tab()
        self.init_pipeline_tab()
        self.init_pathogen_tab()

        self.load_saved_settings()

        self.apply_button_styles()

    def apply_button_styles(self):
        """仅修改按钮背景颜色，保持原生形状"""
        blue_style = """
            QPushButton {
                background-color: #2E86C1;
                color: white;
            }
            QPushButton:hover {
                background-color: #3498DB;
            }
            QPushButton:pressed {
                background-color: #1A5276;
            }
        """
        # 为所有按钮应用样式
        self.list_dir_btn.setStyleSheet(blue_style)
        self.mkdir_btn.setStyleSheet(blue_style)
        self.toggle_btn.setStyleSheet(blue_style)
        self.run_btn.setStyleSheet(blue_style)
        
        # Pathogen tab buttons
        self.pathogen_toggle_btn.setStyleSheet(blue_style)
        self.pathogen_run_btn.setStyleSheet(blue_style)
        self.view_results_btn.setStyleSheet(blue_style)
        self.export_csv_btn.setStyleSheet(blue_style)

    def load_local_template(self, pth):
        """读取指定路径的 YAML 模板"""
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

        # 1. 服务器验证
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

        # 2. 搜索目录与项目创建 (Search Directory)
        dir_group = QGroupBox("Search Directory & Project Setup")
        dir_form = QFormLayout(dir_group)
        
        self.search_dir_edit = QLineEdit("/home/data/mNGS_projects")
        
        btn_layout = QHBoxLayout()
        self.list_dir_btn = QPushButton("🔍 List Files")
        self.list_dir_btn.clicked.connect(self.on_list_remote)
        self.mkdir_btn = QPushButton("📁 Create/Select Project Folder")
        self.mkdir_btn.setStyleSheet("background-color: #D4E6F1; font-weight: bold;")
        self.mkdir_btn.clicked.connect(self.on_create_project)
        
        btn_layout.addWidget(self.list_dir_btn)
        btn_layout.addWidget(self.mkdir_btn)

        self.sync_info_label = QLabel("Active Project: None")
        self.sync_info_label.setStyleSheet("color: #2980B9; font-style: italic;")

        dir_form.addRow("Target Path:", self.search_dir_edit)
        dir_form.addRow(btn_layout)
        dir_form.addRow(self.sync_info_label)
        layout.addWidget(dir_group)

        # 3. 日志与状态
        status_layout = QHBoxLayout()
        self.status_label = QLabel("● Disconnected")
        self.status_label.setStyleSheet("color: #7F8C8D; font-weight: bold;")
        status_layout.addStretch()
        status_layout.addWidget(self.status_label)
        layout.addLayout(status_layout)

        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setStyleSheet("background-color: #0c0c0c; color: #00FF41; font-family: 'Consolas';")
        layout.addWidget(self.log_text)
        
        self.tabs.addTab(tab, "1. SSH & Files")

    def init_pipeline_tab(self):
        scroll = QScrollArea()
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setSpacing(10)
        layout.setContentsMargins(10, 10, 10, 10)

        # Project info display (read-only)
        self.pipeline_project_label = QLabel("📁 Project: None")
        self.pipeline_project_label.setStyleSheet("color: #2980B9; font-weight: bold;")
        layout.addWidget(self.pipeline_project_label)

        # 输入设置
        io_group = QGroupBox("Input Data (Remote Paths)")
        io_form = QFormLayout(io_group)
        self.fq1_edit = QLineEdit()
        self.fq2_edit = QLineEdit()
        io_form.addRow("Reads R1:", self.fq1_edit)
        io_form.addRow("Reads R2:", self.fq2_edit)
        layout.addWidget(io_group)

        # 2. 运行模式设置 (Direct vs Slurm) - 新增功能
        exec_group = QGroupBox("Execution Mode")
        exec_form = QFormLayout(exec_group)
        self.run_mode_combo = QComboBox()
        self.run_mode_combo.addItems(["Direct Python", "Slurm (srun)"])
        exec_form.addRow("Run Mode:", self.run_mode_combo)
        layout.addWidget(exec_group)

        # 3. 高级参数切换按钮
        self.toggle_btn = QPushButton("Advanced Parameters ▼")
        self.toggle_btn.setCheckable(True)
        layout.addWidget(self.toggle_btn)

        # 4. 联动区域
        self.adv_widget = QWidget()
        adv_vbox = QVBoxLayout(self.adv_widget)
        adv_vbox.setContentsMargins(0, 0, 0, 0)
        self.ui_param_map = {}

        # 遍历所有键值对并展示
        for section, data in self.pipeline_cfg.items():
            if section == 'global' and 'threads' not in data: continue
            
            sec_group = QGroupBox(f"Setting: {section}")
            sec_form = QFormLayout(sec_group)
            
            def add_rows(item_data, prefix_keys=[]):
                for k, v in item_data.items():
                    if isinstance(v, dict):
                        add_rows(v, prefix_keys + [k])
                    else:
                        edit = QLineEdit(str(v))
                        label_text = ".".join(prefix_keys + [k])
                        sec_form.addRow(f"{label_text}:", edit)
                        map_key = tuple([section] + prefix_keys + [k])
                        self.ui_param_map[map_key] = edit

            if isinstance(data, dict):
                add_rows(data)
                adv_vbox.addWidget(sec_group)
        
        layout.addWidget(self.adv_widget)
        self.adv_widget.setVisible(False)
        self.toggle_btn.clicked.connect(lambda: self.adv_widget.setVisible(self.toggle_btn.isChecked()))

        # 5. 保存设置选项
        settings_group = QGroupBox("Settings")
        settings_form = QFormLayout(settings_group)
        self.save_pipeline_cb = QCheckBox("Save Settings")
        self.save_pipeline_cb.setChecked(True)
        self.pipeline_verbose_cb = QCheckBox("Verbose Logging")
        self.pipeline_verbose_cb.setChecked(False)
        self.pipeline_resume_cb = QCheckBox("Resume from checkpoint (if exists)")
        self.pipeline_resume_cb.setChecked(True)
        settings_form.addRow("", self.save_pipeline_cb)
        settings_form.addRow("", self.pipeline_verbose_cb)
        settings_form.addRow("", self.pipeline_resume_cb)
        layout.addWidget(settings_group)

        # 6. 启动与日志
        self.run_btn = QPushButton("🚀 RUN PIPELINE")
        self.run_btn.clicked.connect(lambda: self.on_launch(step_name="preprocess"))
        layout.addWidget(self.run_btn)

        self.pipe_log_text = QTextEdit()
        self.pipe_log_text.setReadOnly(True)
        self.pipe_log_text.setStyleSheet("background-color: #0c0c0c; color: #00FF41; font-family: 'Consolas';")
        self.pipe_log_text.setMinimumHeight(400)
        layout.addWidget(self.pipe_log_text, stretch=1)
        layout.addStretch()

        scroll.setWidget(content)
        scroll.setWidgetResizable(True)
        self.tabs.addTab(scroll, "2. Pre-processing")

    def init_pathogen_tab(self):
        scroll = QScrollArea()
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setSpacing(10)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # Project info display (read-only)
        self.pathogen_project_label = QLabel("📁 Project: None")
        self.pathogen_project_label.setStyleSheet("color: #2980B9; font-weight: bold;")
        layout.addWidget(self.pathogen_project_label)
        
        # Input Settings
        input_group = QGroupBox("Input Settings")
        input_form = QFormLayout(input_group)
        
        # TaxID file dropdown (editable for Custom)
        self.taxid_file_combo = QComboBox()
        self.taxid_file_combo.setEditable(True)
        self.taxid_file_combo.addItems(list(PRESET_TAXID_FILES.keys()))
        self.taxid_file_combo.lineEdit().setPlaceholderText("Enter custom TaxID file path...")
        self.taxid_file_combo.currentTextChanged.connect(self.on_taxid_file_changed)
        input_form.addRow("TaxID List File:", self.taxid_file_combo)
        
        # Virus Index dropdown (editable for Custom)
        self.virus_index_combo = QComboBox()
        self.virus_index_combo.setEditable(True)
        self.virus_index_combo.addItems(list(PRESET_VIRUS_INDEXES.keys()))
        self.virus_index_combo.lineEdit().setPlaceholderText("Enter custom virus index path...")
        self.virus_index_combo.currentTextChanged.connect(self.on_virus_index_changed)
        input_form.addRow("Virus Index Path:", self.virus_index_combo)
        layout.addWidget(input_group)
        
        # Execution Mode
        exec_group = QGroupBox("Execution Mode")
        exec_form = QFormLayout(exec_group)
        self.pathogen_run_mode_combo = QComboBox()
        self.pathogen_run_mode_combo.addItems(["Direct Python", "Slurm (srun)"])
        exec_form.addRow("Run Mode:", self.pathogen_run_mode_combo)
        layout.addWidget(exec_group)
        
        # Advanced Parameters Toggle
        self.pathogen_toggle_btn = QPushButton("▼ Advanced Parameters")
        self.pathogen_toggle_btn.setCheckable(True)
        layout.addWidget(self.pathogen_toggle_btn)
        
        # Advanced Parameters Area
        self.pathogen_adv_widget = QWidget()
        adv_vbox = QVBoxLayout(self.pathogen_adv_widget)
        adv_vbox.setContentsMargins(0, 0, 0, 0)
        
        # Filter Parameters (Diamond)
        filter_group = QGroupBox("Filter Parameters (Diamond)")
        filter_form = QFormLayout(filter_group)
        self.diamond_pident_edit = QLineEdit("80")
        self.diamond_qcov_edit = QLineEdit("50")
        self.diamond_evalue_edit = QLineEdit("1e-10")
        filter_form.addRow("Identity (%):", self.diamond_pident_edit)
        filter_form.addRow("Coverage (%):", self.diamond_qcov_edit)
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
        
        # Tool Paths
        tools_group = QGroupBox("Tool Paths")
        tools_form = QFormLayout(tools_group)
        
        self.bowtie2_bin_edit = QLineEdit("bowtie2")
        self.samtools_bin_edit = QLineEdit("samtools")
        self.bedtools_bin_edit = QLineEdit("bedtools")
        self.bedgraph_to_bigwig_edit = QLineEdit("/public21/home/sc90258/huangying/tools/bedGraphToBigWig")
        self.blastn_bin_edit = QLineEdit("blastn")
        
        tools_form.addRow("bowtie2:", self.bowtie2_bin_edit)
        tools_form.addRow("samtools:", self.samtools_bin_edit)
        tools_form.addRow("bedtools:", self.bedtools_bin_edit)
        tools_form.addRow("bedGraphToBigWig:", self.bedgraph_to_bigwig_edit)
        tools_form.addRow("blastn:", self.blastn_bin_edit)
        
        adv_vbox.addWidget(tools_group)
        
        layout.addWidget(self.pathogen_adv_widget)
        self.pathogen_adv_widget.setVisible(False)
        self.pathogen_toggle_btn.clicked.connect(
            lambda: self.pathogen_adv_widget.setVisible(self.pathogen_toggle_btn.isChecked())
        )
        
        # Settings
        settings_group = QGroupBox("Settings")
        settings_form = QFormLayout(settings_group)
        self.save_pathogen_cb = QCheckBox("Save Settings")
        self.save_pathogen_cb.setChecked(True)
        self.verbose_cb = QCheckBox("Verbose Logging")
        self.verbose_cb.setChecked(False)
        self.resume_cb = QCheckBox("Resume from checkpoint (if exists)")
        self.resume_cb.setChecked(True)
        settings_form.addRow("", self.save_pathogen_cb)
        settings_form.addRow("", self.verbose_cb)
        settings_form.addRow("", self.resume_cb)
        layout.addWidget(settings_group)
        
        # Run Button
        self.pathogen_run_btn = QPushButton("🚀 RUN DETECTION")
        self.pathogen_run_btn.clicked.connect(self.on_pathogen_run)
        layout.addWidget(self.pathogen_run_btn)
        
        # Log Output - 紧贴 RUN 按钮和底部按钮
        self.pathogen_log_text = QTextEdit()
        self.pathogen_log_text.setReadOnly(True)
        self.pathogen_log_text.setStyleSheet("background-color: #0c0c0c; color: #00FF41; font-family: 'Consolas';")
        layout.addWidget(self.pathogen_log_text, stretch=1)
        
        # Result Buttons - 最底部
        result_btn_layout = QHBoxLayout()
        self.view_results_btn = QPushButton("💾 View Results")
        self.view_results_btn.clicked.connect(self.on_view_pathogen_results)
        
        self.export_csv_btn = QPushButton("💾 Export CSV")
        self.export_csv_btn.clicked.connect(self.on_export_pathogen_csv)
        
        result_btn_layout.addWidget(self.view_results_btn)
        result_btn_layout.addWidget(self.export_csv_btn)
        result_btn_layout.addStretch()
        layout.addLayout(result_btn_layout)
        
        scroll.setWidget(content)
        scroll.setWidgetResizable(True)
        self.tabs.addTab(scroll, "3. Pathogen Detection")

    def on_list_remote(self):
        path = self.search_dir_edit.text().strip()
        self.run_ssh_task(command=f"ls -lh {path}", log_widget=self.log_text)

    def on_create_project(self):
        base = self.search_dir_edit.text().strip().rstrip('/')
        proj_name, ok = QInputDialog.getText(self, "Project Setup", "Enter New/Exist Project Name:", QLineEdit.Normal)
        if ok and proj_name:
            target_path = f"{base}/{proj_name}"
            self.run_ssh_task(target_dir=target_path, task_type="mkdir", log_widget=self.log_text)

    def on_task_finished(self, success, msg):
        task_type = getattr(self.worker, 'task_type_ext', 'default')
        
        if task_type == "pathogen":
            if success:
                self.pathogen_log_text.append("=" * 50)
                self.pathogen_log_text.append("✅ Pathogen detection completed!")
                self.pathogen_log_text.append(f"Results saved to: {self.pathogen_output_folder}")
                self.pathogen_log_text.append("=" * 50)
            else:
                self.pathogen_log_text.append(f"❌ Error: {msg}")
            return
        
        if success and msg == "CREATED":
            self.current_project_dir = self.worker.target_dir
            self.sync_info_label.setText(f"Active Project: {self.current_project_dir}")
            self.pipeline_project_label.setText(f"📁 Project: {self.current_project_dir}")
            self.pathogen_project_label.setText(f"📁 Project: {self.current_project_dir}")
            QMessageBox.information(self, "Success", f"Workspace ready:\n{self.current_project_dir}")
        elif msg == "EXISTS":
            QMessageBox.warning(self, "Conflict", "Directory already exists!")

    def on_launch(self, step_name="preprocess"):
        if not self.current_project_dir:
            QMessageBox.warning(self, "Path Error", "No project folder selected!")
            return
        
        # 确定远程子目录
        target_sub_dir = f"{self.current_project_dir}/{step_name}"
        sync_yaml = self.get_ui_as_yaml()
        
        # 组合原始运行命令
        if step_name == "preprocess":
            # 构建额外参数
            extra_args = ""
            if hasattr(self, 'pipeline_verbose_cb') and self.pipeline_verbose_cb.isChecked():
                extra_args += " --verbose"
            if hasattr(self, 'pipeline_resume_cb') and not self.pipeline_resume_cb.isChecked():
                extra_args += " --no-resume"
            
            base_cmd = f"mNGS_pip -1 {self.fq1_edit.text()} -2 {self.fq2_edit.text()} -o results{extra_args}"
        
        # 处理运行模式
        if self.run_mode_combo.currentText() == "Slurm (srun)":
            cmd = f"srun --job-name=mNGS_{step_name} {base_cmd}"
        else:
            cmd = base_cmd
        
        # 记录并运行，日志定向到当前页面的 pipe_log_text
        self.pipe_log_text.clear()
        self.run_ssh_task(command=cmd, upload_cfg=sync_yaml, target_dir=target_sub_dir, log_widget=self.pipe_log_text)

    def on_pathogen_run(self):
        if not self.current_project_dir:
            QMessageBox.warning(self, "Path Error", "No project folder selected!")
            return
        
        # Validate inputs - 使用新的下拉菜单获取路径
        taxid_file = self.get_taxid_file_path()
        virus_index = self.get_virus_index_path()
        
        if not taxid_file:
            QMessageBox.warning(self, "Input Error", "Please select or enter TaxID List File path!")
            return
        if not virus_index:
            QMessageBox.warning(self, "Input Error", "Please select or enter Virus Index Path!")
            return
        
        # Get output folder
        self.pathogen_output_folder = f"{self.current_project_dir}/{self.pathogen_output_folder_edit.text().strip()}"
        
        # Build command
        pident = self.diamond_pident_edit.text().strip()
        qcov = self.diamond_qcov_edit.text().strip()
        evalue = self.diamond_evalue_edit.text().strip()
        threads = "64"  # 默认使用所有 threads
        blastn_id = self.blastn_identity_edit.text().strip()
        blastn_evalue = self.blastn_evalue_edit.text().strip()
        blastn_cov = self.blastn_coverage_edit.text().strip()
        
        # Tab 2 的输出在 preprocess 子目录中
        tab2_step = "preprocess"
        diamond_results = f"{self.current_project_dir}/{tab2_step}/results.diamond.tsv"
        megahit_results = f"{self.current_project_dir}/{tab2_step}/results_megahit_result"
        
        # 从 Tab 2 获取原始 fastq 路径
        fq1 = self.fq1_edit.text()
        fq2 = self.fq2_edit.text()
        
        # 构建额外参数
        extra_args = ""
        if self.verbose_cb.isChecked():
            extra_args += " --verbose"
        if not self.resume_cb.isChecked():
            extra_args += " --no-resume"
        
        # 使用 pathogen_detect 命令（通过 pip install -e . 安装）
        base_cmd = f"pathogen_detect \
            --diamond {diamond_results} \
            --fq1 {fq1} --fq2 {fq2} \
            --megahit {megahit_results} \
            --taxids {taxid_file} \
            --virus-index {virus_index} \
            --output {self.pathogen_output_folder} \
            --bowtie2-bin {self.bowtie2_bin_edit.text()} \
            --samtools-bin {self.samtools_bin_edit.text()} \
            --bedtools-bin {self.bedtools_bin_edit.text()} \
            --bedgraph-to-bigwig-bin {self.bedgraph_to_bigwig_edit.text()} \
            --blastn-bin {self.blastn_bin_edit.text()} \
            --pident {pident} --qcov {qcov} --evalue {evalue} \
            --blastn-identity {blastn_id} --blastn-evalue {blastn_evalue} --blastn-coverage {blastn_cov} \
            --threads {threads}{extra_args}"
        
        # Handle run mode
        if self.pathogen_run_mode_combo.currentText() == "Slurm (srun)":
            cmd = f"srun --job-name=pathogen {base_cmd}"
        else:
            cmd = base_cmd
        
        # Run
        self.pathogen_log_text.clear()
        self.pathogen_log_text.append(f"Starting pathogen detection...")
        self.pathogen_log_text.append(f"Project: {self.current_project_dir}")
        self.pathogen_log_text.append(f"Output: {self.pathogen_output_folder}")
        self.pathogen_log_text.append("-" * 50)
        
        # Run pathogen detection (directory will be created by the script)
        self.run_ssh_task(
            command=cmd,
            log_widget=self.pathogen_log_text,
            task_type="pathogen",
            finished_callback=self.on_task_finished
        )

    def on_view_pathogen_results(self):
        if not self.current_project_dir:
            QMessageBox.warning(self, "Path Error", "No project folder selected!")
            return
        
        results_file = f"{self.current_project_dir}/{self.pathogen_output_folder_edit.text().strip()}/detection_results.tsv"
        
        # Try to read results from remote
        try:
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(self.host_edit.text().strip(), port=int(self.port_edit.text()),
                       username=self.user_edit.text().strip(), password=self.pwd_edit.text(), timeout=10)
            
            stdin, stdout, stderr = ssh.exec_command(f"cat {results_file}")
            content = stdout.read().decode('utf-8')
            ssh.close()
            
            # 检查是否为空或只有表头（支持新旧三种格式）
            old_header_v1 = "Virus_Name\tTaxID\tCoverage\tDepth\tReads_Count\tContigs_Count\tContig_IDs"
            old_header_v2 = "Virus_Name\tTaxID\tBest_Hit_Seq\tCoverage(%)\tDepth\tContig_ID"
            new_header = "Virus_Seq\tBest_Hit_Seq\tCoverage(%)\tDepth\tContig_ID"
            if not content or content.strip() == old_header_v1 or content.strip() == old_header_v2 or content.strip() == new_header:
                QMessageBox.information(self, "No Results", "No pathogen detection results found.")
                return
            
            # Parse results
            results = []
            lines = content.strip().split('\n')
            
            # Detect format from header
            header_line = lines[0] if lines else ""
            has_taxid = "TaxID" in header_line
            
            if len(lines) > 1:
                for line in lines[1:]:
                    fields = line.strip().split('\t')
                    if has_taxid:
                        # Old format: Virus_Name, TaxID, Best_Hit_Seq, Coverage, Depth, Contig_ID
                        if len(fields) >= 6:
                            results.append({
                                'Virus_Seq': fields[0],
                                'TaxID': fields[1],
                                'Best_Hit_Seq': fields[2] if len(fields) > 2 else '',
                                'Coverage': fields[3] if len(fields) > 3 else '',
                                'Depth': fields[4] if len(fields) > 4 else '',
                                'Contig_ID': fields[5] if len(fields) > 5 else ''
                            })
                    else:
                        # New format: Virus_Seq, Best_Hit_Seq, Coverage, Depth, Contig_ID
                        if len(fields) >= 5:
                            results.append({
                                'Virus_Seq': fields[0],
                                'Best_Hit_Seq': fields[1] if len(fields) > 1 else '',
                                'Coverage': fields[2] if len(fields) > 2 else '',
                                'Depth': fields[3] if len(fields) > 3 else '',
                                'Contig_ID': fields[4] if len(fields) > 4 else ''
                            })
                
                # Show dialog
                dialog = ResultsDialog(results, self.current_project_dir, has_taxid, self)
                dialog.exec_()
            else:
                QMessageBox.information(self, "No Results", "No pathogen detected.")
                
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to read results:\n{str(e)}")

    def get_remote_script_base(self):
        """获取远程服务器上的脚本基础路径"""
        import re
        
        try:
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(self.host_edit.text().strip(), port=int(self.port_edit.text()),
                       username=self.user_edit.text().strip(), password=self.pwd_edit.text(), timeout=10)
            
            # 使用 bash -l -c 并将 stderr 重定向到 stdout（get_script_base 输出到 stderr）
            stdin, stdout, stderr = ssh.exec_command("bash -l -c 'get_script_base' 2>&1")
            stdout_output = stdout.read().decode('utf-8', errors='ignore')
            ssh.close()
            
            # 清理 ANSI 颜色代码，提取包含 hvFinder 的路径
            clean_output = re.sub(r'\x1b\[[0-9;]*m', '', stdout_output)
            lines = [line.strip() for line in clean_output.strip().split('\n') if line.strip()]
            
            # 查找包含 hvFinder 的路径行
            for line in reversed(lines):
                if 'hvFinder' in line and line.startswith('/'):
                    return line
            
            return ''
        except Exception as e:
            return ''

    def on_load_default_taxids(self):
        """保留兼容性，现在用于刷新脚本基础路径"""
        self._script_base = self.get_remote_script_base()
        if not self._script_base:
            QMessageBox.warning(self, "Error", "Failed to get script base path\n\nPlease ensure hvFinder is installed on remote server with: pip install -e .")

    def on_load_default_virus_index(self):
        """保留兼容性，现在用于刷新脚本基础路径"""
        self._script_base = self.get_remote_script_base()
        if not self._script_base:
            QMessageBox.warning(self, "Error", "Failed to get script base path")

    def on_virus_index_changed(self, text):
        """Virus Index 下拉菜单选择变化"""
        # 更新内部存储的脚本基础路径
        if not hasattr(self, '_script_base') or not self._script_base:
            self._script_base = self.get_remote_script_base()
        
        # 如果选择 Custom，清空输入框内容（placeholder 会显示）
        if text == "Custom":
            self.virus_index_combo.lineEdit().clear()

    def on_taxid_file_changed(self, text):
        """TaxID File 下拉菜单选择变化"""
        # 更新内部存储的脚本基础路径
        if not hasattr(self, '_script_base') or not self._script_base:
            self._script_base = self.get_remote_script_base()
        
        # 如果选择 Custom，清空输入框内容（placeholder 会显示）
        if text == "Custom":
            self.taxid_file_combo.lineEdit().clear()

    def resolve_preset_path(self, relative_path):
        """将相对路径解析为绝对路径"""
        if not relative_path:
            return ""
        if not hasattr(self, '_script_base') or not self._script_base:
            self._script_base = self.get_remote_script_base()
        if self._script_base:
            return f"{self._script_base}/{relative_path}"
        return relative_path

    def get_virus_index_path(self):
        """获取当前选择的 Virus Index 路径"""
        current_text = self.virus_index_combo.currentText()
        # 如果是预设选项（非 Custom），返回解析后的路径
        if current_text in PRESET_VIRUS_INDEXES and current_text != "Custom":
            return self.resolve_preset_path(PRESET_VIRUS_INDEXES[current_text])
        else:
            # Custom 或用户输入的路径
            return current_text

    def get_taxid_file_path(self):
        """获取当前选择的 TaxID File 路径"""
        current_text = self.taxid_file_combo.currentText()
        # 如果是预设选项（非 Custom），返回解析后的路径
        if current_text in PRESET_TAXID_FILES and current_text != "Custom":
            return self.resolve_preset_path(PRESET_TAXID_FILES[current_text])
        else:
            return current_text

    def on_export_pathogen_csv(self):
        self.on_view_pathogen_results()

    def run_ssh_task(self, command=None, upload_cfg=None, target_dir=None, task_type="cmd", log_widget=None, finished_callback=None):
        self.save_settings()
        self.worker = SSHWorker(
            self.host_edit.text().strip(), int(self.port_edit.text()),
            self.user_edit.text().strip(), self.pwd_edit.text(),
            command=command, upload_cfg=upload_cfg, target_dir=target_dir, task_type=task_type
        )
        # 动态绑定日志信号到指定的组件
        if log_widget:
            self.worker.log_signal.connect(lambda msg: self.safe_log(msg, log_widget))
        else:
            self.worker.log_signal.connect(lambda msg: self.safe_log(msg, self.log_text))
            
        self.worker.status_signal.connect(self.update_status_ui)
        
        if finished_callback:
            self.worker.finished_signal.connect(finished_callback)
        else:
            self.worker.finished_signal.connect(self.on_task_finished)
        
        self.worker.task_type_ext = task_type
        self.worker.start()

    def get_ui_as_yaml(self):
        """将 UI 参数正确回填至嵌套字典并生成 YAML 字符串"""
        import copy
        new_cfg = copy.deepcopy(self.pipeline_cfg)
        
        for path_tuple, edit in self.ui_param_map.items():
            val = edit.text()
            # 自动类型转换
            if val.lower() == 'true': val = True
            elif val.lower() == 'false': val = False
            elif val.isdigit(): val = int(val)
            
            # 根据元组路径导航并更新嵌套字典
            temp = new_cfg
            for key in path_tuple[:-1]:
                temp = temp.setdefault(key, {})
            temp[path_tuple[-1]] = val
            
        return yaml.dump(new_cfg, sort_keys=False, allow_unicode=True)

    def update_status_ui(self, code):
        styles = {0: ("● Disconnected", "#7F8C8D"), 1: ("● Online", "#27AE60"), 2: ("● Running Task", "#E67E22")}
        text, color = styles.get(code, styles[0])
        self.status_label.setText(text)
        self.status_label.setStyleSheet(f"color: {color}; font-weight: bold;")

    def safe_log(self, text, widget):
        widget.append(text)
        widget.moveCursor(QTextCursor.End)

    def save_settings(self):
        config = configparser.ConfigParser()
        config['SSH'] = {
            'host': self.host_edit.text(), 'port': self.port_edit.text(),
            'user': self.user_edit.text(), 'pwd': self.pwd_edit.text(),
            'base': self.search_dir_edit.text()
        }
        # 保存 Pipeline 设置（如果勾选了保存选项）
        if hasattr(self, 'save_pipeline_cb') and self.save_pipeline_cb.isChecked():
            config['PIPELINE'] = {
                'fq1': self.fq1_edit.text(),
                'fq2': self.fq2_edit.text(),
                'run_mode': self.run_mode_combo.currentText(),
                'verbose': str(self.pipeline_verbose_cb.isChecked()),
                'resume': str(self.pipeline_resume_cb.isChecked())
            }
            # 保存高级参数
            config['ADVANCED'] = {}
            for path_tuple, edit in self.ui_param_map.items():
                key = ".".join(path_tuple)
                config['ADVANCED'][key] = edit.text()
        
        # 保存 Pathogen Detection 设置
        if hasattr(self, 'save_pathogen_cb') and self.save_pathogen_cb.isChecked():
            config['PATHOGEN'] = {
                'taxid_file_preset': self.taxid_file_combo.currentText(),
                'virus_index_preset': self.virus_index_combo.currentText(),
                'run_mode': self.pathogen_run_mode_combo.currentText(),
                'output_folder': self.pathogen_output_folder_edit.text(),
                'diamond_pident': self.diamond_pident_edit.text(),
                'diamond_qcov': self.diamond_qcov_edit.text(),
                'diamond_evalue': self.diamond_evalue_edit.text(),
                'blastn_identity': self.blastn_identity_edit.text(),
                'blastn_evalue': self.blastn_evalue_edit.text(),
                'blastn_coverage': self.blastn_coverage_edit.text(),
                'bowtie2_bin': self.bowtie2_bin_edit.text(),
                'samtools_bin': self.samtools_bin_edit.text(),
                'bedtools_bin': self.bedtools_bin_edit.text(),
                'bedgraph_to_bigwig_bin': self.bedgraph_to_bigwig_edit.text(),
                'blastn_bin': self.blastn_bin_edit.text(),
                'verbose': str(self.verbose_cb.isChecked()),
                'resume': str(self.resume_cb.isChecked())
            }
        
        with open(self.settings_file, 'w') as f: config.write(f)

    def load_saved_settings(self):
        if os.path.exists(self.settings_file):
            config = configparser.ConfigParser(); config.read(self.settings_file)
            if 'SSH' in config:
                self.host_edit.setText(config['SSH'].get('host', ''))
                self.port_edit.setText(config['SSH'].get('port', '22'))
                self.user_edit.setText(config['SSH'].get('user', ''))
                self.pwd_edit.setText(config['SSH'].get('pwd', ''))
                self.search_dir_edit.setText(config['SSH'].get('base', ''))
            
            # 加载 Pipeline 设置
            if 'PIPELINE' in config:
                self.fq1_edit.setText(config['PIPELINE'].get('fq1', ''))
                self.fq2_edit.setText(config['PIPELINE'].get('fq2', ''))
                run_mode = config['PIPELINE'].get('run_mode', 'Direct Python')
                idx = self.run_mode_combo.findText(run_mode)
                if idx >= 0:
                    self.run_mode_combo.setCurrentIndex(idx)
                
                # 加载复选框状态
                if hasattr(self, 'pipeline_verbose_cb') and 'verbose' in config['PIPELINE']:
                    self.pipeline_verbose_cb.setChecked(config['PIPELINE']['verbose'] == 'True')
                if hasattr(self, 'pipeline_resume_cb') and 'resume' in config['PIPELINE']:
                    self.pipeline_resume_cb.setChecked(config['PIPELINE']['resume'] == 'True')
            
            # 加载高级参数
            if 'ADVANCED' in config:
                for key, value in config['ADVANCED'].items():
                    path_tuple = tuple(key.split("."))
                    if path_tuple in self.ui_param_map:
                        self.ui_param_map[path_tuple].setText(value)
            
            # 加载 Pathogen Detection 设置
            if 'PATHOGEN' in config:
                # 加载 TaxID 下拉菜单
                taxid_preset = config['PATHOGEN'].get('taxid_file_preset', '')
                if taxid_preset:
                    idx = self.taxid_file_combo.findText(taxid_preset)
                    if idx >= 0:
                        self.taxid_file_combo.setCurrentIndex(idx)
                    else:
                        self.taxid_file_combo.setCurrentText(taxid_preset)
                else:
                    # 兼容旧配置：如果有 taxid_file，设置为 Custom
                    old_taxid = config['PATHOGEN'].get('taxid_file', '')
                    if old_taxid:
                        self.taxid_file_combo.setCurrentText(old_taxid)
                
                # 加载 Virus Index 下拉菜单
                virus_preset = config['PATHOGEN'].get('virus_index_preset', '')
                if virus_preset:
                    idx = self.virus_index_combo.findText(virus_preset)
                    if idx >= 0:
                        self.virus_index_combo.setCurrentIndex(idx)
                    else:
                        self.virus_index_combo.setCurrentText(virus_preset)
                else:
                    # 兼容旧配置：如果有 virus_index，设置为 Custom
                    old_virus = config['PATHOGEN'].get('virus_index', '')
                    if old_virus:
                        self.virus_index_combo.setCurrentText(old_virus)
                
                self.pathogen_output_folder_edit.setText(config['PATHOGEN'].get('output_folder', 'pathogen_results'))
                self.diamond_pident_edit.setText(config['PATHOGEN'].get('diamond_pident', '80'))
                self.diamond_qcov_edit.setText(config['PATHOGEN'].get('diamond_qcov', '50'))
                self.diamond_evalue_edit.setText(config['PATHOGEN'].get('diamond_evalue', '1e-10'))
                self.blastn_identity_edit.setText(config['PATHOGEN'].get('blastn_identity', '70'))
                self.blastn_evalue_edit.setText(config['PATHOGEN'].get('blastn_evalue', '1e-4'))
                self.blastn_coverage_edit.setText(config['PATHOGEN'].get('blastn_coverage', '70'))
                self.bowtie2_bin_edit.setText(config['PATHOGEN'].get('bowtie2_bin', 'bowtie2'))
                self.samtools_bin_edit.setText(config['PATHOGEN'].get('samtools_bin', 'samtools'))
                self.bedtools_bin_edit.setText(config['PATHOGEN'].get('bedtools_bin', 'bedtools'))
                self.bedgraph_to_bigwig_edit.setText(config['PATHOGEN'].get('bedgraph_to_bigwig_bin', '/public21/home/sc90258/huangying/tools/bedGraphToBigWig'))
                self.blastn_bin_edit.setText(config['PATHOGEN'].get('blastn_bin', 'blastn'))
                
                run_mode = config['PATHOGEN'].get('run_mode', 'Direct Python')
                idx = self.pathogen_run_mode_combo.findText(run_mode)
                if idx >= 0:
                    self.pathogen_run_mode_combo.setCurrentIndex(idx)
                
                # 加载复选框状态
                if 'verbose' in config['PATHOGEN']:
                    self.verbose_cb.setChecked(config['PATHOGEN']['verbose'] == 'True')
                if 'resume' in config['PATHOGEN']:
                    self.resume_cb.setChecked(config['PATHOGEN']['resume'] == 'True')

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())