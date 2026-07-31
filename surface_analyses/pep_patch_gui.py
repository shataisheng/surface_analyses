#!/usr/bin/env python3
"""
PEP-Patch 批量分析 GUI
======================
支持批量处理 PDB 文件的静电势和疏水性分析。
双击 .py 文件或命令行运行: pep_patch_gui
"""

import sys, os, subprocess, threading, csv, queue, glob as globmod, io
from datetime import datetime
from collections import OrderedDict
from enum import Enum
from dataclasses import dataclass, field
from contextlib import redirect_stdout, redirect_stderr

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

# ─── 平台适配 ───────────────────────────────────────────
try:
    from surface_analyses.platform_config import get_config
    CFG = get_config()
except ImportError:
    CFG = None

# ─── 预设 ───────────────────────────────────────────────
SCALE_PRESETS = OrderedDict([
    ("Wimley-White (残基级)", "test/trastuzumab/wimley-white-scaled.csv"),
    ("Eisenberg (原子类型)", "eisenberg"),
    ("Crippen (原子 logP)", "crippen"),
])
SURFACE_TYPES = ["sas", "ses", "gauss"]
SURFTYPES = ["normal", "sc_norm", "atom_norm"]


class JobStatus(Enum):
    QUEUED = "⬜"
    RUNNING_ES = "⚡"
    RUNNING_HB = "💧"
    DONE = "✅"
    FAILED = "❌"
    SKIPPED = "⏭"


@dataclass
class JobItem:
    pdb_path: str
    stem: str = ""
    run_es: bool = True
    run_hb: bool = True
    status: JobStatus = JobStatus.QUEUED
    error: str = ""
    es_outputs: dict = field(default_factory=dict)
    hb_outputs: dict = field(default_factory=dict)

    def __post_init__(self):
        if not self.stem:
            self.stem = os.path.splitext(os.path.basename(self.pdb_path))[0]
            for sfx in ["_fixed", "_out", "_output"]:
                if self.stem.endswith(sfx):
                    self.stem = self.stem[:-len(sfx)]
                    break


# ─── 主窗口 ─────────────────────────────────────────────
class PepPatchGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("PEP-Patch — Batch Surface Analysis")
        self.root.geometry("1150x820")
        self.root.minsize(950, 650)

        self.root.configure(bg=self.COLORS["bg"])
        self._setup_styles()

        # ── 变量 ──
        self.jobs: list[JobItem] = []
        self.job_index = -1
        self.running = False
        self.stop_requested = False
        self.process = None
        self.msg_queue = queue.Queue()

        self.run_es = tk.BooleanVar(value=True)
        self.run_hb = tk.BooleanVar(value=True)

        # ES params
        self.es_apbs_dir = tk.StringVar(value="")
        self.es_ph = tk.StringVar(value="5.5")
        self.es_ph_enabled = tk.BooleanVar(value=True)
        self.es_surface_type = tk.StringVar(value="sas")
        self.es_probe_radius = tk.StringVar(value="0.14")
        self.es_patch_cutoff_pos = tk.StringVar(value="2.0")
        self.es_patch_cutoff_neg = tk.StringVar(value="-2.0")
        self.es_check_cdr = tk.BooleanVar(value=False)

        # HB params
        self.hb_scale = tk.StringVar(value="test/trastuzumab/wimley-white-scaled.csv")
        self.hb_surftype = tk.StringVar(value="normal")
        self.hb_compute_potential = tk.BooleanVar(value=True)
        self.hb_compute_patches = tk.BooleanVar(value=True)
        self.hb_compute_sap = tk.BooleanVar(value=False)
        self.hb_compute_sh = tk.BooleanVar(value=False)
        self.hb_blur_rad = tk.StringVar(value="0.5")
        self.hb_sh_rad = tk.StringVar(value="0.8")
        self.hb_solv_rad = tk.StringVar(value="0.14")
        self.hb_grid_spacing = tk.StringVar(value="0.05")
        self.hb_rcut = tk.StringVar(value="0.5")
        self.hb_alpha = tk.StringVar(value="15.0")
        self.hb_patch_min = tk.StringVar(value="0.12")

        # Output
        self.gen_ply = tk.BooleanVar(value=True)
        self.gen_detailed_csv = tk.BooleanVar(value=True)

        self._build_ui()
        self._auto_detect_apbs()
        self._poll_queue()

    # ───────────────────────────────────────────────────────
    #  现代主题
    # ───────────────────────────────────────────────────────
    COLORS = {
        "bg":           "#f4f6fb",
        "panel":        "#ffffff",
        "panel_alt":    "#f1f5f9",
        "border":       "#e2e8f0",
        "text":         "#1e293b",
        "muted":        "#64748b",
        "accent":       "#4f46e5",
        "accent_hover": "#4338ca",
        "accent_soft":  "#eef2ff",
        "success":      "#16a34a",
        "warning":      "#f59e0b",
        "danger":       "#ef4444",
        "header_bg":    "#1e293b",
        "console_bg":   "#0f172a",
        "console_fg":   "#e2e8f0",
    }

    def _setup_styles(self):
        style = ttk.Style()
        style.theme_use("clam")
        c = self.COLORS

        style.configure(".", background=c["bg"], foreground=c["text"],
                        font=("Segoe UI", 10), borderwidth=0, relief="flat")
        style.configure("TFrame", background=c["bg"])
        style.configure("TLabel", background=c["bg"], foreground=c["text"],
                        font=("Segoe UI", 10))
        style.configure("Title.TLabel", font=("Segoe UI", 16, "bold"),
                        foreground=c["text"], background=c["bg"])
        style.configure("Subtitle.TLabel", font=("Segoe UI", 9),
                        foreground=c["muted"], background=c["bg"])
        style.configure("Section.TLabel", font=("Segoe UI", 11, "bold"),
                        foreground=c["accent"], background=c["panel"])

        # 卡片式 LabelFrame
        style.configure("Card.TLabelframe", background=c["panel"],
                        borderwidth=1, relief="solid", bordercolor=c["border"], padding=10)
        style.configure("Card.TLabelframe.Label", background=c["panel"],
                        foreground=c["accent"], font=("Segoe UI", 10, "bold"))

        # 按钮
        style.configure("TButton", font=("Segoe UI", 10), padding=(10, 6),
                        background=c["panel"], foreground=c["text"],
                        borderwidth=1, relief="solid", bordercolor=c["border"])
        style.map("TButton", background=[("active", c["accent_soft"])],
                  bordercolor=[("active", c["accent"])])

        style.configure("Accent.TButton", font=("Segoe UI", 11, "bold"),
                        padding=(16, 9), background=c["accent"],
                        foreground="white", borderwidth=0, relief="flat")
        style.map("Accent.TButton",
                  background=[("active", c["accent_hover"]), ("disabled", "#cbd5e1")],
                  foreground=[("disabled", "#f1f5f9")])

        style.configure("Ghost.TButton", font=("Segoe UI", 10), padding=(10, 6),
                        background=c["panel"], foreground=c["muted"],
                        borderwidth=1, relief="solid", bordercolor=c["border"])
        style.map("Ghost.TButton", background=[("active", c["panel_alt"])])

        # 进度条
        style.configure("TProgressbar", thickness=8, background=c["accent"],
                        troughcolor=c["border"], borderwidth=0, relief="flat")

        # 笔记本 / 标签页
        style.configure("TNotebook", background=c["bg"], borderwidth=0)
        style.configure("TNotebook.Tab", font=("Segoe UI", 10),
                        background=c["panel_alt"], foreground=c["muted"],
                        padding=(18, 8), borderwidth=0)
        style.map("TNotebook.Tab",
                  background=[("selected", c["panel"])],
                  foreground=[("selected", c["accent"])])

        # 树形视图
        style.configure("Treeview", background=c["panel"], foreground=c["text"],
                        font=("Segoe UI", 9), rowheight=26,
                        fieldbackground=c["panel"], borderwidth=0, relief="flat")
        style.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"),
                        background=c["accent"], foreground="white",
                        relief="flat", borderwidth=0)
        style.map("Treeview.Heading", background=[("active", c["accent_hover"])])
        style.map("Treeview", background=[("selected", c["accent_soft"])],
                  foreground=[("selected", c["accent"])])

        # 输入控件
        style.configure("TEntry", font=("Segoe UI", 10), padding=5,
                        fieldbackground=c["panel"], foreground=c["text"],
                        borderwidth=1, relief="solid", bordercolor=c["border"])
        style.configure("TCombobox", font=("Segoe UI", 10), padding=5,
                        fieldbackground=c["panel"], foreground=c["text"],
                        borderwidth=1, relief="solid", bordercolor=c["border"])
        style.map("TCombobox", fieldbackground=[("readonly", c["panel"])])

        # 分隔线 / 面板
        style.configure("TSeparator", background=c["border"])
        style.configure("TPanedwindow", background=c["bg"], borderwidth=0)

    # ═══════════════════════════════════════════════════════
    #  UI 构建
    # ═══════════════════════════════════════════════════════
    def _build_ui(self):
        c = self.COLORS
        # 顶部标题栏
        header = tk.Frame(self.root, bg=c["header_bg"], height=62)
        header.pack(fill=tk.X)
        header.pack_propagate(False)
        tk.Label(header, text="⚡ PEP-Patch", bg=c["header_bg"], fg="white",
                 font=("Segoe UI", 16, "bold")).pack(side=tk.LEFT, padx=18, pady=10)
        tk.Label(header, text="批量表面分析 · 静电势 / 疏水性",
                 bg=c["header_bg"], fg="#cbd5e1",
                 font=("Segoe UI", 9)).pack(side=tk.LEFT, padx=(2, 0), pady=(16, 0))

        main_pw = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_pw.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 左侧控制面板
        left = ttk.Frame(main_pw, width=420)
        main_pw.add(left, weight=1)
        canvas = tk.Canvas(left, highlightthickness=0, bg=c["bg"])
        sb = ttk.Scrollbar(left, orient=tk.VERTICAL, command=canvas.yview)
        self.scroll_frame = ttk.Frame(canvas)
        self.scroll_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=sb.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(int(-1*(e.delta/120)), "units"))

        self._build_file_list()
        self._build_mode()
        self._build_es_params()
        self._build_hb_params()
        self._build_output()
        self._build_run()

        # 右侧面板
        right = ttk.Frame(main_pw)
        main_pw.add(right, weight=2)
        self._build_right(right)

    def _build_file_list(self):
        f = ttk.LabelFrame(self.scroll_frame, text="📁 待分析文件 (Batch)", style="Card.TLabelframe")
        f.pack(fill=tk.X, padx=4, pady=4)
        tb = ttk.Frame(f)
        tb.pack(fill=tk.X, pady=(0, 4))
        ttk.Button(tb, text="+ 添加文件", command=self._add_files, width=10, style="Ghost.TButton").pack(side=tk.LEFT, padx=2)
        ttk.Button(tb, text="+ 添加目录", command=self._add_dir, width=10, style="Ghost.TButton").pack(side=tk.LEFT, padx=2)
        ttk.Button(tb, text="✕ 清空", command=self._clear_jobs, width=6, style="Ghost.TButton").pack(side=tk.LEFT, padx=2)
        self.job_count_lbl = ttk.Label(tb, text="0 files", foreground=self.COLORS["muted"])
        self.job_count_lbl.pack(side=tk.RIGHT, padx=4)

        lf = ttk.Frame(f)
        lf.pack(fill=tk.BOTH, expand=True)
        self.job_tree = ttk.Treeview(lf, columns=("status", "es", "hb", "name"), show="headings", height=8)
        for col, w, txt in [("status", 35, ""), ("es", 30, "ES"), ("hb", 30, "HB"), ("name", 250, "PDB File")]:
            self.job_tree.heading(col, text=txt, anchor="center" if col != "name" else "w")
            self.job_tree.column(col, width=w, anchor="center" if col != "name" else "w", stretch=(col == "name"))
        self.job_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        ts = ttk.Scrollbar(lf, orient=tk.VERTICAL, command=self.job_tree.yview)
        ts.pack(side=tk.RIGHT, fill=tk.Y)
        self.job_tree.configure(yscrollcommand=ts.set)
        c = self.COLORS
        for st, col in [("DONE", c["success"]), ("FAILED", c["danger"]),
                        ("RUNNING_ES", c["accent"]), ("RUNNING_HB", c["accent"]),
                        ("SKIPPED", c["warning"]), ("QUEUED", c["muted"])]:
            self.job_tree.tag_configure(f"st_{st}", foreground=col)
        self.job_tree.bind("<Delete>", lambda e: self._remove_selected())
        self.job_tree.bind("<Double-1>", lambda e: self._open_job_pdb())

    def _build_mode(self):
        f = ttk.LabelFrame(self.scroll_frame, text="🔬 分析类型", style="Card.TLabelframe")
        f.pack(fill=tk.X, padx=4, pady=4)
        ttk.Checkbutton(f, text="⚡ 静电势 (Electrostatic)", variable=self.run_es).grid(row=0, column=0, sticky="w", padx=4)
        ttk.Checkbutton(f, text="💧 疏水性 (Hydrophobic)", variable=self.run_hb).grid(row=0, column=1, sticky="w", padx=12)
        ttk.Label(f, text="对每个 PDB 执行选中的分析", font=("", 8, "italic"), foreground=self.COLORS["muted"]).grid(row=1, column=0, columnspan=2, sticky="w", padx=4)
        ttk.Separator(f, orient=tk.HORIZONTAL).grid(row=2, column=0, columnspan=2, sticky="ew", pady=4)
        ttk.Label(f, text="表面类型:").grid(row=3, column=0, sticky="w", padx=4)
        ttk.Combobox(f, textvariable=self.es_surface_type, values=SURFACE_TYPES, state="readonly", width=8).grid(row=3, column=0, sticky="e", padx=4)
        ttk.Label(f, text="探针半径 (nm):").grid(row=3, column=1, sticky="w", padx=4)
        ttk.Entry(f, textvariable=self.es_probe_radius, width=7).grid(row=3, column=1, sticky="e", padx=4)

    def _build_es_params(self):
        f = ttk.LabelFrame(self.scroll_frame, text="⚡ 静电势参数", style="Card.TLabelframe")
        f.pack(fill=tk.X, padx=4, pady=4)
        ttk.Label(f, text="APBS 目录:").grid(row=0, column=0, sticky="w")
        f0 = ttk.Frame(f)
        f0.grid(row=1, column=0, columnspan=2, sticky="ew", pady=2)
        ttk.Entry(f0, textvariable=self.es_apbs_dir, width=30).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(f0, text="浏览", command=self._browse_apbs, width=5).pack(side=tk.LEFT, padx=2)
        ttk.Label(f, text="留空 = 自动检测", font=("", 7, "italic")).grid(row=2, column=0, columnspan=2, sticky="w")
        fph = ttk.Frame(f)
        fph.grid(row=3, column=0, columnspan=2, sticky="ew", pady=4)
        ttk.Checkbutton(fph, text="pH:", variable=self.es_ph_enabled).pack(side=tk.LEFT)
        ttk.Entry(fph, textvariable=self.es_ph, width=5).pack(side=tk.LEFT, padx=4)
        ttk.Label(f, text="正阈值:").grid(row=4, column=0, sticky="w")
        ttk.Entry(f, textvariable=self.es_patch_cutoff_pos, width=7).grid(row=4, column=0, sticky="e", padx=4)
        ttk.Label(f, text="负阈值:").grid(row=4, column=1, sticky="w")
        ttk.Entry(f, textvariable=self.es_patch_cutoff_neg, width=7).grid(row=4, column=1, sticky="e", padx=4)
        ttk.Checkbutton(f, text="检测 CDR (抗体)", variable=self.es_check_cdr).grid(row=5, column=0, columnspan=2, sticky="w")

    def _build_hb_params(self):
        f = ttk.LabelFrame(self.scroll_frame, text="💧 疏水性参数", style="Card.TLabelframe")
        f.pack(fill=tk.X, padx=4, pady=4)
        ttk.Label(f, text="疏水表量表:").grid(row=0, column=0, sticky="w")
        fs = ttk.Frame(f)
        fs.grid(row=1, column=0, columnspan=2, sticky="ew", pady=2)
        ttk.Entry(fs, textvariable=self.hb_scale, width=24).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(fs, text="浏览", command=self._browse_scale, width=5).pack(side=tk.LEFT, padx=2)
        cb = ttk.Combobox(f, values=list(SCALE_PRESETS.keys()), state="readonly", width=28)
        cb.grid(row=2, column=0, columnspan=2, sticky="ew")
        cb.bind("<<ComboboxSelected>>", lambda e: self.hb_scale.set(SCALE_PRESETS[cb.get()]))
        ttk.Label(f, text="归一化:").grid(row=3, column=0, sticky="w", pady=(6, 0))
        ttk.Combobox(f, textvariable=self.hb_surftype, values=SURFTYPES, state="readonly", width=12).grid(row=3, column=1, sticky="w", padx=4, pady=(6, 0))
        of = ttk.Frame(f)
        of.grid(row=4, column=0, columnspan=2, sticky="ew", pady=6)
        for r, (txt, var) in enumerate([("疏水势", self.hb_compute_potential), ("Patch", self.hb_compute_patches),
                                         ("SAP", self.hb_compute_sap), ("SH", self.hb_compute_sh)]):
            ttk.Checkbutton(of, text=txt, variable=var).grid(row=r//2, column=r%2, sticky="w", padx=(0 if r%2==0 else 8))

    def _build_output(self):
        f = ttk.LabelFrame(self.scroll_frame, text="📤 输出", style="Card.TLabelframe")
        f.pack(fill=tk.X, padx=4, pady=4)
        ttk.Checkbutton(f, text="生成 PLY 可视化文件", variable=self.gen_ply).pack(anchor="w")
        ttk.Checkbutton(f, text="生成残基级详细 CSV", variable=self.gen_detailed_csv).pack(anchor="w")

    def _build_run(self):
        f = ttk.Frame(self.scroll_frame)
        f.pack(fill=tk.X, padx=4, pady=6)
        self.run_btn = ttk.Button(f, text="▶ 开始批量分析", style="Accent.TButton", command=self._start_batch)
        self.run_btn.pack(side=tk.LEFT, padx=2)
        self.stop_btn = ttk.Button(f, text="⏹ 停止", command=self._stop_batch, state=tk.DISABLED, style="Ghost.TButton")
        self.stop_btn.pack(side=tk.LEFT, padx=2)
        self.batch_progress = ttk.Progressbar(f, mode="determinate")
        self.batch_progress.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=6)
        self.batch_label = ttk.Label(f, text="")
        self.batch_label.pack(side=tk.RIGHT, padx=4)

    def _build_right(self, parent):
        c = self.COLORS
        nb = ttk.Notebook(parent)
        nb.pack(fill=tk.BOTH, expand=True)

        # Console
        cf = ttk.Frame(nb)
        nb.add(cf, text="📋 日志")
        self.console = scrolledtext.ScrolledText(cf, wrap=tk.WORD, font=("Consolas", 10),
            bg=c["console_bg"], fg=c["console_fg"], insertbackground="white",
            relief="flat", borderwidth=0)
        self.console.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        for tag, color in [("info", "#60a5fa"), ("success", "#4ade80"), ("warning", "#fbbf24"),
                           ("error", "#f87171"), ("header", "#a78bfa")]:
            self.console.tag_config(tag, foreground=color)
        self.console.tag_config("header", font=("Consolas", 10, "bold"))

        # Files
        ff = ttk.Frame(nb)
        nb.add(ff, text="📁 输出文件")
        self.files_list = tk.Listbox(ff, font=("Consolas", 10), bg=c["console_bg"], fg=c["console_fg"],
                                     selectbackground=c["accent"], relief="flat", borderwidth=0,
                                     highlightthickness=0)
        self.files_list.pack(fill=tk.BOTH, expand=True, padx=3, pady=3)
        self.files_list.bind("<Double-Button-1>", self._open_file)

    # ═══════════════════════════════════════════════════════
    #  文件管理
    # ═══════════════════════════════════════════════════════
    def _add_files(self):
        paths = filedialog.askopenfilenames(filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")])
        existing = {j.pdb_path for j in self.jobs}
        for p in paths:
            if p not in existing:
                self.jobs.append(JobItem(pdb_path=p, run_es=self.run_es.get(), run_hb=self.run_hb.get()))
                existing.add(p)
        self._refresh_jobs()

    def _add_dir(self):
        d = filedialog.askdirectory(title="选择包含 PDB 文件的目录")
        if not d:
            return
        existing = {j.pdb_path for j in self.jobs}
        added = 0
        for f in os.listdir(d):
            if f.endswith(".pdb"):
                p = os.path.join(d, f)
                if p not in existing:
                    self.jobs.append(JobItem(pdb_path=p, run_es=self.run_es.get(), run_hb=self.run_hb.get()))
                    existing.add(p)
                    added += 1
        self._refresh_jobs()
        self._log(f"从目录添加了 {added} 个 PDB 文件", "info")

    def _clear_jobs(self):
        if self.running:
            return
        self.jobs.clear()
        self._refresh_jobs()

    def _remove_selected(self):
        if self.running:
            return
        sel = self.job_tree.selection()
        indices = sorted([int(self.job_tree.item(s, "tags")[0]) for s in sel if self.job_tree.item(s, "tags")], reverse=True)
        for i in indices:
            if 0 <= i < len(self.jobs):
                self.jobs.pop(i)
        self._refresh_jobs()

    def _refresh_jobs(self):
        for item in self.job_tree.get_children():
            self.job_tree.delete(item)
        for i, job in enumerate(self.jobs):
            es_mark = "✓" if job.run_es else "—"
            hb_mark = "✓" if job.run_hb else "—"
            self.job_tree.insert("", tk.END,
                                 values=(job.status.value, es_mark, hb_mark, os.path.basename(job.pdb_path)),
                                 tags=(str(i), f"st_{job.status.name}"))
        self.job_count_lbl.config(text=f"{len(self.jobs)} files")
        self.batch_progress["maximum"] = max(len(self.jobs), 1)
        self.batch_progress["value"] = 0

    @staticmethod
    def _open_path(path: str) -> None:
        """跨平台打开文件/目录（Windows 用 startfile，macOS 用 open，Linux 用 xdg-open）。"""
        if sys.platform == "win32":
            os.startfile(path)
        elif sys.platform == "darwin":
            subprocess.run(["open", path], check=False)
        else:
            subprocess.run(["xdg-open", path], check=False)

    def _open_job_pdb(self):
        sel = self.job_tree.selection()
        if sel and self.job_tree.item(sel[0], "tags"):
            i = int(self.job_tree.item(sel[0], "tags")[0])
            if 0 <= i < len(self.jobs):
                self._open_path(self.jobs[i].pdb_path)

    def _auto_detect_apbs(self):
        if CFG:
            self.es_apbs_dir.set(CFG.default_apbs_work_dir)

    def _browse_apbs(self):
        d = filedialog.askdirectory(title="APBS 工作目录")
        if d:
            self.es_apbs_dir.set(d)

    def _browse_scale(self):
        p = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if p:
            self.hb_scale.set(p)

    # ═══════════════════════════════════════════════════════
    #  日志（线程安全）
    # ═══════════════════════════════════════════════════════
    def _log(self, msg, tag=None):
        self.msg_queue.put((msg, tag))

    def _log_now(self, msg, tag=None):
        self.console.insert(tk.END, msg + "\n", tag or "")
        self.console.see(tk.END)

    def _poll_queue(self):
        while True:
            try:
                msg, tag = self.msg_queue.get_nowait()
                self._log_now(msg, tag)
            except queue.Empty:
                break
        self.root.after(100, self._poll_queue)

    # ═══════════════════════════════════════════════════════
    #  批量处理引擎
    # ═══════════════════════════════════════════════════════
    def _start_batch(self):
        if self.running:
            return
        if not self.jobs:
            messagebox.showwarning("无文件", "请先添加 PDB 文件")
            return
        if not self.run_es.get() and not self.run_hb.get():
            messagebox.showwarning("无分析类型", "请至少选择一种分析类型")
            return

        for job in self.jobs:
            job.run_es = self.run_es.get()
            job.run_hb = self.run_hb.get()
            job.status = JobStatus.QUEUED
            job.error = ""
            job.es_outputs.clear()
            job.hb_outputs.clear()

        self.running = True
        self.stop_requested = False
        self.job_index = -1
        self.run_btn.config(state=tk.DISABLED)
        self.stop_btn.config(state=tk.NORMAL)
        # Create timestamped results folder
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.results_dir = os.path.join(str(CFG.project_root) if CFG else os.getcwd(), "results", ts)
        os.makedirs(self.results_dir, exist_ok=True)
        self._refresh_jobs()
        self._log(f"\n{'═'*50}", "header")
        self._log(f"  批量分析: {len(self.jobs)} 文件 | ES={self.run_es.get()} HB={self.run_hb.get()}", "header")
        self._log(f"{'═'*50}\n", "header")
        threading.Thread(target=self._batch_worker, daemon=True).start()

    def _stop_batch(self):
        self.stop_requested = True
        if self.process:
            try:
                self.process.terminate()
            except Exception:
                pass
        self._log("\n⏹ 用户停止", "warning")

    def _update_job_ui(self, job: JobItem):
        def _do():
            for item in self.job_tree.get_children():
                tags = self.job_tree.item(item, "tags")
                if tags and int(tags[0]) == self.jobs.index(job):
                    es, hb = ("✓" if job.run_es else "—"), ("✓" if job.run_hb else "—")
                    self.job_tree.item(item,
                                       values=(job.status.value, es, hb, os.path.basename(job.pdb_path)),
                                       tags=(str(self.jobs.index(job)), f"st_{job.status.name}"))
                    self.job_tree.see(item)
                    break
            self.batch_progress["value"] = self.job_index + 1
            self.batch_label.config(text=f"{self.job_index + 1}/{len(self.jobs)}")
        self.root.after(0, _do)

    def _batch_worker(self):
        for i, job in enumerate(self.jobs):
            if self.stop_requested:
                job.status = JobStatus.SKIPPED
                self._update_job_ui(job)
                continue

            self.job_index = i
            stem = job.stem
            self._log(f"\n{'─'*45}")
            self._log(f"  [{i+1}/{len(self.jobs)}] {stem}", "info")

            ok = True

            if job.run_es:
                job.status = JobStatus.RUNNING_ES
                self._update_job_ui(job)
                self._log(f"  ⚡ Electrostatic...", "info")
                try:
                    job.es_outputs = self._run_es(job)
                    self._log(f"  ⚡ ES done", "success")
                except Exception as e:
                    job.error = f"ES: {e}"
                    self._log(f"  ❌ ES: {e}", "error")
                    ok = False

            if job.run_hb and not self.stop_requested:
                job.status = JobStatus.RUNNING_HB
                self._update_job_ui(job)
                self._log(f"  💧 Hydrophobic...", "info")
                try:
                    job.hb_outputs = self._run_hb(job)
                    self._log(f"  💧 HB done", "success")
                except Exception as e:
                    job.error = (job.error + " | " if job.error else "") + f"HB: {e}"
                    self._log(f"  ❌ HB: {e}", "error")
                    ok = False

            if ok and self.gen_detailed_csv.get() and not self.stop_requested:
                try:
                    self._run_post(job, stem)
                except Exception as e:
                    self._log(f"  ⚠ Post: {e}", "warning")

            job.status = JobStatus.DONE if ok else JobStatus.FAILED
            self._update_job_ui(job)

        self.root.after(0, self._batch_done)

    def _run_post(self, job: JobItem, stem: str):
        try:
            import pandas as pd
            from src.unified_analyzer import analyze
            if job.run_es:
                pcsv = os.path.join(self.results_dir, f"{stem}_es_patches.csv")
                if os.path.exists(pcsv):
                    result, _ = analyze(pdb_path=job.pdb_path, patches_csv=pcsv,
                                        ply_pos=os.path.join(self.results_dir, f"{stem}_es-pos.ply"),
                                        ply_neg=os.path.join(self.results_dir, f"{stem}_es-neg.ply"), stem=f"{stem}_es")
                    out_csv = os.path.join(self.results_dir, f"{stem}_es_residues_detailed.csv")
                    result.to_csv(out_csv, index=False, encoding="utf-8-sig")
                    self._log(f"  [ES] Saved: {os.path.basename(out_csv)}")
                    self._save_patch_summary(result, f"{stem}_es", self.results_dir)
                    self._print_top_patches(result, "ES")
            if job.run_hb:
                npz = os.path.join(self.results_dir, f"{stem}_hb_out.npz")
                if os.path.exists(npz):
                    result, _ = analyze(pdb_path=job.pdb_path, npz_path=npz, stem=f"{stem}_hb")
                    out_csv = os.path.join(self.results_dir, f"{stem}_hb_residues_detailed.csv")
                    result.to_csv(out_csv, index=False, encoding="utf-8-sig")
                    self._log(f"  [HB] Saved: {os.path.basename(out_csv)}")
                    self._save_patch_summary(result, f"{stem}_hb", self.results_dir)
                    self._print_top_patches(result, "HB")
        except Exception as e:
            self._log(f"  Post: {e}", "warning")

    def _save_patch_summary(self, df, prefix: str, out_dir: str):
        """Generate and save patch-level summary CSV."""
        try:
            if "patch_total_area_A2" not in df.columns:
                return
            patch_summary = (
                df.groupby(["patch_nr", "patch_type"])
                .agg(
                    patch_total_area_A2=("patch_total_area_A2", "first"),
                    n_residues=("res_id", "nunique"),
                    top_residue=("res_id", "first"),
                    top_frac=("frac_of_patch", "first"),
                )
                .reset_index()
                .sort_values(["patch_type", "patch_nr"])
            )
            summary_file = os.path.join(out_dir, f"{prefix}_patch_summary.csv")
            patch_summary.to_csv(summary_file, index=False, encoding="utf-8-sig")
            self._log(f"  [{prefix.split('_')[0].upper()}] Saved: {os.path.basename(summary_file)}")
        except Exception:
            pass

    def _print_top_patches(self, df, label: str):
        """Print top-3 patches with type label."""
        if label == "ES":
            groups = [(["positive"], "Positive"), (["negative"], "Negative")]
        else:
            groups = [(["hydrophobic"], "Hydrophobic"), (["hydrophilic"], "Hydrophilic")]
        for ptype_list, pname in groups:
            subset = df[df.patch_type.isin(ptype_list)]
            if len(subset) == 0: continue
            self._log(f"  [{label}] {pname} top:", "header")
            seen, count = set(), 0
            for _, row in subset.iterrows():
                pn = row["patch_nr"]
                if pn in seen: continue
                seen.add(pn); count += 1
                if count > 3: break
                self._log(f"    #{pn}  {row.get('patch_total_area_A2','?')} A2  <- {row.get('chain_id','')}:{row.get('res_id','')} ({row.get('frac_of_patch',0)*100:.0f}%)")
            if len(seen) > 3:
                self._log(f"    ... and {len(seen)-3} more")

    def _run_es(self, job: JobItem) -> dict:
        from surface_analyses.commandline_electrostatic import run_electrostatics
        from surface_analyses.structure import load_trajectory_using_commandline_args
        from surface_analyses.platform_config import get_config as _get_cfg

        stem, pdb = job.stem, job.pdb_path
        es_prefix = f"{stem}_es"
        cfg = _get_cfg()
        cfg.setup_path()  # ensure pdb2pqr/apbs are on PATH
        project_root = str(cfg.project_root) if cfg else os.getcwd()

        class Args: pass
        args = Args()
        args.parm = pdb; args.trajs = [pdb]; args.stride = 1
        args.ref = None; args.protein_ref = None
        args.dx = None
        apbs = self.es_apbs_dir.get().strip()
        args.apbs_dir = apbs if apbs else os.path.join(project_root, f"Tools/apbs_batch_{stem}")
        args.probe_radius = float(self.es_probe_radius.get())
        args.out = os.path.join(self.results_dir, f"{es_prefix}_patches.csv")
        args.resout = None
        args.patch_cutoff = (float(self.es_patch_cutoff_pos.get()), float(self.es_patch_cutoff_neg.get()))
        args.integral_cutoff = (0.3, -0.3)
        args.surface_type = self.es_surface_type.get()
        args.ply_out = os.path.join(self.results_dir, es_prefix) if self.gen_ply.get() else None
        args.pos_patch_cmap = "tab20c"; args.neg_patch_cmap = "tab20c"
        args.ply_cmap = "coolwarm_r"; args.ply_clim = None
        args.check_cdrs = self.es_check_cdr.get()
        args.n_patches = 0; args.size_cutoff = 0.0
        args.gauss_shift = 0.1; args.gauss_scale = 1.0
        args.pH = float(self.es_ph.get()) if self.es_ph_enabled.get() else None
        args.ion_species = None

        traj = load_trajectory_using_commandline_args(args)
        del args.parm, args.trajs, args.stride, args.ref, args.protein_ref

        log_path = os.path.join(self.results_dir, f"{es_prefix}_run.log")
        return self._run_with_log(log_path, run_electrostatics, traj, **vars(args))

    def _run_hb(self, job: JobItem) -> dict:
        from surface_analyses.commandline_hydrophobic import run_hydrophobic
        from surface_analyses.structure import load_trajectory_using_commandline_args
        from surface_analyses.platform_config import get_config as _get_cfg

        stem, pdb = job.stem, job.pdb_path
        hb_prefix = f"{stem}_hb"
        cfg = _get_cfg()
        project_root = str(cfg.project_root) if cfg else os.getcwd()

        scale = self.hb_scale.get().strip()
        if not os.path.isabs(scale):
            scale = os.path.join(project_root, scale)

        class Args: pass
        args = Args()
        args.parm = pdb; args.trajs = [pdb]; args.stride = 1
        args.ref = None; args.protein_ref = None
        args.scale = scale
        args.smiles = None; args.atom_propensities = None
        args.out = os.path.join(self.results_dir, f"{hb_prefix}_out.npz")
        args.surftype = self.hb_surftype.get()
        args.group_heavy = False
        args.surfscore = False
        args.sap = self.hb_compute_sap.get()
        args.blur_rad = float(self.hb_blur_rad.get())
        args.sh = self.hb_compute_sh.get()
        args.sh_rad = float(self.hb_sh_rad.get())
        args.potential = self.hb_compute_potential.get()
        args.rmax = 0.3
        args.solv_rad = float(self.hb_solv_rad.get())
        args.grid_spacing = float(self.hb_grid_spacing.get())
        args.rcut = float(self.hb_rcut.get())
        args.alpha = float(self.hb_alpha.get())
        args.blur_sigma = 0.6
        args.ply_out = os.path.join(self.results_dir, hb_prefix) if self.gen_ply.get() else None
        args.ply_cmap = None; args.ply_clim = None
        args.patches = self.hb_compute_patches.get()
        args.patch_min = float(self.hb_patch_min.get())
        args.verbose = False

        traj = load_trajectory_using_commandline_args(args)
        parm = args.parm
        del args.parm, args.trajs, args.stride, args.ref, args.protein_ref

        log_path = os.path.join(self.results_dir, f"{hb_prefix}_run.log")
        return self._run_with_log(log_path, run_hydrophobic, parm, traj, **vars(args))

    def _run_with_log(self, log_path: str, func, *args, **kwargs) -> dict:
        """Run function with stdout/stderr captured to log file and GUI."""
        os.makedirs(os.path.dirname(log_path), exist_ok=True)

        class TeeWriter:
            def __init__(self, gui_log, file_handle):
                self.gui = gui_log
                self.fh = file_handle
            def write(self, s):
                if s:
                    for line in s.split("\n"):
                        stripped = line.rstrip()
                        if stripped:
                            self.gui("    " + stripped[:200])
                    self.fh.write(s)
            def flush(self):
                self.fh.flush()

        with open(log_path, "w", encoding="utf-8", buffering=1) as log_f:
            tee = TeeWriter(self._log, log_f)
            try:
                with redirect_stdout(tee), redirect_stderr(tee):
                    result = func(*args, **kwargs)
                log_f.write("\n[COMPLETED]\n")
            except Exception as e:
                log_f.write(f"\n[ERROR: {e}]\n")
                raise

        outputs = {"_run.log": log_path}
        for pat in ["*_patches.csv", "*-pos.ply", "*-neg.ply", "*-potential.ply",
                     "*_out.npz", "*_residues_detailed.csv", "*_patch_summary.csv"]:
            for f in globmod.glob(os.path.join(self.results_dir, pat)):
                outputs[os.path.basename(f)] = f
        return outputs

    def _batch_done(self):
        self.running = False
        self.process = None
        self.run_btn.config(state=tk.NORMAL)
        self.stop_btn.config(state=tk.DISABLED)
        self.batch_progress["value"] = self.batch_progress["maximum"]
        self.batch_label.config(text="")
        done = sum(1 for j in self.jobs if j.status == JobStatus.DONE)
        fail = sum(1 for j in self.jobs if j.status == JobStatus.FAILED)
        self._log(f"\n{'='*50}", "header")
        self._log(f"  Done: {done} OK  {fail} FAIL  ({len(self.jobs)} total)", "header")
        self._log(f"{'='*50}\n", "header")
        if done > 1:
            try:
                from src.batch_summary import summarize_electrostatic, summarize_hydrophobic, write_csv
                es_rows, hb_rows = [], []
                for job in self.jobs:
                    if job.status != JobStatus.DONE: continue
                    if job.run_es:
                        s = summarize_electrostatic(job.stem, self.results_dir)
                        if s: es_rows.append(s)
                    if job.run_hb:
                        s = summarize_hydrophobic(job.stem, self.results_dir)
                        if s: hb_rows.append(s)
                if es_rows:
                    write_csv(es_rows, os.path.join(self.results_dir, "batch_summary_es.csv"))
                    self._log(f"  Batch ES: batch_summary_es.csv", "success")
                if hb_rows:
                    write_csv(hb_rows, os.path.join(self.results_dir, "batch_summary_hb.csv"))
                    self._log(f"  Batch HB: batch_summary_hb.csv", "success")
            except Exception as e:
                self._log(f"  Summary: {e}", "warning")
        self._update_files_all()

    def _update_files_all(self):
        self.files_list.delete(0, tk.END)
        seen = set()
        # Scan results directory for all known output patterns
        all_patterns = [
            "*_run.log", "*_patches.csv", "*-pos.ply", "*-neg.ply", "*-potential.ply",
            "*_out.npz", "*_residues_detailed.csv", "*_patch_summary.csv"
        ]
        for pat in all_patterns:
            for f in globmod.glob(os.path.join(self.results_dir, pat)):
                name = os.path.basename(f)
                if f not in seen and os.path.exists(f):
                    seen.add(f)
                    sz = os.path.getsize(f)
                    ss = f"{sz/1e6:.1f}MB" if sz > 1e6 else f"{sz/1e3:.1f}KB" if sz > 1e3 else f"{sz}B"
                    # Try to match stem from filename
                    stem = name.split("_es_")[0].split("_hb_")[0]
                    self.files_list.insert(tk.END, f"[{stem}] {name} ({ss})")

    def _open_file(self, event):
        sel = self.files_list.curselection()
        if sel:
            text = self.files_list.get(sel[0])
            name = text.split("] ", 1)[1].split(" (")[0] if "] " in text else text
            path = os.path.join(self.results_dir, name)
            if os.path.exists(path):
                self._open_path(path)

    # ═══════════════════════════════════════════════════════
    def run(self):
        self.root.mainloop()


def main():
    PepPatchGUI().run()


if __name__ == "__main__":
    main()
