# PEP-Patch 批量表面分析 GUI — 使用说明

PEP-Patch 是一个对 PDB 文件进行**批量表面分析**的图形化工具，支持：

- ⚡ **静电势分析（Electrostatic, ES）**：基于 APBS 计算表面电势，按正负阈值提取 patch，可选检测抗体 CDR 区域。
- 💧 **疏水性分析（Hydrophobic, HB）**：基于原子级疏水量表（Wimley-White / Eisenberg / Crippen）计算疏水势与 patch。
- 📊 **批量结果汇总**：自动生成 `batch_summary_es.csv` / `batch_summary_hb.csv` 与残基级详细 CSV。
- 📁 **PLY 可视化**：可选导出正负电势、疏水势的 PLY 网格文件，便于 PyMOL/VMD 查看。

---

## 1. 启动方式

### Windows（推荐，双击即用）
- 双击 `PEP-Patch_GUI.bat` 启动（自动激活 `.venv` 并运行 GUI）。
- 或双击 `launch.vbs` / 在 PowerShell 运行 `.\launch.ps1`。

### Linux / macOS
`.bat` / `.ps1` / `.vbs` 为 Windows 专属启动器，在 Linux 上不可用。请改用命令行：

```bash
# 1) 激活虚拟环境（仓库根目录已有 activate.sh）
source activate.sh
# 或：source .venv/bin/activate

# 2) 直接运行模块
python -m surface_analyses.pep_patch_gui
```

> 提示：当前仓库**尚未提供** `launch.sh`。如需一键启动，可在仓库根目录自行创建：
> ```bash
> #!/usr/bin/env bash
> cd "$(dirname "$0")"
> source .venv/bin/activate
> python -m surface_analyses.pep_patch_gui
> ```

---

## 2. 界面与参数

界面分为左右两栏：

- **左侧（控制面板）**：添加 PDB 文件/目录 → 选择分析类型 → 配置 ES / HB 参数 → 设置输出 → 点击「开始批量分析」。
- **右侧（标签页）**：`日志`（实时终端输出）、`结果`（CSV 表格预览、可导出）、`输出文件`（生成文件列表，双击打开）。

### 静电势（ES）参数
| 参数 | 说明 |
| --- | --- |
| APBS 目录 | 留空 = 自动检测（`Tools/apbs_work`） |
| pH | 可选的等电点相关设置 |
| 表面类型 | `sas` / `ses` / `gauss` |
| 探针半径 (nm) | 默认 0.14 |
| 正/负阈值 | patch 提取阈值，默认 +2.0 / −2.0 |
| 检测 CDR | 抗体模式，标注 CDR 区域 |

### 疏水性（HB）参数
| 参数 | 说明 |
| --- | --- |
| 疏水表量表 | Wimley-White / Eisenberg / Crippen（含预设下拉） |
| 归一化 | `normal` / `sc_norm` / `atom_norm` |
| 计算项 | 疏水势 / Patch / SAP / SH |
| 模糊半径 / SH 半径 / 溶剂半径 | 网格平滑与溶剂参数 |
| 网格间距 / rcut / alpha / patch 最小面积 | 计算与 patch 提取参数 |

### 输出
- **生成 PLY 可视化文件**：导出网格 PLY。
- **生成残基级详细 CSV**：由 `src.unified_analyzer` 汇总每个 patch 的残基组成。

---

## 3. 平台兼容性（重点）

**结论：核心工具是跨平台的；GUI 此前唯一的 Windows 专属调用（`os.startfile`）已修复，现在 Linux 也可正常运行。**

### 3.1 已确认跨平台的部分
- **核心库** `surface_analyses.*` 与 `src.*`：纯 Python，无平台依赖。
- **平台自动探测** `surface_analyses/platform_config.py`：自动识别 `windows` / `wsl` / `linux`，并据此配置外部工具路径。
- **外部工具调用**：`platform_config` 按「环境变量 `PEP_PATCH_*` > 仓库内置 `Tools/` 二进制 > 系统 `PATH`」顺序解析 `msms` / `apbs` / `pdb2pqr` / `anarci`。**本地化二进制优先**，因此即使 uv 未安装 git 源包装包（如 `msms-wrapper`）也不会找不到工具。
- **GUI 框架**：基于 `tkinter` + `ttk`（clam 主题），本身是跨平台的。

### 3.2 已修复的兼容性问题
GUI 中「双击打开 PDB / 输出文件」原本直接调用 `os.startfile()`，该函数在 Linux / macOS 上**不存在**，会导致 `AttributeError` 崩溃。现已替换为跨平台实现：

```python
# Windows -> os.startfile ; macOS -> open ; Linux -> xdg-open
PEP_Patch_GUI._open_path(path)
```

### 3.3 Linux 使用时的注意事项
1. **启动器**：不要用 `.bat` / `.ps1` / `.vbs`，改用 `source activate.sh && python -m surface_analyses.pep_patch_gui`（见第 1 节）。
2. **外部工具（关键澄清）**：仓库 `Tools/` 内置的是 **Windows 二进制**（`APBS-3.4.1.Windows` 等），在 Linux 上无法运行、会被自动忽略。解析顺序为「环境变量 > `Tools/` 本地二进制 > `PATH`」，因此**本地化副本始终优先**，是否经 uv 安装 Python 包装包不影响可用性：
   - ✅ `pdb2pqr`：已在 `pyproject.toml` 主依赖声明（`pdb2pqr>=3`），`uv sync` 会装入 `.venv`；但运行时优先用 `Tools/pdb2pqr-portable/`，Linux 上则回退到系统 `pdb2pqr`。
   - ⚠️ `msms`：其 Python 包装包 `msms-wrapper`（git 源，**构建风险高**）已改为**可选 extra**，不再阻塞 `uv sync`（`uv sync --extra tools` 才安装）。运行时直接使用 `Tools/msms/msms(.exe)`，无需该包装包。
   - ⚠️ `apbs`：**外部 C/Fortran 二进制，未纳入 pyproject 依赖**，且仓库内置为 Windows 版。这是 Linux 上**真正的缺口**——需另行安装 Linux 原生 APBS（系统包管理器 / conda / `pip install apbs`）并确保其在 `PATH` 中。
   - ⚠️ `anarci`：pip 可装（`pip install anarci`），但未在 pyproject 声明；仅抗体 CDR 编号用到，属**可选项**。需要时在 venv 内 `uv pip install anarci` 即可。
3. **WSL 用户**：已显式支持。WSL 下建议直接在 WSL 内安装上述 Linux 版工具，`platform_config` 会按 `wsl` 分支配置（APBS 工作目录可用环境变量 `PEP_PATCH_APBS_DIR` 覆盖）。

---

## 4. 依赖与外部工具

- Python ≥ 3.11（仓库 `.venv` 为 3.12）。
- 运行时依赖见 `pyproject.toml` / `requirements.txt`（tkinter 通常随 Python 自带；若缺失需 `apt install python3-tk` 等）。
- 外部二进制：`MSMS`、`APBS`、`pdb2pqr`、`ANARCI`，位于 `Tools/`（Windows 版）。Linux 请另行安装对应版本。

### 4.1 获取外部组件（Tools / models / test 数据）

> 为减小仓库体积，以下大文件**不再纳入 git**（见根目录 `.gitignore`），仓库中仅保留**空目录占位**（`.gitkeep`）：
> `models/`（预训练模型）、`test/trastuzumab/`（夹具数据）。
> `Tools/` 中 **`msms`（~2 MB）已随仓库提供**，其余大型工具（APBS、pdb2pqr、ANARCI）仍需手动获取。

#### 4.1.1 外部工具（Tools/）

运行检查脚本可打印就绪状态与获取指引：
```bash
python scripts/setup_tools.py            # 检查并打印获取指引
python scripts/setup_tools.py --install # 额外尝试 pip 安装 pdb2pqr / anarci
```

| 工具 | 期望放置路径（相对仓库根） | 获取方式 |
|---|---|---|
| MSMS | `Tools/msms/msms.exe`（Windows）<br>`Tools/msms/msms`（Linux/macOS） | ✅ **已随仓库提供**（`msms.exe` + `cygwin1.dll`）。Linux/macOS 用户需自行获取对应平台的 `msms` 二进制放入该目录。 |
| APBS | `Tools/APBS-3.4.1.Windows/bin/apbs.exe`（Windows） | Windows: 下载 APBS-3.4.1 Windows 版解压到 `Tools/APBS-3.4.1.Windows/`（`apbs.exe` 在 `bin/` 下），来源 https://github.com/Electrostatics/APBS/releases 。Linux/macOS: `apt install apbs` 或 `conda install -c conda-forge apbs`。 |
| pdb2pqr | `Tools/pdb2pqr-portable/pdb2pqr.exe`（Windows） | Windows: 用便携版 `pdb2pqr-portable`（含 `pdb2pqr.exe`）放到 `Tools/pdb2pqr-portable/`。Linux/macOS: `pip install "pdb2pqr>=3"`。 |
| ANARCI | `Tools/ANARCI/anarci.exe`（Windows） | Windows: 下载 `anarci.exe` 放到 `Tools/ANARCI/`（仅抗体 CDR 编号用到，可选）。Linux/macOS: `pip install anarci`。 |

> 解析顺序：`platform_config` 按「环境变量 `PEP_PATCH_*` > 仓库内置 `Tools/` 二进制 > 系统 `PATH`」解析上述工具，因此本地 `Tools/` 副本始终优先。

#### 4.1.2 预训练模型（models/）

| 数据 | 期望放置路径 | 说明 |
|---|---|---|
| ImmuneBuilder 模型 | `models/immunebuilder/` | 预训练模型 `.pdb` 文件，请本地放置到 `models/immunebuilder/`（或从内部数据源获取）。 |

#### 4.1.3 测试夹具数据（test/trastuzumab/）

| 数据 | 期望放置路径 | 说明 |
|---|---|---|
| trastuzumab 测试数据 | `test/trastuzumab/` | 测试输入/结果（`.npz` / `.dx` / `.pdb` / `.parm7` / `.rst7` / `.save` / `.csv` 等），请本地放置到 `test/trastuzumab/`，或运行测试时由程序自动生成。 |

---

## 5. 已知限制

- 原生 Linux 缺少 `launch.sh`，需手动用模块命令启动（见 3.3）。
- 内置 `Tools/` 为 Windows 二进制，Linux 上被忽略（解析时回退到系统工具）。`pdb2pqr` 经 `uv sync` 安装、`msms` 由本地 `Tools/` 提供；**唯有 `apbs` 需另行安装 Linux 原生版本**，否则静电势（ES）分析会失败。`anarci` 仅抗体编号用到，可选安装。
- 大批量 PDB 分析为串行执行，日志与文件列表在完成后刷新。
