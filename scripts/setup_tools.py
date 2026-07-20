#!/usr/bin/env python3
"""检查并准备 PEP-Patch 的外部工具、模型与测试数据。

这些大体积组件**不再纳入 git 仓库**（见根目录 `.gitignore`），因此 clone
之后需要本脚本确认它们已就绪，或按提示下载 / 本地放置。

用法
----
    python scripts/setup_tools.py            # 仅检查并打印获取指引
    python scripts/setup_tools.py --install # 额外尝试用 pip 安装可装的组件

说明
----
* Windows 下优先使用仓库内置的 `Tools/` 二进制（由 `platform_config` 解析）；
  缺失时按下方指引下载并解压到对应目录即可。
* Linux / macOS 下 `apbs` / `pdb2pqr` / `anarci` / `msms` 通常直接由系统
  包管理器或 pip 提供，`platform_config` 会在 `Tools/` 找不到时回退到 PATH。
* `models/` 与 `test/` 中的大文件为项目数据，请本地放置（或从内部源获取），
  本脚本仅检查其是否存在。
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TOOLS = REPO_ROOT / "Tools"
MODELS = REPO_ROOT / "models"
TEST = REPO_ROOT / "test"

IS_WIN = os.name == "nt"


def _exists(p: Path) -> bool:
    return p.exists()


# 各外部工具的期望路径与获取方式。
# windows / posix 分别给出仓库内置二进制的期望位置；
# pip 为可经 pip 安装的包名（用于 --install）；
# note 为人工获取指引。
TOOL_SPECS = [
    {
        "name": "MSMS",
        "windows": TOOLS / "msms" / "msms.exe",
        "posix": TOOLS / "msms" / "msms",
        "pip": None,
        "note": (
            "Windows: 下载 msms 二进制（msms.exe + cygwin1.dll），"
            "放置到 Tools/msms/。来源: http://www.scripps.edu/~sanner/software/msms_home.html "
            "或 MGLTools。Linux/macOS: `pip install msms-wrapper` 或系统包。"
        ),
    },
    {
        "name": "APBS",
        "windows": TOOLS / "APBS-3.4.1.Windows" / "bin" / "apbs.exe",
        "posix": None,  # Linux/macOS 走系统 PATH
        "pip": None,
        "note": (
            "Windows: 下载 APBS-3.4.1 Windows 版并解压到 Tools/APBS-3.4.1.Windows/ "
            "(apbs.exe 应在 bin/ 下)。来源: https://github.com/Electrostatics/APBS/releases 。"
            "Linux/macOS: `apt install apbs` 或 `conda install -c conda-forge apbs`。"
        ),
    },
    {
        "name": "pdb2pqr",
        "windows": TOOLS / "pdb2pqr-portable" / "pdb2pqr.exe",
        "posix": None,
        "pip": "pdb2pqr>=3",
        "note": (
            "Windows: 使用便携版 pdb2pqr-portable（含 pdb2pqr.exe），放置到 "
            "Tools/pdb2pqr-portable/。Linux/macOS: `pip install \"pdb2pqr>=3\"`。"
        ),
    },
    {
        "name": "ANARCI",
        "windows": TOOLS / "ANARCI" / "anarci.exe",
        "posix": None,
        "pip": "anarci",
        "note": (
            "Windows: 下载 anarci.exe 放置到 Tools/ANARCI/。仅抗体 CDR 编号用到，属可选项。"
            "Linux/macOS: `pip install anarci`。"
        ),
    },
]

DATA_SPECS = [
    {
        "name": "models/immunebuilder",
        "path": MODELS / "immunebuilder",
        "note": "预训练模型数据。请本地放置到 models/immunebuilder/（或从内部数据源获取）。",
    },
    {
        "name": "test 夹具数据",
        "path": TEST / "trastuzumab",
        "note": (
            "测试输入/结果数据（.npz/.dx/.pdb 等）。请本地放置到 test/，"
            "或运行测试时由程序生成。"
        ),
    },
]


def _check_tool(spec: dict) -> tuple[bool, str]:
    expected = spec["windows"] if IS_WIN else spec["posix"]
    if expected is None:
        # 走系统 PATH
        on_path = shutil.which(spec["name"].lower()) is not None
        if on_path:
            return True, f"系统 PATH 中已找到 {spec['name'].lower()}"
        return False, "未在 Tools/ 或系统 PATH 中找到（见下方获取指引）"
    if _exists(expected):
        return True, str(expected)
    return False, f"缺失: {expected}"


def _pip_install(pkg: str) -> bool:
    try:
        subprocess.run(
            [sys.executable, "-m", "pip", "install", pkg],
            check=True,
            capture_output=True,
            text=True,
        )
        return True
    except Exception as exc:  # noqa: BLE001
        print(f"    pip 安装 {pkg} 失败: {exc}")
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description="检查并准备 PEP-Patch 外部组件")
    parser.add_argument(
        "--install",
        action="store_true",
        help="尝试用 pip 安装可装的组件（pdb2pqr / anarci）",
    )
    args = parser.parse_args()

    print("=" * 70)
    print(f"PEP-Patch 组件检查  (仓库根: {REPO_ROOT})")
    print("=" * 70)

    all_ok = True

    print("\n[外部工具]")
    for spec in TOOL_SPECS:
        ok, detail = _check_tool(spec)
        mark = "[+]" if ok else "[-]"
        print(f"  {mark} {spec['name']:<10} {detail}")
        if not ok:
            all_ok = False
            if args.install and spec["pip"]:
                print(f"    -> 尝试 pip 安装 {spec['pip']} ...")
                if _pip_install(spec["pip"]):
                    ok2, _ = _check_tool(spec)
                    if ok2:
                        print("    -> 安装成功。")
                        continue
            print(f"    -> {spec['note']}")

    print("\n[项目数据]")
    for spec in DATA_SPECS:
        ok = _exists(spec["path"])
        mark = "[+]" if ok else "[-]"
        print(f"  {mark} {spec['name']:<22} {'存在' if ok else '缺失'}: {spec['path']}")
        if not ok:
            all_ok = False
            print(f"    -> {spec['note']}")

    print("\n" + "=" * 70)
    if all_ok:
        print("全部组件就绪。")
        return 0
    print("部分组件缺失。请按上述指引下载/本地放置后重试本脚本。")
    print("提示: 本地已存在这些文件时（未纳入 git），无需任何操作即可直接运行。")
    print("=" * 70)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
