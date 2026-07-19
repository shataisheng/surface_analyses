#!/bin/bash
# PEP-Patch 环境激活脚本 (WSL / Linux)
# 用法: source activate.sh
# 自动检测平台并配置 PATH

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 检测平台
if grep -qi microsoft /proc/version 2>/dev/null || grep -qi wsl /proc/version 2>/dev/null; then
    PLATFORM="WSL"
elif [[ "$(uname -s)" == "Linux" ]]; then
    PLATFORM="Linux"
elif [[ "$(uname -s)" == "Darwin" ]]; then
    PLATFORM="macOS"
else
    PLATFORM="Unknown"
fi

echo ""
echo "=== PEP-Patch Environment ==="
echo -e "  Platform : \033[36m${PLATFORM}\033[0m"

# 激活 uv 虚拟环境
if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then
    source "$SCRIPT_DIR/.venv/bin/activate"
    echo -e "  Python   : \033[32m$(python --version 2>&1)\033[0m"
else
    echo -e "  \033[33mWarning: .venv not found. Run 'uv venv && uv sync' first.\033[0m"
fi

# 添加本地 Tools 到 PATH（优先使用系统已安装的）
TOOLS_DIR="$SCRIPT_DIR/Tools"
if [ -d "$TOOLS_DIR" ]; then
    # MSMS
    if [ -d "$TOOLS_DIR/msms" ]; then
        export PATH="$TOOLS_DIR/msms:$PATH"
    fi
    # APBS
    for apbs_dir in "$TOOLS_DIR"/APBS*/bin; do
        if [ -d "$apbs_dir" ]; then
            export PATH="$apbs_dir:$PATH"
        fi
    done
    # pdb2pqr
    if [ -d "$TOOLS_DIR/pdb2pqr-portable" ]; then
        export PATH="$TOOLS_DIR/pdb2pqr-portable:$PATH"
    fi
fi

# 验证工具可用性
check_tool() {
    local name=$1
    local path
    path=$(which "$name" 2>/dev/null)
    if [ -n "$path" ]; then
        echo -e "  ${name} : \033[32m${path}\033[0m"
    else
        echo -e "  ${name} : \033[31mNOT FOUND\033[0m"
    fi
}

check_tool msms
check_tool apbs
check_tool pdb2pqr

echo ""
echo -e "\033[33mCommands: pep_patch_hydrophobic | pep_patch_electrostatic\033[0m"
echo ""

# Export platform info for scripts
export PEP_PATCH_PLATFORM="$PLATFORM"
