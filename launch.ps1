# PEP-Patch 快速启动器 (PowerShell)
# 双击此文件或在终端运行: .\launch.ps1
# 首次使用可能需要: Set-ExecutionPolicy -Scope CurrentUser RemoteSigned

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

# 激活虚拟环境
$activateScript = Join-Path $ScriptDir ".venv\Scripts\Activate.ps1"
if (Test-Path $activateScript) {
    . $activateScript
} else {
    Write-Host "Virtual environment not found. Please run setup first." -ForegroundColor Red
    Write-Host "Run: uv venv && uv sync" -ForegroundColor Yellow
    Read-Host "Press Enter to exit"
    exit 1
}

# 检查关键依赖
python -c "import tkinter" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "tkinter not available. Please install Python with tkinter support." -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

Write-Host "Starting PEP-Patch GUI..." -ForegroundColor Cyan
python -m surface_analyses.pep_patch_gui
