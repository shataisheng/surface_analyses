# PEP-Patch 环境激活脚本 (Windows PowerShell)
# 用法: . .\activate.ps1
# 自动检测平台并配置 PATH

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

# 激活 uv 虚拟环境
$venvActivate = Join-Path $ScriptDir ".venv\Scripts\Activate.ps1"
if (Test-Path $venvActivate) {
    & $venvActivate
} else {
    Write-Warning "Virtual environment not found at .venv\. Run 'uv venv' first."
}

# 自动检测并添加 Tools 子目录到 PATH
$ToolsDir = Join-Path $ScriptDir "Tools"

if (Test-Path $ToolsDir) {
    $pathsToAdd = @()
    
    # MSMS
    $msmsDir = Join-Path $ToolsDir "msms"
    if (Test-Path $msmsDir) { $pathsToAdd += $msmsDir }
    
    # APBS (支持多个版本)
    $apbsDirs = @(Get-ChildItem $ToolsDir -Directory -Filter "APBS*" 2>$null)
    foreach ($d in $apbsDirs) {
        $binDir = Join-Path $d.FullName "bin"
        if (Test-Path $binDir) { $pathsToAdd += $binDir }
    }
    
    # pdb2pqr
    $pdb2pqrDir = Join-Path $ToolsDir "pdb2pqr-portable"
    if (Test-Path $pdb2pqrDir) { $pathsToAdd += $pdb2pqrDir }
    
    # 添加所有找到的路径到 PATH
    if ($pathsToAdd.Count -gt 0) {
        $env:PATH = ($pathsToAdd -join ";") + ";" + $env:PATH
    }
}

# 显示平台信息
Write-Host ""
Write-Host "=== PEP-Patch Environment ===" -ForegroundColor Cyan
Write-Host "  Platform : Windows" -ForegroundColor White
Write-Host "  Python   : $(python --version 2>&1)" -ForegroundColor White

# 验证工具可用性
$msmsOk = Get-Command msms.exe -ErrorAction SilentlyContinue
$apbsOk = Get-Command apbs.exe -ErrorAction SilentlyContinue
$pdb2pqrOk = Get-Command pdb2pqr.exe -ErrorAction SilentlyContinue

if ($msmsOk) { Write-Host "  MSMS     : $($msmsOk.Source)" -ForegroundColor Green }
else { Write-Host "  MSMS     : NOT FOUND" -ForegroundColor Red }

if ($apbsOk) { Write-Host "  APBS     : $($apbsOk.Source)" -ForegroundColor Green }
else { Write-Host "  APBS     : NOT FOUND" -ForegroundColor Red }

if ($pdb2pqrOk) { Write-Host "  pdb2pqr  : $($pdb2pqrOk.Source)" -ForegroundColor Green }
else { Write-Host "  pdb2pqr  : NOT FOUND" -ForegroundColor Red }

Write-Host ""
Write-Host "Commands: pep_patch_hydrophobic | pep_patch_electrostatic" -ForegroundColor Yellow

