' PEP-Patch GUI 无控制台启动器
' 双击此文件即可启动，不会弹出命令行窗口

Set objShell = CreateObject("WScript.Shell")
scriptDir = objShell.CurrentDirectory

' 激活 venv 并启动 GUI
cmd = "cmd /c cd /d """ & scriptDir & """ && call .venv\Scripts\activate.bat && python -m surface_analyses.pep_patch_gui"
objShell.Run cmd, 0, False
