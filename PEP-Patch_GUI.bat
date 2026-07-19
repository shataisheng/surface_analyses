@echo off
cd /d "%~dp0"
call .venv\Scripts\activate.bat
python -m surface_analyses.pep_patch_gui
if errorlevel 1 pause
