@echo off
REM ===================================================================
REM  Cluster Analysis HTML - Windows launcher
REM  Double-click this file to start the app.
REM  IMPORTANT: use the browser tab this opens automatically.
REM ===================================================================
cd /d "%~dp0"

where py >nul 2>nul
if %errorlevel%==0 (
    py server.py
) else (
    python server.py
)

pause
