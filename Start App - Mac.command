#!/usr/bin/env bash
# ====================================================================
#  Cluster Analysis HTML - Mac launcher
#  Double-click this file in Finder to start the app.
#  (First time only: right-click -> Open, to bypass the unsigned-script
#  warning. After that, double-clicking works normally.)
# ====================================================================
cd "$(dirname "$0")" || exit 1

if command -v python3 >/dev/null 2>&1; then
    python3 server.py
elif command -v python >/dev/null 2>&1; then
    python server.py
else
    echo "Python 3 was not found. Install it from https://www.python.org/downloads/"
    read -p "Press Enter to close..."
    exit 1
fi
