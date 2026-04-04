#!/bin/bash
# TADF Rate Constant Calculator — double-click to launch
cd "$(dirname "$0")"
source ~/.zshrc 2>/dev/null
echo "Starting TADF Rate Constant Calculator..."
python3 src/ui/flask_web/app.py
