#!/usr/bin/env python3
"""Build standalone TADF Calculator.app using PyInstaller."""

import PyInstaller.__main__
import os

ROOT = os.path.dirname(os.path.abspath(__file__))

PyInstaller.__main__.run([
    os.path.join(ROOT, 'src', 'ui', 'flask_web', 'app.py'),
    '--name', 'TADF Calculator',
    '--onedir',
    '--windowed',                # .app bundle on macOS, no terminal window
    '--noconfirm',               # overwrite without asking
    # Bundle data files
    '--add-data', f'{ROOT}/src/ui/flask_web/templates:templates',
    '--add-data', f'{ROOT}/src/ui/flask_web/config:config',
    '--add-data', f'{ROOT}/src/engine:src/engine',
    # Hidden imports that PyInstaller might miss
    '--hidden-import', 'numpy',
    '--hidden-import', 'flask',
    '--hidden-import', 'src.engine.calculator',
    '--hidden-import', 'src.engine.rate_equation',
    # Output
    '--distpath', os.path.join(ROOT, 'dist'),
    '--workpath', os.path.join(ROOT, 'build'),
    '--specpath', ROOT,
])
