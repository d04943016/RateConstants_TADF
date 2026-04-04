#!/usr/bin/env python3
# Copyright (c) 2020 Wei-Kai Lee. All rights reserved
# Web interface for TADF Rate Constant Calculator

import io
import json
import os
import sys
import time
import threading
import numpy as np
from flask import Flask, render_template, request, jsonify, send_file

# ── Path resolution (works both in dev and PyInstaller bundle) ─────────────────
def _resolve_base():
    """Return base directory for bundled resources."""
    if getattr(sys, 'frozen', False):
        return sys._MEIPASS  # PyInstaller temp dir
    return os.path.dirname(os.path.abspath(__file__))

_BASE = _resolve_base()

# Resolve project root (4 levels up from src/ui/flask_web/)
if getattr(sys, 'frozen', False):
    ROOT = os.path.dirname(sys.executable)  # .app/Contents/MacOS/
else:
    ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

if getattr(sys, 'frozen', False):
    sys.path.insert(0, _BASE)
    from src.engine import calculator as engine
else:
    from src.engine import calculator as engine

# Load config
_CFG_PATH = os.path.join(_BASE, 'config', 'config.json')
with open(_CFG_PATH) as _f:
    _CFG = json.load(_f)

OUTPUT_DIR = os.path.join(ROOT, _CFG.get('output_dir', 'data/output'))
os.makedirs(OUTPUT_DIR, exist_ok=True)

_TEMPLATE_DIR = os.path.join(_BASE, 'templates')
app = Flask(__name__, template_folder=_TEMPLATE_DIR)

# ── Heartbeat: auto-shutdown when browser closes ──────────────────────────────
_last_heartbeat = time.time()
_HEARTBEAT_TIMEOUT = 10  # seconds without heartbeat before shutdown


def _watchdog():
    """Background thread: shutdown server if no heartbeat received."""
    while True:
        time.sleep(3)
        if time.time() - _last_heartbeat > _HEARTBEAT_TIMEOUT:
            print('\n[TADF] No browser heartbeat — shutting down server.')
            os._exit(0)


@app.route('/api/heartbeat', methods=['POST'])
def heartbeat():
    global _last_heartbeat
    _last_heartbeat = time.time()
    return '', 204


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/calculate', methods=['POST'])
def calculate():
    try:
        data = request.get_json()
        tau_PF, tau_DF, Phi_PF, Phi_DF = _parse_inputs(data)
        n_points = max(10, int(data.get('n_points', 200)))
        results = _run_calculation(tau_PF, tau_DF, Phi_PF, Phi_DF, n_points)
        return jsonify({'status': 'ok', 'data': results})
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)}), 400


@app.route('/api/download', methods=['POST'])
def download():
    try:
        data = request.get_json()
        fname = data.get('filename', 'results').strip() or 'results'
        tau_PF, tau_DF, Phi_PF, Phi_DF = _parse_inputs(data)
        n_points = max(10, int(data.get('n_points', 200)))

        kPF, kDF = engine.tau2k([tau_PF, tau_DF])
        PLQY = Phi_PF + Phi_DF
        phi_Tnr_PL_Array = np.linspace(0, 1 - PLQY, n_points)

        (ks_Array, ksr_Array, ksnr_Array, kisc_Array,
         kt_Array, ktr_Array, ktnr_Array, krisc_Array) = \
            engine.cal_intrinsic_rate_constants(kPF, kDF, Phi_PF, Phi_DF,
                                                phi_Tnr_PL=phi_Tnr_PL_Array)

        phi_sr_Array, phi_snr_Array, phi_isc_Array = \
            engine.cal_phi_sr_snr_isc(ksr_Array, ksnr_Array, kisc_Array)
        phi_tr_Array, phi_tnr_Array, phi_risc_Array = \
            engine.cal_phi_tr_tnr_risc(ktr_Array, ktnr_Array, krisc_Array)

        buf = io.StringIO()
        header = ('{:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}'
                  ' {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}\n')
        buf.write(header.format(
            'ks(1/s)', 'ksr(1/s)', 'ksnr(1/s)', 'kisc(1/s)',
            'kt(1/s)', 'ktr(1/s)', 'ktnr(1/s)', 'krisc(1/s)',
            'phi_sr(%)', 'phi_snr(%)', 'phi_isc(%)',
            'phi_tr(%)', 'phi_tnr(%)', 'phi_risc(%)'))

        row_fmt = ('{:>15.5e} {:>15.5e} {:>15.5e} {:>15.5e}'
                   ' {:>15.5e} {:>15.5e} {:>15.5e} {:>15.5e}'
                   ' {:>15.5f} {:>15.5f} {:>15.5f}'
                   ' {:>15.5f} {:>15.5f} {:>15.5f}\n')
        for ii in range(len(phi_Tnr_PL_Array)):
            buf.write(row_fmt.format(
                ks_Array[ii], ksr_Array[ii], ksnr_Array[ii], kisc_Array[ii],
                kt_Array[ii], ktr_Array[ii], ktnr_Array[ii], krisc_Array[ii],
                phi_sr_Array[ii] * 100, phi_snr_Array[ii] * 100,
                phi_isc_Array[ii] * 100, phi_tr_Array[ii] * 100,
                phi_tnr_Array[ii] * 100, phi_risc_Array[ii] * 100))

        buf.seek(0)
        return send_file(
            io.BytesIO(buf.getvalue().encode()),
            mimetype='text/plain',
            as_attachment=True,
            download_name=f'{fname}.txt'
        )
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)}), 400


def _parse_inputs(data):
    mode_lf_k = int(data.get('mode_lf_k', 0))
    mode_qy_b = int(data.get('mode_qy_b', 0))

    if mode_lf_k == 0:  # lifetime
        tau_PF = float(data['pf_value']) * 1e-9
        tau_DF = float(data['df_value']) * 1e-6
    else:               # rate constant
        kPF = float(data['pf_value']) * 1e8
        kDF = float(data['df_value']) * 1e5
        taus = engine.k2tau([kPF, kDF])
        tau_PF, tau_DF = float(taus[0]), float(taus[1])

    if mode_qy_b == 0:  # quantum yield
        Phi_PF = float(data['phi_pf_value']) * 0.01
        Phi_DF = float(data['phi_df_value']) * 0.01
    else:               # B prefactors + PLQY
        B_PF = float(data['phi_pf_value'])
        B_DF = float(data['phi_df_value'])
        PLQY = float(data['plqy_value']) * 0.01
        Phi_PF, Phi_DF = engine.cal_phi_PF_DF(PLQY, tau_PF, tau_DF, B_PF, B_DF)

    return tau_PF, tau_DF, Phi_PF, Phi_DF


def _run_calculation(tau_PF, tau_DF, Phi_PF, Phi_DF, n_points=200):
    kPF, kDF = engine.tau2k([tau_PF, tau_DF])
    PLQY = Phi_PF + Phi_DF
    phi_Tnr_PL_Array = np.linspace(0, 1 - PLQY, n_points)

    (ks_Array, ksr_Array, ksnr_Array, kisc_Array,
     kt_Array, ktr_Array, ktnr_Array, krisc_Array) = \
        engine.cal_intrinsic_rate_constants(kPF, kDF, Phi_PF, Phi_DF,
                                            phi_Tnr_PL=phi_Tnr_PL_Array)

    phi_sr_Array, phi_snr_Array, phi_isc_Array = \
        engine.cal_phi_sr_snr_isc(ksr_Array, ksnr_Array, kisc_Array)
    phi_tr_Array, phi_tnr_Array, phi_risc_Array = \
        engine.cal_phi_tr_tnr_risc(ktr_Array, ktnr_Array, krisc_Array)

    x = (phi_Tnr_PL_Array * 100).tolist()

    def s(v):
        return float(v)

    return {
        'kPF': s(kPF), 'kDF': s(kDF),
        'PLQY': s(PLQY * 100),
        'Phi_PF': s(Phi_PF * 100), 'Phi_DF': s(Phi_DF * 100),
        'x': x,
        'rate_constants': {
            'ks':    ks_Array.tolist(),
            'ksr':   ksr_Array.tolist(),
            'ksnr':  ksnr_Array.tolist(),
            'kisc':  kisc_Array.tolist(),
            'kt':    kt_Array.tolist(),
            'ktr':   ktr_Array.tolist(),
            'ktnr':  ktnr_Array.tolist(),
            'krisc': krisc_Array.tolist(),
        },
        'quantum_yields': {
            'phi_sr':   (phi_sr_Array   * 100).tolist(),
            'phi_snr':  (phi_snr_Array  * 100).tolist(),
            'phi_isc':  (phi_isc_Array  * 100).tolist(),
            'phi_tr':   (phi_tr_Array   * 100).tolist(),
            'phi_tnr':  (phi_tnr_Array  * 100).tolist(),
            'phi_risc': (phi_risc_Array * 100).tolist(),
        },
        'summary': {
            'phi_tnr_0': {
                'ks':    s(ks_Array[0]),   'ksr':   s(ksr_Array[0]),
                'ksnr':  s(ksnr_Array[0]), 'kisc':  s(kisc_Array[0]),
                'kt':    s(kt_Array[0]),   'ktr':   s(ktr_Array[0]),
                'ktnr':  s(ktnr_Array[0]), 'krisc': s(krisc_Array[0]),
                'phi_sr':   s(phi_sr_Array[0]   * 100),
                'phi_snr':  s(phi_snr_Array[0]  * 100),
                'phi_isc':  s(phi_isc_Array[0]  * 100),
                'phi_tr':   s(phi_tr_Array[0]   * 100),
                'phi_tnr':  s(phi_tnr_Array[0]  * 100),
                'phi_risc': s(phi_risc_Array[0] * 100),
            },
            'phi_tnr_max': {
                'ks':    s(ks_Array[-1]),   'ksr':   s(ksr_Array[-1]),
                'ksnr':  s(ksnr_Array[-1]), 'kisc':  s(kisc_Array[-1]),
                'kt':    s(kt_Array[-1]),   'ktr':   s(ktr_Array[-1]),
                'ktnr':  s(ktnr_Array[-1]), 'krisc': s(krisc_Array[-1]),
                'phi_sr':   s(phi_sr_Array[-1]   * 100),
                'phi_snr':  s(phi_snr_Array[-1]  * 100),
                'phi_isc':  s(phi_isc_Array[-1]  * 100),
                'phi_tr':   s(phi_tr_Array[-1]   * 100),
                'phi_tnr':  s(phi_tnr_Array[-1]  * 100),
                'phi_risc': s(phi_risc_Array[-1] * 100),
            },
        },
    }


if __name__ == '__main__':
    import webbrowser

    host = _CFG.get('host', '0.0.0.0')
    port = _CFG.get('port', 5000)
    debug = _CFG.get('debug', True)
    auto_open = _CFG.get('auto_open', True)

    # Frozen bundle: disable debug + reloader (Werkzeug reloader spawns subprocesses
    # that confuse macOS windowed mode, causing "not responding" dialog)
    if getattr(sys, 'frozen', False):
        debug = False

    # Start watchdog (only when not in reloader child process)
    if auto_open and os.environ.get('WERKZEUG_RUN_MAIN') != 'true':
        threading.Thread(target=_watchdog, daemon=True).start()
        # Poll until server is ready, then open browser
        def _open_browser():
            import urllib.request
            browse_host = '127.0.0.1' if host == '0.0.0.0' else host
            url = f'http://{browse_host}:{port}/'
            for _ in range(30):   # up to ~6 seconds
                try:
                    urllib.request.urlopen(url, timeout=0.5)
                    break
                except Exception:
                    time.sleep(0.2)
            webbrowser.open(url)
        threading.Thread(target=_open_browser, daemon=True).start()

    app.run(debug=debug, host=host, port=port,
            use_reloader=False if getattr(sys, 'frozen', False) else debug)
