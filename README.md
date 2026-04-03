# RateConstants_TADF

Intrinsic rate constants extraction for thermally activated delayed fluorescence (TADF) materials.

<p align="center">
<img src="https://github.com/d04943016/RateConstants_TADF/blob/main/Graph/process.png" width="600">
</p>

## Project Structure

```
RateConstants_TADF/
├── README.md
├── requirements.txt
├── data/                    # Calculation output data
└── src/
    ├── engine/              # Core calculation engine
    │   ├── calculator.py    # Rate constants calculator
    │   └── rate_equation.py # Numerical rate equations
    └── ui/
        ├── flask_web/       # Flask web interface  (port 5000)
        ├── fastapi/         # FastAPI REST API     (port 8000)
        └── pyqt5/           # PyQt5 desktop GUI
```

## Quick Start

Install all dependencies:
```bash
pip install -r requirements.txt
```

Or install per UI (see each ui's requirements.txt).

---

## UI Interfaces

### 1. Flask Web UI

Browser-based interface with interactive plots.

```bash
pip install -r src/ui/flask_web/requirements.txt
python src/ui/flask_web/app.py
# Open: http://localhost:5000
```

Config: `src/ui/flask_web/config/config.json`

---

### 2. FastAPI REST API

REST API backend with auto-generated docs.

```bash
pip install -r src/ui/fastapi/requirements.txt
python src/ui/fastapi/app.py
# API:  http://localhost:8000
# Docs: http://localhost:8000/docs
```

Config: `src/ui/fastapi/config/config.json`

**Endpoints:**

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/` | Welcome message |
| POST | `/tau2k` | Lifetime → rate constant |
| POST | `/k2tau` | Rate constant → lifetime |
| POST | `/cal_phi_PF_DF` | Calculate φ_PF and φ_DF |
| POST | `/cal_determined_intrinsic_rate_constants` | Directly determinable rate constants |
| POST | `/cal_intrinsic_rate_constants` | All intrinsic rate constants |
| POST | `/script` | Full analysis script |

Example:
```bash
python src/ui/fastapi/example_client.py
```

---

### 3. PyQt5 Desktop GUI

Standalone desktop application.

```bash
pip install -r src/ui/pyqt5/requirements.txt
python src/ui/pyqt5/main.py
```

Config: `src/ui/pyqt5/config/config.json`

---

## Core Engine

### Usage

```python
import sys
sys.path.insert(0, '/path/to/RateConstants_TADF')
from src.engine import calculator as engine

# Convert lifetime to rate constant
kPF, kDF = engine.tau2k([15e-9, 2.9e-6])

# Calculate intrinsic rate constants
import numpy as np
phi_Tnr_PL = np.linspace(0, 0.18, 200)
ks, ksr, ksnr, kisc, kt, ktr, ktnr, krisc = engine.cal_intrinsic_rate_constants(
    kPF, kDF, phi_PF=0.70, phi_DF=0.12, phi_Tnr_PL=phi_Tnr_PL
)

# Full script with plots and file output
engine.script(tau_PF=15e-9, tau_DF=2.9e-6, phi_PF=0.70, phi_DF=0.12,
              fpath='data/output', fname='DPAC-TRZ')
```

### Key Symbols

| Symbol | Description |
|--------|-------------|
| PF / DF | Prompt / Delayed fluorescence |
| τ (tau) | Exciton lifetime |
| k | Rate constant |
| ksr | Singlet radiative rate constant |
| ksnr | Singlet non-radiative rate constant |
| kisc | Intersystem crossing (S1→T1) |
| krisc | Reverse ISC (T1→S1) |
| PLQY | Photoluminescence quantum yield |

## Output

Results are saved to `data/output/` by default (configurable per UI in `config/config.json`).

## References

1. https://onlinelibrary.wiley.com/doi/abs/10.1002/adfm.201602501
2. https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201601675
3. https://pubs.rsc.org/no/content/articlelanding/2015/cc/c5cc05022g/unauth
4. https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201704961
