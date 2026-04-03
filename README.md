# RateConstants_TADF

Intrinsic rate constants extraction for thermally activated delayed fluorescence (TADF) materials with two delayed components.

This module calculates intrinsic rate constants between energy states: S0 (ground state), S1 (1st excited singlet), and T1 (1st excited triplet state).

<p align="center">
<img src="https://github.com/d04943016/RateConstants_TADF/blob/main/Graph/process.png" width="600">
</p>

## Project Structure

```
RateConstants_TADF/
├── src/                          # Core calculation module (new)
├── RateConstantCalculator.py     # Core calculation module (legacy)
├── ui/
│   ├── flask_web/                # Flask web UI
│   ├── fastapi/                  # FastAPI REST API
│   └── pyqt5/                    # PyQt5 desktop GUI
└── requirements.txt              # Combined dependencies
```

---

## User Interfaces

Three independent UIs are available — each can be run separately.

---

### 1. Flask Web UI

A browser-based interface built with Flask.

**Install dependencies:**
```bash
pip install -r ui/flask_web/requirements.txt
```

**Run:**
```bash
cd ui/flask_web
python app.py
```

**Open browser at:** `http://localhost:5000`

- Input PF/DF lifetimes (or rate constants) and quantum yields
- View rate constants and quantum yield plots interactively
- Download results as `.txt`

---

### 2. FastAPI REST API

A REST API backend built with FastAPI. Suitable for programmatic access or integration with other frontends.

**Install dependencies:**
```bash
pip install -r ui/fastapi/requirements.txt
```

**Run:**
```bash
cd ui/fastapi
python app.py
```

**API available at:** `http://localhost:8000`

**Interactive docs at:** `http://localhost:8000/docs`

**Main endpoints:**

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/tau2k` | Convert lifetime to rate constant |
| POST | `/k2tau` | Convert rate constant to lifetime |
| POST | `/cal_phi_PF_DF` | Calculate phi_PF and phi_DF from PLQY |
| POST | `/cal_determined_intrinsic_rate_constants` | Calculate determinable rate constants |
| POST | `/cal_intrinsic_rate_constants` | Calculate all intrinsic rate constants |
| POST | `/script` | Run full analysis script |

**Example usage:**
```bash
python ui/fastapi/example_client.py
```

---

### 3. PyQt5 Desktop GUI

A standalone desktop application built with PyQt5.

**Install dependencies:**
```bash
pip install -r ui/pyqt5/requirements.txt
```

**Run:**
```bash
cd ui/pyqt5
python GUI_PyQt5_main.py
```

**Features:**
- Select input mode: lifetime or rate constant
- Input PF/DF values and quantum yields (or prefactors + PLQY)
- Choose save path and filename
- View summarized results in terminal
- Export rate constants and quantum yield plots

<p align="center">
<img src="https://github.com/d04943016/RateConstants_TADF/blob/main/Graph/Panel.png" width="1200">
</p>

---

## Core Module

### Symbols

| Symbol | Description |
|--------|-------------|
| PF | Prompt fluorescence |
| DF | Delayed fluorescence |
| tau | Exciton lifetime |
| k | Rate constant |
| PLQY | Photoluminescence quantum yield |

### Key functions (`RateConstantCalculator`)

```python
import RateConstantCalculator as RCC

# Convert lifetime <-> rate constant
RCC.tau2k(tau)
RCC.k2tau(k)

# Calculate phi_PF and phi_DF
RCC.phi_PF_DF(PLQY, tauPF, tauDF, B_PF, B_DF)

# Calculate intrinsic rate constants
RCC.IntrinsicRateConstants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL)

# Full analysis scripts
RCC.script_for_100_PLQY(tau_PF, tau_DF, phi_PF, phi_DF, name='')
RCC.script(tau_PF, tau_DF, phi_PF, phi_DF, fpath='', fname='')
```

### Example

```python
import numpy as np
import RateConstantCalculator as RCC

# DPAC-TRZ
RCC.script(tau_PF=15e-9, tau_DF=2.9e-6, phi_PF=0.70, phi_DF=0.12,
           fname='DPAC-TRZ', fpath='./data')
```

---

## References

1. https://onlinelibrary.wiley.com/doi/abs/10.1002/adfm.201602501
2. https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201601675
3. https://pubs.rsc.org/no/content/articlelanding/2015/cc/c5cc05022g/unauth
4. https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201704961
