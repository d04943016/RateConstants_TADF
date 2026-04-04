from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import List, Optional
import json
import numpy as np
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.engine import calculator as engine

_CFG_PATH = os.path.join(os.path.dirname(__file__), 'config', 'config.json')
with open(_CFG_PATH) as _f:
    _CFG = json.load(_f)

OUTPUT_DIR = os.path.join(ROOT, _CFG.get('output_dir', 'data/output'))
os.makedirs(OUTPUT_DIR, exist_ok=True)

app = FastAPI(
    title="TADF Rate Constants Calculator API",
    description="API for calculating rate constants in Thermally Activated Delayed Fluorescence (TADF) materials"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class LifetimeInput(BaseModel):
    tau: List[float]

class RateConstantInput(BaseModel):
    k: List[float]

class PhiPFDFInput(BaseModel):
    PLQY: float
    tauPF: float
    tauDF: float
    B_PF: float
    B_DF: float

class DeterminedRateConstantsInput(BaseModel):
    kPF: float
    kDF: float
    phi_PF: float
    phi_DF: float

class IntrinsicRateConstantsInput(BaseModel):
    kPF: float
    kDF: float
    phi_PF: float
    phi_DF: float
    phi_Tnr_PL: List[float]

class ScriptInput(BaseModel):
    tau_PF: float
    tau_DF: float
    phi_PF: float
    phi_DF: float
    fname: Optional[str] = None


_TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), 'templates')

@app.get("/")
def serve_ui():
    return FileResponse(os.path.join(_TEMPLATE_DIR, 'index.html'))

@app.get("/api")
def read_root():
    return {"message": "Welcome to the TADF Rate Constants Calculator API"}

@app.post("/tau2k")
def convert_tau_to_k(input_data: LifetimeInput):
    """Convert lifetime to rate constant"""
    try:
        return {"k": engine.tau2k(input_data.tau).tolist()}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/k2tau")
def convert_k_to_tau(input_data: RateConstantInput):
    """Convert rate constant to lifetime"""
    try:
        return {"tau": engine.k2tau(input_data.k).tolist()}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/cal_phi_PF_DF")
def calculate_phi_pf_df(input_data: PhiPFDFInput):
    """Calculate the absolute ratio of prompt and delayed fluorescence"""
    try:
        phi_PF, phi_DF = engine.cal_phi_PF_DF(
            input_data.PLQY, input_data.tauPF, input_data.tauDF,
            input_data.B_PF, input_data.B_DF
        )
        return {"phi_PF": phi_PF, "phi_DF": phi_DF}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/cal_determined_intrinsic_rate_constants")
def calculate_determined_rate_constants(input_data: DeterminedRateConstantsInput):
    """Calculate rate constants determinable directly from experimental data"""
    try:
        ksr, ks, kt, kisckrisc = engine.cal_determined_intrinsic_rate_constants(
            input_data.kPF, input_data.kDF, input_data.phi_PF, input_data.phi_DF
        )
        return {"ksr": float(ksr), "ks": float(ks), "kt": float(kt), "kisckrisc": float(kisckrisc)}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/cal_intrinsic_rate_constants")
def calculate_intrinsic_rate_constants(input_data: IntrinsicRateConstantsInput):
    """Calculate all intrinsic rate constants"""
    try:
        ks, ksr, ksnr, kisc, kt, ktr, ktnr, krisc = engine.cal_intrinsic_rate_constants(
            input_data.kPF, input_data.kDF, input_data.phi_PF,
            input_data.phi_DF, input_data.phi_Tnr_PL
        )
        return {
            "ks": ks.tolist(), "ksr": ksr.tolist(), "ksnr": ksnr.tolist(), "kisc": kisc.tolist(),
            "kt": kt.tolist(), "ktr": ktr.tolist(), "ktnr": ktnr.tolist(), "krisc": krisc.tolist()
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/script")
def run_script(input_data: ScriptInput):
    """Run the rate constants calculator script"""
    try:
        fname = input_data.fname if input_data.fname else ""
        result = engine.script(
            input_data.tau_PF, input_data.tau_DF,
            input_data.phi_PF, input_data.phi_DF,
            fpath=OUTPUT_DIR, fname=fname
        )
        return {k: float(v) if hasattr(v, "item") else v for k, v in result.items()}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host=_CFG.get('host', '0.0.0.0'), port=_CFG.get('port', 8000))
