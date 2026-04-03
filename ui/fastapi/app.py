from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
import numpy as np
import os
import sys

# Add project root to path so src/ can be imported
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src import rate_constants_calculator as rcc

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


@app.get("/")
def read_root():
    return {"message": "Welcome to the TADF Rate Constants Calculator API"}

@app.post("/tau2k")
def convert_tau_to_k(input_data: LifetimeInput):
    """Convert lifetime to rate constant"""
    try:
        return {"k": rcc.tau2k(input_data.tau).tolist()}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/k2tau")
def convert_k_to_tau(input_data: RateConstantInput):
    """Convert rate constant to lifetime"""
    try:
        return {"tau": rcc.k2tau(input_data.k).tolist()}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/cal_phi_PF_DF")
def calculate_phi_pf_df(input_data: PhiPFDFInput):
    """Calculate the absolute ratio of prompt and delayed fluorescence"""
    try:
        phi_PF, phi_DF = rcc.cal_phi_PF_DF(
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
        ksr, ks, kt, kisckrisc = rcc.cal_determined_intrinsic_rate_constants(
            input_data.kPF, input_data.kDF, input_data.phi_PF, input_data.phi_DF
        )
        return {"ksr": float(ksr), "ks": float(ks), "kt": float(kt), "kisckrisc": float(kisckrisc)}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/cal_intrinsic_rate_constants")
def calculate_intrinsic_rate_constants(input_data: IntrinsicRateConstantsInput):
    """Calculate all intrinsic rate constants"""
    try:
        ks, ksr, ksnr, kisc, kt, ktr, ktnr, krisc = rcc.cal_intrinsic_rate_constants(
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
        result = rcc.script(
            input_data.tau_PF, input_data.tau_DF,
            input_data.phi_PF, input_data.phi_DF,
            fpath="", fname=fname
        )
        return {k: float(v) if hasattr(v, "item") else v for k, v in result.items()}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
