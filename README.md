# RateConstants_TADF

A module for extracting intrinsic rate constants in thermally activated delayed fluorescence (TADF) materials with two delayed components.

Prerequisites:
```
pip install numpy matplotlib
```

References:<br/>
[1] [Efficient and Tunable Thermally Activated Delayed Fluorescence Emitters Having Orientation-Adjustable CN-Substituted Pyridine and Pyrimidine Acceptor Units
](https://onlinelibrary.wiley.com/doi/abs/10.1002/adfm.201602501)<br/>
[2] [Sky-Blue Organic Light Emitting Diode with 37% External Quantum Efficiency Using Thermally Activated Delayed Fluorescence from Spiroacridine-Triazine Hybrid
](https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201601675)<br/>
[3] [A versatile thermally activated delayed fluorescence emitter for both highly efficient doped and non-doped organic light emitting devices](https://pubs.rsc.org/no/content/articlelanding/2015/cc/c5cc05022g/unauth#!divAbstract)<br/>
[4] [Achieving Nearly 30% External Quantum Efficiency for Orange–Red Organic Light Emitting Diodes by Employing Thermally Activated Delayed Fluorescence Emitters Composed of 1,8-Naphthalimide-Acridine Hybrids
](https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201704961)<br/>

This module calculates intrinsic rate constants between energy states:
- S0 (ground state)
- S1 (1st excited singlet)
- T1 (1st excited triplet state)

Including transitions:
- S1→S0 (fluorescence and heat loss)
- S1→T1 (intersystem crossing)
- T1→S1 (reverse intersystem crossing)
- T1→S0 (phosphorescence and heat loss)

<p align="center">
<img src="./docs/Graph/process.png" width="600">
</p>

The module also provides TADF material efficiency and exciton concentration calculations.

## Symbols
- PF: Prompt fluorescence
- DF: Delayed fluorescence
- tau: Exciton lifetime
- k: Rate constant

#### Measurables:
- tauPF/tauDF: PF/DF exciton lifetime from transient data
- kPF/kDF: PF/DF rate constant (calculated from tauPF/tauDF)
- PLQY: Photoluminescence quantum yield (thin film)

#### Intrinsic rate constants in S1
- ks: Total rate constant in S1 (ks = ksr + ksnr + kisc)
- ksr: Radiative rate constant in S1
- ksnr: Non-radiative rate constant in S1
- kisc: Intersystem rate constant (S1→T1)

#### Intrinsic rate constants in T1
- kt: Total rate constant in T1 (kt = ktr + ktnr + krisc)
- ktr: Radiative rate constant in T1
- ktnr: Non-radiative rate constant in T1
- krisc: Reverse intersystem rate constant (T1→S1)

#### Quantum yield (Q.Y.)
- phi_PF: PLQY component from PF (PL excitation, thin film)
- phi_DF: PLQY component from DF (PL excitation, thin film)

- phi_sr: Singlet radiative efficiency
- phi_snr: Singlet non-radiative efficiency
- phi_isc: Intersystem crossing efficiency

- phi_tr: Triplet radiative efficiency
- phi_tnr: Triplet non-radiative efficiency
- phi_risc: Reverse intersystem crossing efficiency

- phi_Tnr_PL: Total loss efficiency from triplet state in PL excitation
- phi_Snr_PL: Total loss efficiency from singlet state in PL excitation

#### Prefactor
- B_PF: Prefactor of PF
- B_DF: Prefactor of DF

#### Purcell Effect
F: Purcell factor
IQE = (PLQY * F/((1-PLQY)+PLQY * F)) * (alpha + phi_risc * (1-alpha))
where:
- alpha = 1.00 for PL excitation
- alpha = 0.25 for EL excitation

## RateConstantsCalculator Module

```python
import RateConstantsCalculator as RCC
```

#### Utility Functions
```python
RCC.tau2k(tau)          # Calculate rate constant from lifetime
RCC.k2tau(k)           # Calculate lifetime from rate constant
RCC.exponential_ratio(tau_Array, B_Array)  # Calculate relative ratio of each exponent
RCC.phi_PF_DF(PLQY, tauPF, tauDF, B_PF, B_DF)  # Calculate quantum efficiency from PF and DF
```

#### intrinsic rate constant calculator
`RCC.IntrinsicRateConstants_Determined(kPF, kDF, phi_PF, phi_DF)`<br/>
a function to calculate the intrinsice rate constants that can be determined by kPF, kDF, phi_PF, and phi_DF. <br/>
These four values (kPF, kDF, phi_PF, phi_DF) can be directly extracted from the transient data and PLQY (phi_PF, phi_DF)<br/>

    kPF=RCC.tau2k(tauPF)
    kDF=RCC.tau2k(tauDF) 
    
    [phi_PF, phi_DF] = RCC.phi_PF_DF(PLQY, tauPF, tauDF, B_PF, B_DF)
    ksr, ks, kt, kisckrisc = RCC.IntrinsicRateConstants_Determined(kPF, kDF, phi_PF, phi_DF)

output : ksr, ks, kt, kisckrisc (a product of the intersystem rate constant and the reverse intersystem crossing rate constant)<br/>

`RCC.IntrinsicRateConstants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL)`<br/>
a function to calculate the rate constants which cannot be calculated directly from the experimental data and should considering a loss ratio from triplet state by PL excitation (phi_Tnr_PL), which is a value between 0 to 1-PLQY<br/>
phi_Tnr_PL=0 indicates all the loss is from singlet state; in constrast, phi_Tnr_PL=1-PLQY indicates all the loss is from triplet state.
    
    import numpy as np
    phi_Tnr_PL = np.linspace(0, 1-PLQY, 200)
    ks, ksr, ksnr, kisc, kt, ktr, ktnr, krisc = RCC.IntrinsicRateConstants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL)

#### phi function
`RCC.phi_sr_snr_isc(ksr, ksnr, kisc)`<br/>
phi_sr_snr_isc is a function to calculate the efficiency of each process in S1<br/>
output : phi_sr, phi_snr, phi_isc

`RCC.phi_tr_tnr_risc(ktr, ktnr, krisc)`<br/>
phi_tr_tnr_risc is a function to calculate the efficiency of each process in T1<br/>
output : phi_tr, phi_tnr, phi_risc

#### Internal Quantum Efficiency (IQE)
`RCC.IQE_PurcellEffect(IQE, F=1.0)`<br/>
a function to calculate the internal quantum yield considering Purcell effect. <br/>
IQE : internal quantum yield (PL excitation)<br/>
F : Purcell factor<br/>

`RCC.IQE_RateConstants(ksr, kt, krisc, kPF, kDF, alpha=1.0, F=1.0)`<br/>
a function to calculate the internal quantum yield considering Purcell effect. <br/>

`RCC.IQE_Phi(phi_sr, phi_isc, phi_risc, alpha=1.0, F=1.0)`<br/>
a function to calculate the internal quantum efficiency given with intrinsic state efficiency<br/>

`RCC.IQE_PLQY(phi_risc, PLQY, alpha=1.0, F=1.0)`<br/>
 a function to calculate the internal quantum efficiency given with PLQY<br/>

`RCC.PLQY_phi(phi_sr, phi_isc, phi_risc)`<br/>
a function to calculate the PLQY from intrinsic state efficiency<br/>

#### Differential Equations
`RCC.k_total(kr, knr, kisc)`<br/>
a function to calculate the total rate constant in a state.<br/> 

`RCC.kPF_kDF(ksr, ksnr, kisc, ktr, ktnr, krisc, F=1.0)`<br/>
a function to calculate kPF and kDF from intrinsic rate constant considering Purcell Effect.<br/>
output : kPF, kDF<br/>

`RCC.Concentration_SteadyState(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0)`<br/>
a function to calculate the steady state S1 and T1 exciton concentration.<br/>
G : total exciton generation rate (i.e. S1+T1)<br/>
output : ratio

`RCC.ConcentrationRatio_SteadyState(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0)`<br/>
a function to calculate the steady state concentration ratio between T1 and S1 (T1/S1).<br/>
output : S1_SS, T1_SS<br/>

`RCC.PulseResponse(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0)`<br/>
a function to calculte the time evolution of S1 and T1 concentration.<br/>
output : S1_t, T1_t, caches<br/>
caches = (ks, kt, ks, kt, S10, T10, V_PF, V_DF)<br/>
V_PF and V_DF are the characteristic vector for PF and DF, respectively.<br/>

#### script
`RCC.script_for_100_PLQY(tau_PF, tau_DF, phi_PF, phi_DF, name='')`<br/>
analyzer for PLQY=100%<br/>

`RCC.script(tau_PF, tau_DF, phi_PF, phi_DF, fpath='', fname='')`<br/>
analyzer for PLQY not equal to 100%<br/>

`RCC.pulseresponse_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, name='')`<br/>
time evolution calculator<br/>

Example:<br/>

    # case: calculate intrinsic rate constants
    # SpiroAC-TRZ
    RCC.script_for_100_PLQY(tau_PF=17e-9, tau_DF=2.1e-6, phi_PF=0.79, phi_DF=0.21, name='SpiroAC-TRZ')
    
    # DPAC-TRZ
    RCC.script(tau_PF=15e-9, tau_DF=2.9e-6, phi_PF=0.70, phi_DF=0.12, fname='DPAC-TRZ', fpath='./data')
    
    # DMAC-TRZ
    RCC.script(tau_PF=20e-9, tau_DF=1.9e-6, phi_PF=0.59, phi_DF=0.31, fname='DMAC-TRZ', fpath='./data')
    
    # case : calculate pulse response
    
    # SpiroAC-TRZ
    t = np.linspace( 0 , 10e-6, 1024 )
    ksr, ksnr, kisc  = 4.66e7, 0, 1.21e7
    ktr, ktnr, krisc =      0, 0, 6.01e5 
    
    name = r'SpiroAC-TRZ $\alpha  = 1.0$'
    RCC.pulseresponse_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1e4, name=name)
    
    name = r'SpiroAC-TRZ $\alpha  = 0.25$'
    RCC.pulseresponse_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=0.25, G=1e4, name=name)

## GUI 

The graphic user interface is built by PyQt5, so please install PyQt5 first.<br/>

    pip install PyQt5

And the GUI can be executed as<br/>

    python GUI_PyQt5_main.py

The panel would be like:
<p align="center">
<img src="./docs/Graph/Panel.png" width="1200">
</p>


### 1. First
select lifetime or rate constant<br/>

### 2. input measured parameters
a) lifetimes (tau)/rate constants(k) of prompt fluorescence(PF) and delayed fluorescence (DF)<br/>
   Please notice the unit. <br/>
b) quantum yield of prompt fluorescence(Phi PF) and delayed fluorescence (Phi DF)<br/>
   or prefactors of prompt fluorescence(PF) and delayed fluorescence (DF) and photoluminescence quantum yield (PLQY)

### 3. calculate
select save file path and name<br/>
<p align="center">
<img src="./docs/Graph/CalculateData.png" width="800">
</p>


### 4. result
Two different cases (all loss from S1 and all loss from T1) are summarized on the terminal.<br/>
<p align="center">
<img src="./docs/Graph/Result.png" width="800">
</p>

The detailed data are in the save file path.<br/>
a) txt<br/>
b) rate constants v.s. loss from T1<br/>

<p align="center">
<img src="./examples/DPAC-TRZ/DPAC-TRZ_rate_constants.png" width="800">
</p>

c) quantum yield v.s. loss from T1
<p align="center">
<img src="./examples/DPAC-TRZ/DPAC-TRZ_quantum_yield.png" width="800">
</p>

## Project Structure
```
RateConstants_TADF/
├── src/                    # Source code
│   ├── gui/               # GUI-related files
│   ├── RateConstantCalculator.py
│   └── RateEquation.py
├── docs/                  # Documentation and images
│   └── Graph/            # Images and diagrams
├── examples/             # Example data and results
│   └── DPAC-TRZ/        # Example calculations
├── tests/               # Test files (future use)
└── README.md           # Project documentation
```