# RateConstants_TADF
Intrinsic rate constants extraction for thermally activated delayed fluorescence (two delayed components)

This module would use numpy and matplotlib, so please install these two modules first.

    pip install numpy and matplotlib

Ref.<br/>
[1] https://onlinelibrary.wiley.com/doi/abs/10.1002/adfm.201602501 <br/>
[2] https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201601675 <br/>
[3] https://pubs.rsc.org/no/content/articlelanding/2015/cc/c5cc05022g/unauth#!divAbstract <br/>
[4] https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201704961 <br/>

This module provides functions to calculate the intrinsic rate constants between S0 (ground state), S1 (1st excited singlet), and T1 (1st excited triplet state), including S1->S0 (fluorescence and heat loss), S1->T1 (intersystem crossing), T1->S1 (reverse intersystem crossing), and T1->S0 (phosphorescence and heat loss).

Besides, the module also provides the efficiency of TADF material and the exciton concentration.

## RateConstantsCalculator Module

    import RateConstantsCalculator as RCC
    
#### Utilty function
`RCC.tau2k(tau)`<br/>
a function to calculate rate constant (k) from lifetime (tau).

`RCC.k2tau(k)`<br/>
a function to calculate lifetime from (tau) rate constant (k).

`RCC.exponential_ratio(tau_Array, B_Array)`<br/>
a function to calculate the relative ratio of each exponent

`RCC.phi_PF_DF(PLQY, tauPF, tauDF, B_PF, B_DF)`<br/>
a function to calculate the quantum efficiency contributed from prompt fluorescence and the delayed fluorescence

#### intrinsic rate constant calculator
`RCC.IntrinsicRateConstants_Determined(kPF, kDF, phi_PF, phi_DF)`<br/>
a function to calculate the intrinsice rate constants that can be determined by kPF (rate constant of prompt fluorescence), kDF (rate constant of delayed fluorescence), phi_PF (prompt fluorescence quantum yield) and phi_DF (delayed fluorescence quantum yield). <br/>
These four values (kPF, kDF, phi_PF, phi_DF) can be directly extracted from the transient data 

    kPF=RCC.tau2k(tauPF)
    kDF=RCC.tau2k(tauDF) 
    
and photoluminescence quantum yield (PLQY, phi_PF, phi_DF).

    [phi_PF, phi_DF] = RCC.phi_PF_DF(PLQY, tauPF, tauDF, B_PF, B_DF)
    ksr, ks, kt, kisckrisc = RCC.IntrinsicRateConstants_Determined(kPF, kDF, phi_PF, phi_DF)
    
output : ksr (radiative rate constant from singlet state), ks (total rate constant from singlet state), kt (total rate constant from triplet state), kisckrisc (a product of the intersystem rate constant and the reverse intersystem crossing rate constant)

`RCC.IntrinsicRateConstants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL)`<br/>


#### phi function
`RCC.phi_sr_snr_isc(ksr, ksnr, kisc)`<br/>
`RCC.phi_tr_tnr_risc(ktr, ktnr, krisc)`<br/>

#### Internal Quantum Efficiency (IQE)
`RCC.IQE_PurcellEffect(IQE, F=1.0)`<br/>
`RCC.IQE_RateConstants(ksr, kt, krisc, kPF, kDF, alpha=1.0, F=1.0)`<br/>
`RCC.IQE_Phi(phi_sr, phi_isc, phi_risc, alpha=1.0, F=1.0)`<br/>
`RCC.IQE_PLQY(phi_risc, PLQY, alpha=1.0, F=1.0)`<br/>

#### Differential Equations
`RCC.PLQY_phi(phi_sr, phi_isc, phi_risc)`<br/>
`RCC.k_total(kr, knr, kisc)`<br/>
`RCC.kPF_kDF(ksr, ksnr, kisc, ktr, ktnr, krisc, F=1.0)`<br/>
`RCC.Concentration_SteadyState(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0)`<br/>
`RCC.ConcentrationRatio_SteadyState(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0)`<br/>
`RCC.PulseResponse(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0)`<br/>

#### script
`RCC.script_for_100_PLQY(tau_PF, tau_DF, phi_PF, phi_DF, name='')`<br/>
`RCC.script(tau_PF, tau_DF, phi_PF, phi_DF, name='')`<br/>
`RCC.pulseresponse_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, name='')`<br/>



