# RateConstants_TADF
Intrinsic rate constants extraction for thermally activated delayed fluorescence (two delayed components)

This module would use numpy and matplotlib, so please install these two modules first.

    pip install numpy and matplotlib

Ref.<br/>
[1] https://onlinelibrary.wiley.com/doi/abs/10.1002/adfm.201602501 <br/>
[2] https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201601675 <br/>
[3] https://pubs.rsc.org/no/content/articlelanding/2015/cc/c5cc05022g/unauth#!divAbstract <br/>

This module provides functions to calculate the intrinsic rate constants between S0 (ground state), S1 (1st excited singlet), and T1 (1st excited triplet state), including S1->S0 (fluorescence and heat loss), S1->T1 (intersystem crossing), T1->S1 (reverse intersystem crossing), and T1->S0 (phosphorescence and heat loss).

Besides, the module also provides the efficiency of TADF material and the exciton concentration.

    import RateConstantsCalculator as RCC
    
### Utilty function
`RCC.tau2k`<br/>
`RCC.k2tau`<br/>
`RCC.exponential_ratio`<br/>
`RCC.phi_PF_DF`<br/>

### intrinsic rate constant calculator
`RCC.IntrinsicRateConstants_Determined`<br/>
`RCC.IntrinsicRateConstants`<br/>

### phi function
`RCC.phi_sr_snr_isc`<br/>
`RCC.phi_tr_tnr_risc`<br/>

### Internal Quantum Efficiency (IQE)
`RCC.IQE_PurcellEffect`<br/>
`RCC.IQE_RateConstants`<br/>
`RCC.IQE_Phi`<br/>
`RCC.IQE_PLQY`<br/>

### Differential Equations
`RCC.PLQY_phi`<br/>
`RCC.k_total`<br/>
`RCC.kPF_kDF`<br/>
`RCC.Concentration_SteadyState`<br/>
`RCC.ConcentrationRatio_SteadyState`<br/>
`RCC.PulseResponse`<br/>

### script
`RCC.script_for_100_PLQY`<br/>
`RCC.script`<br/>
`RCC.pulseresponse_script`<br/>



