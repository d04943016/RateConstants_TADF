# RateConstants_TADF
Intrinsic rate constants extraction for thermally activated delayed fluorescence (two delayed components)

This module would use numpy and matplotlib, so please install these two modules first.

    pip install numpy and matplotlib

Ref.

[1] https://onlinelibrary.wiley.com/doi/abs/10.1002/adfm.201602501 

[2] https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201601675

[3] https://pubs.rsc.org/no/content/articlelanding/2015/cc/c5cc05022g/unauth#!divAbstract

This module provides functions to calculate the intrinsic rate constants between S0 (ground state), S1 (1st excited singlet), and T1 (1st excited triplet state), including S1->S0 (fluorescence and heat loss), S1->T1 (intersystem crossing), T1->S1 (reverse intersystem crossing), and T1->S0 (phosphorescence and heat loss).

Besides, the module also provides the efficiency of TADF material and the exciton concentration.

### Utilty function
`tau2k`
`k2tau`
`exponential_ratio`
`phi_PF_DF`

### intrinsic rate constant calculator
`IntrinsicRateConstants_Determined`
`IntrinsicRateConstants`

### phi function
`phi_sr_snr_isc`
`phi_tr_tnr_risc`

### Internal Quantum Efficiency (IQE)
`IQE_PurcellEffect`
`IQE_RateConstants`
`IQE_Phi`
`IQE_PLQY`

### Differential Equations
`PLQY_phi`
`k_total`
`kPF_kDF`
`Concentration_SteadyState`
`ConcentrationRatio_SteadyState`
`PulseResponse`

### script
`script_for_100_PLQY`
`script`
`pulseresponse_script`



