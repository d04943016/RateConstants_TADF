#!/usr/bin/env python3
# Copyright (c) 2020/2025 Wei-Kai Lee. All rights reserved

# coding=utf-8
# -*- coding: utf8 -*-
import os
import numpy as np
from matplotlib import pyplot as plt

"""
RateConstantCalculator is a module for computing the exciton dynamics in a three-level
photophysical system, which includes two excited states (singlet S1 and triplet T1) and one ground state S0.
This module provides functions for determining various rate constants and efficiencies
associated with prompt fluorescence, delayed fluorescence, intersystem crossing, and non-radiative decay.

Key rate constants:
    kPF  : Rate constant of the prompt fluorescence emission.
    kDF  : Rate constant of the delayed fluorescence emission.
    ks   : Total rate constant from the excited singlet state (S1).
    ksr  : Radiative rate constant for the singlet state (S1).
    ksnr : Non-radiative rate constant for the singlet state (S1).
    kisc : Intersystem crossing rate constant from S1 to T1.
    kt   : Total rate constant from the excited triplet state (T1).
    ktnr : Non-radiative rate constant for the triplet state (T1).
    krisc: Reverse intersystem crossing rate constant from T1 back to S1.

Additional parameters:
    alpha: Fraction of excitons created in the singlet state (typically 1.0 for photoluminescence (PL)
           and 0.25 for electroluminescence (EL)).
    G    : Total exciton generation rate.
    F    : Purcell factor, representing the enhancement of radiative decay (F = 1.0 in vacuum conditions).

"""

""" untility function """
def tau2k(tau: np.ndarray) -> np.ndarray:
    """ 
    tau2k is a function to calculate the lifetime from rate constant

    Args:
        tau: lifetime

    Returns:
        rate constant
    
    Usage:
        k = tau2k([tau_PF, tau_DF])

    """

    return 1/np.asarray( tau )

def k2tau(k: np.ndarray) -> np.ndarray:
    """ 
    k2tau is a function to calculate the rate constant from lifetime 

    Args:
        k: rate constant

    Returns:
        lifetime
    
    Usage:
        tau = k2tau(k)
    
    """
    return 1/np.asarray( k )

def cal_exponential_ratio(tau_Array: np.ndarray, B_Array: np.ndarray) -> np.ndarray:
    """ 
    calculate the integration ratio of each exponetial component. 
    
    tau_Array contains lifetimes of exponential components, with B_Array containing their corresponding strengths.
    
    i.e. f(t) = sum( B_Array[i] * exp( -t/tau_Array[i]) )

    """

    # convert to numpy array
    tau_Array = np.asarray( tau_Array )
    B_Array = np.asarray( B_Array ) 

    # calculate the total
    total = np.sum( tau_Array*B_Array )

    # return the integration ratio
    return tau_Array*B_Array/total

def cal_phi_PF_DF(PLQY: float, tauPF: float, tauDF: float, B_PF: float, B_DF: float) -> tuple[float, float]:
    """ 
    calculate the absolute ratio (considering PLQY) of the prompt fluorescence and the delayed fluorescence.
    
    B_PF and B_DF are the corresponding strength of the exponents.
    #   i.e. f(t) = B_PF * exp( -t/tauPF ) + B_DF * exp( -t/tauDF )

    Args:
        PLQY: PLQY / photoluminescence quantum yield
        tauPF: lifetime of the prompt fluorescence
        tauDF: lifetime of the delayed fluorescence
        B_PF: strength of the prompt fluorescence
        B_DF: strength of the delayed fluorescence

    Returns:
        phi_PF: absolute ratio of the prompt fluorescence
        phi_DF: absolute ratio of the delayed fluorescence

    Usage:
        phi_PF, phi_DF = cal_phi_PF_DF(PLQY, tauPF, tauDF, B_PF, B_DF)
        
    """

    ratio = cal_exponential_ratio([tauPF, tauDF], [B_PF, B_DF])
    phi_PF = PLQY * ratio[0]
    phi_DF = PLQY * ratio[1]
    return phi_PF, phi_DF

""" intrinsic rate constant calculator """ 
def cal_determined_intrinsic_rate_constants(kPF:float, kDF:float, phi_PF:float, phi_DF:float) -> tuple[float, float, float, float]:
    """ 
    calculate the rate constants which can be calculated directly from the experimental data.
    
    (i.e. ksr, kt, kisc*krisc)

    Args:
        kPF: prompt fluorescence rate constant
        kDF: delayed fluorescence rate constant
        phi_PF: prompt fluorescence ratio
        phi_DF: delayed fluorescence ratio

    Returns:
        ksr: radiative rate constant in S1
        ks: total rate constant from excited singlet state (S1)
        kt: total rate constant from excited triplet state (T1)
        kisckrisc: the product of intersystem crossing rate constant and reversed intersystem crossing rate constant

    """

    ksr = (kPF*phi_PF + kDF*phi_DF)
    ks  = ( (kPF**2)*phi_PF + (kDF**2)*phi_DF) / ksr
    kt  = ( phi_PF + phi_DF )*kPF*kDF/ksr
    kisckrisc = ( phi_PF*phi_DF*kPF*kDF*((kPF-kDF)**2) )/(ksr**2)

    return ksr, ks, kt, kisckrisc

def cal_intrinsic_rate_constants(kPF:float, kDF:float, phi_PF:float, phi_DF:float, phi_Tnr_PL:np.ndarray) -> tuple[float, float, float, float, float, float, float, float]:
    """ 
    calculate the rate constants which cannot be calculated directly from the experimental data and should considering
    a the loss ratio from triplet state by PL excitation. 

    Args:
        kPF: prompt fluorescence rate constant
        kDF: delayed fluorescence rate constant
        phi_PF: prompt fluorescence ratio
        phi_DF: delayed fluorescence ratio
        phi_Tnr_PL: the loss ratio from triplet state by PL excitation [array]

    Returns:
        ks: total rate constant from excited singlet state (S1)
        ksr: radiative rate constant in S1
        ksnr: non-radiative rate constant in S1
        kisc: intersystem crossing rate constant
        kt: total rate constant from excited triplet state (T1)
        ktr: radiative rate constant in T1
        ktnr: non-radiative rate constant in T1

    """

    # convert to numpy array
    phi_Tnr_PL_Array = np.asarray( phi_Tnr_PL, dtype=np.float_ )

    # calculate the determined rate constants
    ksr, ks, kt, kisckrisc = cal_determined_intrinsic_rate_constants(kPF, kDF, phi_PF, phi_DF)
    
    ksr = ksr * np.ones( phi_Tnr_PL_Array.size, dtype=np.float_ )
    ks  = ks  * np.ones( phi_Tnr_PL_Array.size, dtype=np.float_ )
    kt  = kt  * np.ones( phi_Tnr_PL_Array.size, dtype=np.float_ )

    # calculate the non-determined rate constants
    kisc = (kisckrisc + phi_Tnr_PL_Array*kPF*kDF)/kt
    krisc = kisckrisc/kisc
    ksnr = ks - kisc - ksr

    ktr = np.zeros( phi_Tnr_PL_Array.size, dtype=np.float_ )
    ktnr = kt - krisc - ktr

    # return the rate constants
    return ks, ksr, ksnr, kisc, kt, ktr, ktnr, krisc

""" phi/efficiency function"""
def cal_phi_sr_snr_isc(ksr:float, ksnr:float, kisc:float) -> tuple[float, float, float]:
    """
    calculate the ratio of each process in S1.

    Args:
        ksr: radiative rate constant in S1
        ksnr: non-radiative rate constant in S1
        kisc: intersystem crossing rate constant

    Returns:
        phi_sr: radiative ratio in S1
        phi_snr: non-radiative ratio in S1
        phi_isc: intersystem crossing ratio in S1

    Usage:
        phi_sr, phi_snr, phi_isc = cal_phi_sr_snr_isc(ksr, ksnr, kisc)

    """

    ks = cal_k_total(ksr, ksnr, kisc)
    return ksr/ks, ksnr/ks, kisc/ks

def cal_phi_tr_tnr_risc(ktr:float, ktnr:float, krisc:float) -> tuple[float, float, float]:
    """
    calculate the ratio of each process in T1.

    Args:
        ktr: radiative rate constant in T1
        ktnr: non-radiative rate constant in T1
        krisc: reversed intersystem crossing rate constant
    
    Returns:
        phi_tr: radiative ratio in T1
        phi_tnr: non-radiative ratio in T1
        phi_risc: reversed intersystem crossing ratio in T1

    Usage:
        phi_tr, phi_tnr, phi_risc = cal_phi_tr_tnr_risc(ktr, ktnr, krisc)

    """
    return cal_phi_sr_snr_isc(ktr, ktnr, krisc)

""" Internal Quantum Efficiency (IQE) """
def cal_IQE_enhanced_by_Purcell_effect(IQE:float, F:float=1.0) -> float:
    """
    calculate the internal quantum efficiency (IQE) considering the Purcell effect.

    Args:
        IQE: internal quantum efficiency 
        F: Purcell factor

    """
    return IQE*F/( (1-IQE) + IQE*F )

def cal_IQE_from_rate_constants(ksr:float, kt:float, krisc:float, kPF:float, kDF:float, alpha:float=1.0, F:float=1.0) -> float:
    """ 
    calculate the internal quantum efficiency (IQE) by rate constants.

    Args:
        ksr: radiative rate constant in S1
        kt: total rate constant from excited triplet state (T1)
        krisc: reversed intersystem crossing rate constant
        kPF: prompt fluorescence rate constant
        kDF: delayed fluorescence rate constant
        alpha: the ratio of excited singlet state
        F: Purcell factor

    Returns:
        IQE: internal quantum efficiency

    Usage:
        IQE = cal_IQE_from_rate_constants(ksr, kt, krisc, kPF, kDF, alpha, F)

    """

    IQE = ksr * ( alpha*kt + (1-alpha)*krisc )/( kPF*kDF )
    return cal_IQE_enhanced_by_Purcell_effect(IQE, F=F)

def cal_IQE_from_phi(phi_sr:float, phi_isc:float, phi_risc:float, alpha:float=1.0, F:float=1.0) -> float:
    """
    calculate the internal quantum efficiency (IQE) by the ratio of each process.

    Args:
        phi_sr: radiative ratio in S1
        phi_isc: intersystem crossing ratio in S1
        phi_risc: reversed intersystem crossing ratio in S1
        alpha: the ratio of excited singlet state
        F: Purcell factor

    Returns:
        IQE: internal quantum efficiency

    Usage:
        IQE = cal_IQE_from_phi(phi_sr, phi_isc, phi_risc, alpha, F)
    
    """
    PLQY = cal_PLQY_from_phi(phi_sr, phi_isc, phi_risc)
    return cal_IQE_from_PLQY(phi_risc = phi_risc, PLQY = PLQY, alpha=alpha, F=F)

def cal_IQE_from_PLQY(phi_risc:float, PLQY:float, alpha:float=1.0, F:float=1.0) -> float:
    """
    calculate the internal quantum efficiency (IQE) considering the PLQY.

    Args:
        phi_risc: reversed intersystem crossing ratio in T1
        PLQY: photoluminescence quantum yield
        alpha: the ratio of excited singlet state
        F: Purcell factor

    Returns:
        IQE: internal quantum efficiency

    Usage:
        IQE = cal_IQE_from_PLQY(phi_risc, PLQY, alpha, F)

    """

    IQE = ( alpha + (1-alpha)*phi_risc ) * PLQY
    return cal_IQE_enhanced_by_Purcell_effect(IQE, F=F)

def cal_PLQY_from_phi(phi_sr:float, phi_isc:float, phi_risc:float) -> float:
    """
    calculate the PLQY considering the ratio of each process.

    Args:
        phi_sr: radiative ratio in S1
        phi_isc: intersystem crossing ratio in S1
        phi_risc: reversed intersystem crossing ratio in T1

    Returns:
        PLQY: photoluminescence quantum yield

    Usage:
        PLQY = cal_PLQY_from_phi(phi_sr, phi_isc, phi_risc)

    """

    return phi_sr/( 1-phi_isc*phi_risc )

""" Differential Equations """
def cal_k_total(kr:float, knr:float, kisc:float) -> float:
    """ 
    calculate the total rate constant.
   
    Args:
        kr: radiative rate constant
        knr: non-radiative rate constant
        kisc: intersystem crossing rate constant

    Returns:
        k_total: total rate constant

    Usage:
        k_total = cal_k_total(kr, knr, kisc)

    """
    return kr+knr+kisc

def cal_kPF_kDF(ksr:float, ksnr:float, kisc:float, ktr:float, ktnr:float, krisc:float, F:float=1.0) -> tuple[float, float]:
    """
    calculate the prompt and delayed rate constant from intrinsic rate constants.

    Args:
        ksr: radiative rate constant in S1
        ksnr: non-radiative rate constant in S1
        kisc: intersystem crossing rate constant
        ktr: radiative rate constant in T1
        ktnr: non-radiative rate constant in T1
        krisc: reversed intersystem crossing rate constant
        F: Purcell factor

    Returns:
        kPF: prompt fluorescence rate constant
        kDF: delayed fluorescence rate constant

    Usage:
        kPF, kDF = cal_kPF_kDF(ksr, ksnr, kisc, ktr, ktnr, krisc, F)
        
    """

    ks, kt = cal_k_total(F*ksr, ksnr, kisc), cal_k_total(ktr, ktnr, krisc)
    B = ks + kt
    Delta = np.sqrt( (ks-kt)**2 + 4*kisc*krisc )
    kPF = ( B + Delta )/2
    kDF = ( B - Delta )/2
    return kPF, kDF

def cal_steady_state_exciton_concentration(ksr:float, ksnr:float, kisc:float, ktr:float, ktnr:float, krisc:float, alpha:float=1.0, G:float=1.0, F:float=1.0):
    """
    calcualte the concentration of S1 and T1 in steady state.

    Args:
        ksr: radiative rate constant in S1
        ksnr: non-radiative rate constant in S1
        kisc: intersystem crossing rate constant
        ktr: radiative rate constant in T1
        ktnr: non-radiative rate constant in T1
        krisc: reversed intersystem crossing rate constant
        alpha: the ratio of excited singlet state
        G: total exciton generation rate
        F: Purcell factor
    
    Returns:
        S1_SS: concentration of S1 in steady state
        T1_SS: concentration of T1 in steady state
    
    Usage:
        S1_SS, T1_SS = cal_steady_state_exciton_concentration(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha, G, F)

    """

    ks, kt = cal_k_total(F*ksr, ksnr, kisc), cal_k_total(ktr, ktnr, krisc)
    kPF, kDF = cal_kPF_kDF(ksr, ksnr, kisc, ktr, ktnr, krisc, F=F)

    S1_SS = ((alpha * kt + (1-alpha)*krisc )/( kPF*kDF )) * G
    T1_SS = ((alpha * kisc + (1-alpha)*ks ) /( kPF*kDF )) * G
    return S1_SS, T1_SS

def cal_steady_state_exciton_concetration_ratio(ksr:float, ksnr:float, kisc:float, ktr:float, ktnr:float, krisc:float, alpha:float=1.0, G:float=1.0, F:float=1.0):
    """
    calculate the concentration ratio between T1 and S1 in steady state.

    Args:
        ksr: radiative rate constant in S1
        ksnr: non-radiative rate constant in S1
        kisc: intersystem crossing rate constant
        ktr: radiative rate constant in T1
        ktnr: non-radiative rate constant in T1
        krisc: reversed intersystem crossing rate constant
        alpha: the ratio of excited singlet state
        G: total exciton generation rate
        F: Purcell factor
    
    Returns:
        ratio: concentration ratio between T1 and S1 in steady state

    Usage:
        ratio = cal_steady_state_exciton_concetration_ratio(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha, G, F)

    """
    S1_SS, T1_SS = cal_steady_state_exciton_concentration(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=alpha, G=G, F=F)
    return T1_SS/S1_SS

def cal_pulse_response(t:np.ndarray, ksr:float, ksnr:float, kisc:float, ktr:float, ktnr:float, krisc:float, alpha:float=1.0, G:float=1.0, F:float=1.0):
    """
    calculate the pulse response of three-level system.

    Args:
        t: time series
        ksr: radiative rate constant in S1
        ksnr: non-radiative rate constant in S1
        kisc: intersystem crossing rate constant
        ktr: radiative rate constant in T1
        ktnr: non-radiative rate constant in T1
        krisc: reversed intersystem crossing rate constant
        alpha: the ratio of excited singlet state
        G: total exciton generation rate
        F: Purcell factor
    
    Returns:
        S1_t: concentration of S1 in time series
        T1_t: concentration of T1 in time series

    """
    # total rate constants
    ks, kt = cal_k_total(F*ksr, ksnr, kisc), cal_k_total(ktr, ktnr, krisc)
    kPF, kDF = cal_kPF_kDF(ksr, ksnr, kisc, ktr, ktnr, krisc, F=F)

    # initial concentration
    S10, T10 = alpha*G, (1-alpha)*G

    # eigen vector
    V_PF = np.array( [ (ks-kDF) * S10 -   krisc  * T10,
                        -kisc   * S10 + (kt-kDF) * T10] , dtype=np.float_ )/(kPF-kDF)
    V_DF = np.array( [ (kPF-ks) * S10 +   krisc  * T10,
                         kisc   * S10 + (kPF-kt) * T10] , dtype=np.float_ )/(kPF-kDF)
    
    # time series
    S1_t = V_PF[0] * np.exp( -kPF*t ) + V_DF[0] * np.exp( -kDF*t )
    T1_t = V_PF[1] * np.exp( -kPF*t ) + V_DF[1] * np.exp( -kDF*t )

    # caches
    caches = (ks, kt, ks, kt, S10, T10, V_PF, V_DF)

    return S1_t, T1_t, caches

""" script """
def script_for_100_PLQY(tau_PF, tau_DF, phi_PF, phi_DF, name=''):
    kPF, kDF = tau2k([tau_PF, tau_DF])
    # phi_Tnr_PL_Array = np.linspace(0, 1-PLQY, int(200) )
    ks_Array, ksr_Array, ksnr_Array, kisc_Array, kt_Array, ktr_Array, ktnr_Array, krisc_Array \
         = cal_intrinsic_rate_constants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL=0)

    phi_sr_Array, phi_snr_Array, phi_isc_Array = cal_phi_sr_snr_isc(ksr_Array, ksnr_Array, kisc_Array)
    phi_tr_Array, phi_tnr_Array, phi_risc_Array= cal_phi_tr_tnr_risc(ktr_Array, ktnr_Array, krisc_Array)

    print('=============================== {0} ==============================='.format(name) )
    print('* Rate Constants [1/s] :')
    print('  kPF= {0:>5.2e}, kDF = {1:>5.2e}'.format(kPF, kDF))
    print('  ks = {0:>5.2e}, ksr = {1:>5.2e}, ksnr = {2:>5.2e}, kisc  = {3:>5.2e}'.format(ks_Array[0], ksr_Array[0], ksnr_Array[0], kisc_Array[0]) )
    print('  kt = {0:>5.2e}, ktr = {1:>5.2e}, ktnr = {2:>5.2e}, krisc = {3:>5.2e}'.format(kt_Array[0], ktr_Array[0], ktnr_Array[0], krisc_Array[0]) )
    print('* Process Efficiency [%] :')
    print('  phi_sr = {0:>6.2f}, phi_snr = {1:>6.2f}, phi_isc = {2:>6.2f}'.format(phi_sr_Array[0]*100, phi_snr_Array[0]*100, phi_isc_Array[0]*100) )
    print('  phi_tr = {0:>6.2f}, phi_tnr = {1:>6.2f}, phi_risc= {2:>6.2f}'.format(phi_tr_Array[0]*100, phi_tnr_Array[0]*100, phi_risc_Array[0]*100) )
    print('')

def save_data(fpath, fname, phi_Tnr_PL_Array, ks_Array, ksr_Array, ksnr_Array, kisc_Array, kt_Array, ktr_Array, ktnr_Array, krisc_Array, phi_sr_Array, phi_snr_Array, phi_isc_Array, phi_tr_Array, phi_tnr_Array, phi_risc_Array):
    # write
    if not os.path.isdir(fpath):
        os.makedirs(fpath)
    with open( os.path.join(fpath, fname+'.txt'), 'w') as file:
        file.write('{0:>15s} {1:>15s} {2:>15s} {3:>15s} {4:>15s} {5:>15s} {6:>15s} {7:>15s} {8:>15s} {9:>15s} {10:>15s} {11:>15s} {12:>15s} {13:>15s}\n'.format('ks(1/s)', 'ksr(1/s)', 'ksnr(1/s)', 'kisc(1/s)', 'kt(1/s)', 'ktr(1/s)', 'ktnr(1/s)', 'krisc(1/s)', 'phi_sr(%)', 'phi_snr(%)', 'phi_isc(%)', 'phi_tr(%)', 'phi_tnr(%)', 'phi_risc(%)') )
        for ii in range( len(phi_Tnr_PL_Array) ):
            file.write('{0:>15.5e} {1:>15.5e} {2:>15.5e} {3:>15.5e} {4:>15.5e} {5:>15.5e} {6:>15.5e} {7:>15.5e} {8:>15.5f} {9:>15.5f} {10:>15.5f} {11:>15.5f} {12:>15.5f} {13:>15.5f}\n'.format(ks_Array[ii], ksr_Array[ii], ksnr_Array[ii], kisc_Array[ii], kt_Array[ii], ktr_Array[ii], ktnr_Array[ii], krisc_Array[ii], phi_sr_Array[ii]*100, phi_snr_Array[ii]*100, phi_isc_Array[ii]*100, phi_tr_Array[ii]*100, phi_tnr_Array[ii]*100, phi_risc_Array[ii]*100) )

def plot_data(fpath, fname, phi_Tnr_PL_Array, ks_Array, ksr_Array, ksnr_Array, kisc_Array, kt_Array, ktr_Array, ktnr_Array, krisc_Array, phi_sr_Array, phi_snr_Array, phi_isc_Array, phi_tr_Array, phi_tnr_Array, phi_risc_Array):
    if not os.path.isdir(fpath):
        os.makedirs(fpath)

    # rate constants    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(phi_Tnr_PL_Array*100, ks_Array   , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, ksr_Array  , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, ksnr_Array , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, kisc_Array , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, kt_Array   , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, ktnr_Array , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, krisc_Array, linewidth=2)
    plt.yscale('log')
    plt.xlabel(r'$\Phi_{Tnr}^{PL} [\%]$', fontsize=20)
    plt.ylabel(r'rate constants [1/s]', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(fname, fontsize=20)
    ax.grid(True)
    plt.legend([r'$k_s$', r'$k_{sr}$', r'$k_{snr}$', r'$k_{isc}$', r'$k_t$', r'$k_{tnr}$', r'$k_{risc}$'], fontsize=15)
    ax.set_ylim(1e4,1e8)
    fig.tight_layout()
    fig.savefig(  os.path.join(fpath,fname + '_rate_constants') )
    plt.close(fig)

    # qauntum yield
    # rate constants    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(phi_Tnr_PL_Array*100, phi_sr_Array*100   , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, phi_snr_Array*100  , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, phi_isc_Array*100  , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, phi_tr_Array*100   , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, phi_tnr_Array*100  , linewidth=2)
    plt.plot(phi_Tnr_PL_Array*100, phi_risc_Array*100 , linewidth=2)
    plt.xlabel(r'$\Phi_{Tnr}^{PL} [\%]$', fontsize=20)
    plt.ylabel(r'Q.Y. [\%]', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(fname, fontsize=20)
    ax.grid(True)
    plt.legend([r'$\Phi_{sr}$', r'$\Phi_{snr}$', r'$\Phi_{isc}$', r'$\Phi_{tr}$', r'$\Phi_{tnr}$', r'$\Phi_{risc}$'], fontsize=15)
    fig.tight_layout()
    fig.savefig(  os.path.join(fpath,fname + '_quantum_yield') )
    plt.close(fig)

def script(tau_PF, tau_DF, phi_PF, phi_DF, fpath='', fname=''):
    PLQY = phi_PF+phi_DF

    kPF, kDF = tau2k([tau_PF, tau_DF])
    phi_Tnr_PL_Array = np.linspace(0, 1-PLQY, int(200) )
    ks_Array, ksr_Array, ksnr_Array, kisc_Array, kt_Array, ktr_Array, ktnr_Array, krisc_Array \
         = cal_intrinsic_rate_constants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL=phi_Tnr_PL_Array)

    phi_sr_Array, phi_snr_Array, phi_isc_Array = cal_phi_sr_snr_isc(ksr_Array, ksnr_Array, kisc_Array)
    phi_tr_Array, phi_tnr_Array, phi_risc_Array= cal_phi_tr_tnr_risc(ktr_Array, ktnr_Array, krisc_Array)

    print('=============================== {0} ==============================='.format(fname) )
    print('* Rate Constants [1/s] :')
    print('  kPF= {0:>5.2e}, kDF = {1:>5.2e}'.format(kPF, kDF))
    for idx in [0,-1]:
        print('** phi_Tnr_PL = {0:>6.2f}%'.format(phi_Tnr_PL_Array[idx]*100) )
        print('  ks = {0:>5.2e}, ksr = {1:>5.2e}, ksnr = {2:>5.2e}, kisc  = {3:>5.2e}'.format(ks_Array[idx], ksr_Array[idx], ksnr_Array[idx], kisc_Array[idx]) )
        print('  kt = {0:>5.2e}, ktr = {1:>5.2e}, ktnr = {2:>5.2e}, krisc = {3:>5.2e}'.format(kt_Array[idx], ktr_Array[idx], ktnr_Array[idx], krisc_Array[idx]) )
    print('')
    print('* Process Efficiency [%] :')
    for idx in [0,-1]:
        print('** phi_Tnr_PL = {0:>6.2f}%'.format(phi_Tnr_PL_Array[idx]*100) )
        print('  phi_sr = {0:>6.2f}, phi_snr = {1:>6.2f}, phi_isc = {2:>6.2f}'.format(phi_sr_Array[idx]*100, phi_snr_Array[idx]*100, phi_isc_Array[idx]*100) )
        print('  phi_tr = {0:>6.2f}, phi_tnr = {1:>6.2f}, phi_risc= {2:>6.2f}'.format(phi_tr_Array[idx]*100, phi_tnr_Array[idx]*100, phi_risc_Array[idx]*100) )
    print('')

    # save data 
    if (fpath!='') and (fname!=''):
        save_data(fpath, fname, phi_Tnr_PL_Array, 
                  ks_Array, ksr_Array, ksnr_Array, kisc_Array, 
                  kt_Array, ktr_Array, ktnr_Array, krisc_Array, 
                  phi_sr_Array, phi_snr_Array, phi_isc_Array, 
                  phi_tr_Array, phi_tnr_Array, phi_risc_Array)
        
        plot_data(fpath, fname, phi_Tnr_PL_Array,
                  ks_Array, ksr_Array, ksnr_Array, kisc_Array, 
                  kt_Array, ktr_Array, ktnr_Array, krisc_Array, 
                  phi_sr_Array, phi_snr_Array, phi_isc_Array, 
                  phi_tr_Array, phi_tnr_Array, phi_risc_Array)
        
def pulse_response_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, name=''):
    S1_t, T1_t, caches = cal_pulse_response(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=alpha, G=G)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(t*1e6, S1_t   , linewidth=2)
    plt.plot(t*1e6, T1_t   , linewidth=2)
    plt.yscale('log')
    plt.xlabel(r'time $[\mu s]$')
    plt.ylabel(r'$[S_1], [T_1]$')
    plt.title(name)
    ax.grid(True)
    plt.xlim([np.min(t), np.max(t)*1e6])
    plt.ylim([G*1e-5, G*10])

    plt.legend([r'$[S_1](t)$', r'$[T_1](t)$'])
    plt.show()

if __name__ == '__main__':

    # case: calculate intrinsic rate constants
    # SpiroAC-TRZ
    script_for_100_PLQY(tau_PF=17e-9, tau_DF=2.1e-6, phi_PF=0.79, phi_DF=0.21, name='SpiroAC-TRZ')

    # DPAC-TRZ
    script(tau_PF=15e-9, tau_DF=2.9e-6, phi_PF=0.70, phi_DF=0.12, fname='DPAC-TRZ', fpath='../data')

    # DMAC-TRZ
    script(tau_PF=20e-9, tau_DF=1.9e-6, phi_PF=0.59, phi_DF=0.31, fname='DMAC-TRZ', fpath='../data')

    # case : calculate pulse response

    # SpiroAC-TRZ
    t = np.linspace( 0 , 10e-6, 1024 )
    ksr, ksnr, kisc  = 4.66e7, 0, 1.21e7
    ktr, ktnr, krisc =      0, 0, 6.01e5 

    name = r'SpiroAC-TRZ $\alpha  = 1.0$'
    pulse_response_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1e4, name=name)

    name = r'SpiroAC-TRZ $\alpha  = 0.25$'
    pulse_response_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=0.25, G=1e4, name=name)













