#!/usr/bin/env python3
# Copyright (c) 2020 Wei-Kai Lee. All rights reserved

# coding=utf-8
# -*- coding: utf8 -*-
import os
import numpy as np
from matplotlib import pyplot as plt

# RateConstantCalculator is a module to calculate the exciton dynamics of a
#   three-level system (two excited states and single ground state).

# kPF  : rate constant of the prompt fluorescence
# kDF  : rate constant of the delayed fluorescence
# ks   : total rate constant from excited singlet state (S1)
# ksr  : radiative rate constant in S1
# ksnr : non-radiative rate constant in S1
# kisc : intersystem rate constant 
# kt   : total rate constant from excited triplet state (T1)
# ktnr : non-radiative rate constant in T1
# krisc: reversed intersystem rate constant 

# alpha: S1 exciton ratio (1.0 for PL, 0.25 for EL)
# G    : total exciton generation
# F    : Purcell Factor (1.0 for vacuum)


##############################################################################
# untility function
##############################################################################
def tau2k(tau):
    # tau2k is a function to calculate the lifetime from rate constant
    return 1/np.array( tau )
def k2tau(k):
    # k2tau is a function to calculate the rate constant from lifetime
    return 1/np.array( k )
def exponential_ratio(tau_Array, B_Array):
    # exponential_ratio is a function to calculate the integration ratio
    # of each exponetial component. tauArray is a list/numpy array containing
    # lifetime of each exponent and the corresponding strength is in B_Array.
    #   i.e. f(t) = sum( B_Array[i] * exp( -t/tau_Array[i]) )

    tau_Array = np.array( tau_Array )
    B_Array = np.array( B_Array )
    total = np.sum( tau_Array*B_Array )
    return tau_Array*B_Array/total
def phi_PF_DF(PLQY, tauPF, tauDF, B_PF, B_DF):
    # phi_PF_DF is a function to calculte the absolute ratio (considering PLQY)
    #   of the prompt fluorescence and the delayed fluorescence.
    # B_PF and B_DF are the corresponding strength of the exponents.
    #   i.e. f(t) = B_PF * exp( -t/tauPF ) + B_DF * exp( -t/tauDF )

    ratio = exponential_ratio([tauPF, tauDF], [B_PF, B_DF])
    phi_PF = PLQY * ratio[0]
    phi_DF = PLQY * ratio[1]
    return phi_PF, phi_DF
##############################################################################
# intrinsic rate constant calculator
##############################################################################
def IntrinsicRateConstants_Determined(kPF, kDF, phi_PF, phi_DF):
    # IntrinsicRateConstants_Determined is a function to calculate the rate 
    #   constants which can be calculated directly from the experimental data.
    #   (i.e. ksr, kt, kisc*krisc)

    ksr = (kPF*phi_PF + kDF*phi_DF)
    ks  = ( (kPF**2)*phi_PF + (kDF**2)*phi_DF) / ksr
    kt  = ( phi_PF + phi_DF )*kPF*kDF/ksr
    kisckrisc = ( phi_PF*phi_DF*kPF*kDF*((kPF-kDF)**2) )/(ksr**2)

    return ksr, ks, kt, kisckrisc
def IntrinsicRateConstants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL):
    # IntrinsicRateConstants is a function to calculate the rate constants which 
    #   cannot be calculated directly from the experimental data and should considering
    #   a the loss ratio from triplet state by PL excitation. 
    #   

    phi_Tnr_PL_Array = np.array( phi_Tnr_PL, dtype=np.float )

    ksr, ks, kt, kisckrisc = IntrinsicRateConstants_Determined(kPF, kDF, phi_PF, phi_DF)
    
    ksr = ksr * np.ones( phi_Tnr_PL_Array.size, dtype=np.float )
    ks  = ks  * np.ones( phi_Tnr_PL_Array.size, dtype=np.float )
    kt  = kt  * np.ones( phi_Tnr_PL_Array.size, dtype=np.float )

    kisc = (kisckrisc + phi_Tnr_PL_Array*kPF*kDF)/kt
    krisc = kisckrisc/kisc
    ksnr = ks - kisc - ksr

    ktr = np.zeros( phi_Tnr_PL_Array.size, dtype=np.float )
    ktnr = kt - krisc - ktr

    return ks, ksr, ksnr, kisc, kt, ktr, ktnr, krisc
##############################################################################
# phi function
##############################################################################
def phi_sr_snr_isc(ksr, ksnr, kisc):
    # phi_sr_snr_isc is a function to calculate the efficiency of each process
    #   in S1.
    ks = k_total(ksr, ksnr, kisc)
    return ksr/ks, ksnr/ks, kisc/ks
def phi_tr_tnr_risc(ktr, ktnr, krisc):
    # phi_tr_tnr_risc is a function to calculate the efficiency of each process
    #   in T1.
    return phi_sr_snr_isc(ktr, ktnr, krisc)
##############################################################################
# Internal Quantum Efficiency (IQE)
##############################################################################
def IQE_PurcellEffect(IQE, F=1.0):
    return IQE*F/( (1-IQE) + IQE*F )
def IQE_RateConstants(ksr, kt, krisc, kPF, kDF, alpha=1.0, F=1.0):
    # IQE_RateConstants is a function to calculate the internal quantum efficiency
    #   given with intrinsic rate constants.
    IQE = ksr * ( alpha*kt + (1-alpha)*krisc )/( kPF*kDF )
    return IQE_PurcellEffect(IQE, F=F)
def IQE_Phi(phi_sr, phi_isc, phi_risc, alpha=1.0, F=1.0):
    # IQE_RateConstants is a function to calculate the internal quantum efficiency
    #   given with intrinsic state efficiency.
    IQE = ( alpha + (1-alpha)*phi_risc ) * PLQY_phi(phi_sr, phi_isc, phi_risc)
    return IQE_PurcellEffect(IQE, F=F)
def IQE_PLQY(phi_risc, PLQY, alpha=1.0, F=1.0):
    # IQE_RateConstants is a function to calculate the internal quantum efficiency
    #   given with PLQY.
    IQE = ( alpha + (1-alpha)*phi_risc ) * PLQY
    return IQE_PurcellEffect(IQE, F=F)
##############################################################################
# Differential Equations
##############################################################################
def PLQY_phi(phi_sr, phi_isc, phi_risc):
    return phi_sr/( 1-phi_isc*phi_risc )
def k_total(kr, knr, kisc):
    # k_total is a function to calculate the total rate constant
    return kr+knr+kisc
def kPF_kDF(ksr, ksnr, kisc, ktr, ktnr, krisc, F=1.0):
    # kPF_kDF is a function to calculate the prompt and delayed rate constant
    #   from intrinsic rate constants.

    ks, kt = k_total(F*ksr, ksnr, kisc), k_total(ktr, ktnr, krisc)
    B = ks + kt
    Delta = np.sqrt( (ks-kt)**2 + 4*kisc*krisc )
    kPF = ( B + Delta )/2
    kDF = ( B - Delta )/2
    return kPF, kDF
def Concentration_SteadyState(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0):
    # Concentration_SteadyState is a function to calculate the concentration 
    #   of S1 and T1 in steady state, where G is the total exciton generation rate.

    ks, kt = k_total(F*ksr, ksnr, kisc), k_total(ktr, ktnr, krisc)
    kPF, kDF = kPF_kDF(ksr, ksnr, kisc, ktr, ktnr, krisc, F=F)

    S1_SS = ((alpha * kt + (1-alpha)*krisc )/( kPF*kDF )) * G
    T1_SS = ((alpha * kisc + (1-alpha)*ks ) /( kPF*kDF )) * G
    return S1_SS, T1_SS
def ConcentrationRatio_SteadyState(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0):
    # ConcentrationRatio_SteadyState is a function to calculate the concentration ratio
    #   between T1 and S1 in steady state, where G is the total exciton generation rate.

    S1_SS, T1_SS = Concentration_SteadyState(ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=alpha, G=G, F=F)
    return T1_SS/S1_SS
def PulseResponse(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, F=1.0):
    # PulseResponse is a function to caluclate the pulse response of two-level system.

    # total rate constants
    ks, kt = k_total(F*ksr, ksnr, kisc), k_total(ktr, ktnr, krisc)
    kPF, kDF = kPF_kDF(ksr, ksnr, kisc, ktr, ktnr, krisc, F=F)

    # initial concentration
    S10, T10 = alpha*G, (1-alpha)*G

    # eigen vector
    V_PF = np.array( [ (ks-kDF) * S10 -   krisc  * T10,
                        -kisc   * S10 + (kt-kDF) * T10] , dtype=np.float )/(kPF-kDF)
    V_DF = np.array( [ (kPF-ks) * S10 +   krisc  * T10,
                         kisc   * S10 + (kPF-kt) * T10] , dtype=np.float )/(kPF-kDF)
    
    # time series
    S1_t = V_PF[0] * np.exp( -kPF*t ) + V_DF[0] * np.exp( -kDF*t )
    T1_t = V_PF[1] * np.exp( -kPF*t ) + V_DF[1] * np.exp( -kDF*t )

    # caches
    caches = (ks, kt, ks, kt, S10, T10, V_PF, V_DF)

    return S1_t, T1_t, caches
##############################################################################
# script
##############################################################################
def script_for_100_PLQY(tau_PF, tau_DF, phi_PF, phi_DF, name=''):
    kPF, kDF = tau2k([tau_PF, tau_DF])
    # phi_Tnr_PL_Array = np.linspace(0, 1-PLQY, int(200) )
    ks_Array, ksr_Array, ksnr_Array, kisc_Array, kt_Array, ktr_Array, ktnr_Array, krisc_Array \
         = IntrinsicRateConstants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL=0)

    phi_sr_Array, phi_snr_Array, phi_isc_Array = phi_sr_snr_isc(ksr_Array, ksnr_Array, kisc_Array)
    phi_tr_Array, phi_tnr_Array, phi_risc_Array= phi_tr_tnr_risc(ktr_Array, ktnr_Array, krisc_Array)

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
    plt.plot(phi_Tnr_PL_Array*100, krisc_Array, linewidth=2)
    plt.yscale('log')
    plt.xlabel(r'$\Phi_{Tnr}^{PL} [\%]$')
    plt.ylabel(r'rate constants [1/s]')
    plt.title(fname)
    ax.grid(True)
    plt.legend([r'$k_s$', r'$k_{sr}$', r'$k_{snr}$', r'$k_{isc}$', r'$k_t$', r'$k_{risc}$'])
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
    plt.xlabel(r'$\Phi_{Tnr}^{PL} [\%]$')
    plt.ylabel(r'Q.Y. [\%]')
    plt.title(fname)
    ax.grid(True)
    plt.legend([r'$\Phi_{sr}$', r'$\Phi_{snr}$', r'$\Phi_{isc}$', r'$\Phi_{tr}$', r'$\Phi_{tnr}$', r'$\Phi_{risc}$'])
    fig.tight_layout()
    fig.savefig(  os.path.join(fpath,fname + '_quantum_yield') )
    plt.close(fig)
def script(tau_PF, tau_DF, phi_PF, phi_DF, fpath='', fname=''):
    PLQY = phi_PF+phi_DF

    kPF, kDF = tau2k([tau_PF, tau_DF])
    phi_Tnr_PL_Array = np.linspace(0, 1-PLQY, int(200) )
    ks_Array, ksr_Array, ksnr_Array, kisc_Array, kt_Array, ktr_Array, ktnr_Array, krisc_Array \
         = IntrinsicRateConstants(kPF, kDF, phi_PF, phi_DF, phi_Tnr_PL=phi_Tnr_PL_Array)

    phi_sr_Array, phi_snr_Array, phi_isc_Array = phi_sr_snr_isc(ksr_Array, ksnr_Array, kisc_Array)
    phi_tr_Array, phi_tnr_Array, phi_risc_Array= phi_tr_tnr_risc(ktr_Array, ktnr_Array, krisc_Array)

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
def pulseresponse_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1.0, name=''):
    S1_t, T1_t, caches = PulseResponse(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=alpha, G=G)

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
    script(tau_PF=15e-9, tau_DF=2.9e-6, phi_PF=0.70, phi_DF=0.12, fname='DPAC-TRZ', fpath='./data')

    # DMAC-TRZ
    script(tau_PF=20e-9, tau_DF=1.9e-6, phi_PF=0.59, phi_DF=0.31, fname='DMAC-TRZ', fpath='./data')

    # case : calculate pulse response

    # SpiroAC-TRZ
    t = np.linspace( 0 , 10e-6, 1024 )
    ksr, ksnr, kisc  = 4.66e7, 0, 1.21e7
    ktr, ktnr, krisc =      0, 0, 6.01e5 

    name = r'SpiroAC-TRZ $\alpha  = 1.0$'
    pulseresponse_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=1.0, G=1e4, name=name)

    name = r'SpiroAC-TRZ $\alpha  = 0.25$'
    pulseresponse_script(t, ksr, ksnr, kisc, ktr, ktnr, krisc, alpha=0.25, G=1e4, name=name)













