#!/usr/bin/env python3
# Copyright (c) 2020/2025 Wei-Kai Lee. All rights reserved

# coding=utf-8
# -*- coding: utf8 -*-

import numpy as np

def generate_zero_time_series(t: np.ndarray, dtype: np.dtype = np.float_) -> np.ndarray:
    """ 
    generate zeros time series

    Args:
        t: time array

    Returns:
        zero array

    Usage:
        G_S1 = generate_zero_time_series(t)

    """
    return np.zeros( t.size, dtype=dtype )

def solve_finite_difference_rate_equation(
        t_array: np.ndarray, 
        ksr: float, 
        ksnr: float, 
        kisc: float, 
        ktr: float, 
        ktnr: float, 
        krisc: float, 
        S10: float = 1.0, 
        T10: float = 0.0, 
        G_S1: callable = generate_zero_time_series, 
        G_T1: callable = generate_zero_time_series, # initial condition and generation function
        kSSA: float = 0.0, 
        kTTA: float = 0.0, 
        kSTA: float = 0.0,
    ) -> tuple[np.ndarray, np.ndarray]:
    """
    calculate the singlet and triplet exciton concentrations

    Args:
        t_array: time array
        ksr: radition rate constant from S1 to S0
        ksnr: non-radition rate constant from S1 to S0
        kisc: intersystem crossing rate constant from S1 to T1
        ktr: radition rate constant from T1 to S0
        ktnr: non-radition rate constant from T1 to S0
        krisc: reverse intersystem crossing rate constant from T1 to S1
        S10: initial singlet exciton concentration
        T10: initial triplet exciton concentration
        G_S1: generation function of singlet exciton | G_S1(t)
        G_T1: generation function of triplet exciton | G_T1(t)
        kSSA: singlet-singlet annihilation rate constant
        kTTA: triplet-triplet annihilation rate constant
        kSTA: singlet-triplet annihilation rate constant

    Returns:
        S1_t: singlet exciton concentration array
        T1_t: triplet exciton concentration array

    Usage:
        S1_t, T1_t = solve_finite_difference_rate_equation(t, ksr, ksnr, kisc, ktr, ktnr, krisc, S10, T10, G_S1, G_T1, kSSA, kTTA, kSTA)
    """
    
    # total rate constants from S1 and T1
    ks, kt = ksr+ksnr+kisc, ktr+ktnr+krisc

    # initial conditions
    t_array = np.array(t_array, dtype=np.float_)
    S1_t, T1_t = np.zeros( t_array.size, dtype=np.float_), np.zeros( t_array.size, dtype=np.float_)
    S1_t[0], T1_t[0] = S10, T10

    # finite difference
    for ii, t in enumerate(t_array[1:]):
        dt = t[ii] - t[ii-1]

        # Method 1
        S1_t[ii] = S1_t + ( -ks  *S1_t[ii-1] + krisc*T1_t[ii-1] + G_S1(t) - kSSA*(S1_t[ii-1])**2 - kSTA*S1_t[ii-1]*T1_t[ii-1] ) * dt
        T1_t[ii] = T1_t + (  kisc*S1_t[ii-1] - kt   *T1_t[ii-1] + G_T1(t) - kTTA*(T1_t[ii-1])**2 - kSTA*S1_t[ii-1]*T1_t[ii-1] ) * dt

        
    return S1_t, T1_t













