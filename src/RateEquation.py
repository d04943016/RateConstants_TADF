#!/usr/bin/env python3
# Copyright (c) 2020 Wei-Kai Lee. All rights reserved

# coding=utf-8
# -*- coding: utf8 -*-

import numpy as np

def zero_generation(t):
    return np.zeros( t.size, dtype=np.float )
def time_series(t_start, t_end, dt = 1e-9):
    t_num = np.ceil( (t_start-t_end)/dt )
    return np.linspace( t_start, t_end, t_num )
def FiniteDifference(t_array, ksr, ksnr, kisc, ktr, ktnr, krisc, # first order rate constants
                     S10 = 1.0, T10 = 0.0, G_S1 = zero_generation, G_T1 = zero_generation, # initial condition and generation function
                     kSTA = 0.0, kTTA = 0.0, kSTA = 0.0): # higer order rate constants
    
    # total rate constants from S1 and T1
    ks, kt = ksr+ksnr+kisc, ktr+ktnr+krisc

    # initial conditions
    t_array = np.array(t_array, dtype=np.float)
    S1_t, T1_t = np.zeros( t_array.size, dtype=np.float), np.zeros( t_array.size, dtype=np.float)
    S1_t[0], T1_t[0] = S10, T10

    # finite difference
    for ii, t in enumerate(t_array[1:]):
        dt = t[ii] - t[ii-1]

        # Method 1
        S1_t[ii] = S1_t + ( -ks  *S1_t[ii-1] + krisc*T1_t[ii-1] + G_S1(t) - kSSA*(S1_t[ii-1])**2 - kSTA*S1_t[ii-1]*T1_t[ii-1] ) * dt
        T1_t[ii] = T1_t + (  kisc*S1_t[ii-1] - kt   *T1_t[ii-1] + G_T1(t) - kTTA*(T1_t[ii-1])**2 - kSTA*S1_t[ii-1]*T1_t[ii-1] ) * dt

        # Method 2
        
        
    return S1_t, T1_t













