#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 08:47:59 2017

@author: Hyeon-Jae
"""

from scipy.integrate import odeint
import numpy as np

# define constants

# transcription / mRNA degradation
k_t
k_td = exp(r_t * t)
r_t
# protein production/degradation
k_A
k_Ad = exp(r_A * t) # not sure if this is necessary
r_A
k_C
k_Cd = exp(r_C * t)
r_C
k_EG
k_EGd = exp(r_EG * t)
r_EG
# secretion 
k_IM
k_OM
# polymerization 
k_nu_1 = 3.74e-2 # h^-1 mM^-1
k_nu_1d = 1.01e-3 # h^-1
k_ni_2
k_nu_2d
k_fb = 2.75e6 # h^-1 mM^-1
k_fbd = 1.84e3 #h^-1

# define a function to solve ODE's
def GeneCircuit(state, t):
    
    # unpack state vector
    g = state[0]
    mRNA = state[1]
    csgA = state[2]
    csgC = state[3]
    csgEG = state[4]
    csgA_per = state[5]
    csgA_per_2 = state[6]
    csgC_per = state[7]
    csgEG_per = state[8]
    csgA__1 = state[9]
    csgA_2 = state[10]
    F = state[11]
    
    
    # compute derivatives
    
    dmRNA = (k_t * g) - (k_td * mRNA)
    # cytoplasm
    dcsgA = (k_A * mRNA) - (k_Ad * csgA) - (k_IM * csgA) - 
            (k_nu_1 * (csgA ** 2.0)) + (k_nu_1d * csgA_2)
    dcsgA_2 = (k_nu_1 * pow(csgA, 2.0)) - (k_nu_1d * csgA_2)
    dcsgC = (k_C * mRNA) - (k_Cd * csgC) - (k_IM * csgC)
    dcsgEG = (k_EG * mRNA) - (k_EGd * csgEG) - (k_IM * csgEG)
    # periplasm
    dcsgA_per = (k_IM * csgA) - (k_Ad * csgA_per) - (k_OM * csgA_per) - 
                (k_nu_1 * pow(csgA_per, 2.)) + (k_nu_1d * csgA_per_2)
    dcsgC_per = (k_IM * csgC) - (k_Cd * csgC_per)
    dcsgEG_per = (k_IM * csgEG) - (k_EGd * csgEG_per)
    # ECM / polymerization
    dcsgA_1 = (kOM * csgA_per) - (k_Ad * csgA_1) - (k_nu_1 * (csgA_1 ** 2.0)) + 
              (k_nu_1d * csgA_2)
    dcsgA_2 
    dF
    
    
    # return derivatives
    return [dmRNA, dcsgA, dcsgA_2, dcsgC, dcsgEG, dcsgA_per, dcsgC_per, 
            dcsgEG_per, dcsgA_1, dcsgA_2, dF]

# define initial conditions
g_0
mRNA_0 = 0.0
csgA_0 = 0.0
csgC_0 = 0.0
csgEG_0 = 0.0
csgA_per_0 = 0.0
csgC_per_0 = 0.0
csgEG_per_0 = 0.0
csgA_1_0 = 0.0
csgA_2_0 = 0.0
F_0 = 0.0

state_0 = [g_0, mRNA_0, csgA_0, csgC_0, csgEG_0, csgA_per_0, csgC_per_0, 
           csgEG_per_0, csgA_1_0, csgA_2_0, F_0]


# define time interval
t = np.arange(0.0, 200.0, 0.1)

# Run simulation using GeneCircuit function
state = odeint(GeneCircuit, state_0, t)

# Plot simulation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = figure()
ax = fig.gca(projection = '3d')
F = state[:, 10]
ax.plot_surface(k_A, k_EG, F)
ax.set_xlabel('rate of csgA translation')
ax.set_ylabel('rate of csgEG translation')
ax.set_zlabel('concentration of fibrils')
