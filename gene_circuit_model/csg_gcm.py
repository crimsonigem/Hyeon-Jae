#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 08:47:59 2017

@author: Hyeon-Jae
"""

from scipy.integrate import odeint
import numpy as np
import pandas as pd
import math
import time

### TIME INTERVAL
t = np.arange(0.0, 60*60, 1) 

### CONSTANTS
# Transcription
alpha1 = 0.0921 # csgA (seq/sec)
alpha2 = 0.0214 # csgBCEFG (seq/sec)

# mRNA Degradation
half_life = 408 # seconds
zeta1 = math.log(2)/half_life
zeta2 = math.log(2)/half_life

# Translation
beta0 = 0.1 # (protein/sec)
brange = np.arange(0.0, 2.1, 0.1)
beta = beta0 * brange

# Periplasmic Export
gamma1 = 0.25 # uM^(-4)sec^(-1)
gamma1d = 0.025 # sec^(-1)
gamma2 = 0.0085 # uM^(-4)sec^(-1)
gamma3 = 1.00 # sec^(-1)

# Extracellular Secretion
delta1 = 0.76 
delta1d = 28e-3 # s^-1
delta3 = 0.76 
delta3d = 28e-3 # s^-1
delta4 = 0.0384 # uM^(-1)sec^(-1)
delta4d = 0.76 # also try 0.0385 ????? (homotetramer formation or dimerization?)
delta5 = 0.25 # uM^(-4)sec^(-1)

# Diffusion
N_A = 6.022e23
omega = 10e-9
D = 1e-10
SA = 4.42e-12

# Aggregation / Polymerization 
epsilon1 = 1.038e-8 # s^(-1)uM^(-1)
epsilon1d = 2.805e-9 # s^(-1)uM^(-1)
epsilon2 = 0.764 # s^(-1)uM^(-1)
epsilon2d = 5.111e-4 # s^(-1)uM^(-1)

### ODE SOLVER
def CsgPathway(state, t, beta1, beta2, beta3, beta4, beta5, beta6):
    
    # unpack state vector
    g_csgA = state[0]
    g_csgBCEFG = state[1]
    mRNA_csgA = state[2]
    mRNA_csgBCEFG = state[3]
    csgA_cyt = state[4]
    secBcsgA = state[5]
    secABYEGcsgA = state[6]
    F_cyt = state[7]
    csgB_cyt = state[8]
    secBcsgB = state[9]
    secABYEGcsgB = state[10]
    csgC_cyt = state[11]
    secBcsgC = state[12]
    secABYEGcsgC = state[13]
    csgE_cyt = state[14]
    secBcsgE = state[15]
    secABYEGcsgE = state[16]
    csgF_cyt = state[17]
    secBcsgF = state[18]
    secABYEGcsgF = state[19]
    csgG_cyt = state[20]
    secBcsgG = state[21]
    secABYEGcsgG = state[22]
    secB = state[23]
    secA = state[24]
    secYEG = state[25]
    csgE_9 = state[26]
    csgG_9 = state[27]
    csgGEF = state[28]
    csgA_per = state[29]
    F_per = state[30]
    csgB_per = state[31]
    csgC_per = state[32]
    csgCcsgA = state[33]
    csgCcsgB = state[34]
    csgE_per = state[35]
    csgF_per = state[36]
    csgG_per = state[37]
    csgF_ECM = state[38]
    csgA_ECM = state[39]
    F_ECM = state[40]
    csgB_ECM = state[41]
    
    # compute derivatives
    dmRNA_csgA = (alpha1 * g_csgA) - (zeta1 * mRNA_csgA)
    dmRNA_csgBCEFG = (alpha2 * g_csgBCEFG) - (zeta2 * mRNA_csgBCEFG)
    
    dcsgA_cyt = (beta1 * mRNA_csgA) - (gamma1 * csgA_cyt * (secB ** 4)) + (gamma1d * secBcsgA) - (epsilon1 * (csgA_cyt ** 2)) + (epsilon1d * F_cyt) - (epsilon2 * F_cyt * csgA_cyt) + (epsilon2d * F_cyt)
    dsecBcsgA = (gamma1 * csgA_cyt * (secB **4)) - (gamma1d * secBcsgA) - (gamma2 * secBcsgA * (secA ** 2) * (secYEG ** 2))
    dsecABYEGcsgA = (gamma2 * secBcsgA * (secA ** 2) * (secYEG ** 2)) - (gamma3 * secABYEGcsgA)
    dF_cyt = (epsilon1 * (csgA_cyt ** 2)) - (epsilon1d * F_cyt)
    dcsgB_cyt = (beta2 * mRNA_csgBCEFG) - (gamma1 * csgB_cyt * (secB ** 4)) + (gamma1d * secBcsgB)
    dsecBcsgB = (gamma1 * csgB_cyt * (secB **4)) - (gamma1d * secBcsgB) - (gamma2 * secBcsgB * (secA ** 2) * (secYEG ** 2))
    dsecABYEGcsgB = (gamma2 * secBcsgB * (secA ** 2) * (secYEG ** 2)) - (gamma3 * secABYEGcsgB)
    dcsgC_cyt = (beta3 * mRNA_csgBCEFG) - (gamma1 * csgC_cyt * (secB ** 4)) + (gamma1d * secBcsgC)
    dsecBcsgC = (gamma1 * csgC_cyt * (secB **4)) - (gamma1d * secBcsgB) - (gamma2 * secBcsgB * (secA ** 2) * (secYEG ** 2))
    dsecABYEGcsgC = (gamma2 * secBcsgC * (secA ** 2) * (secYEG ** 2)) - (gamma3 * secABYEGcsgC)
    dcsgE_cyt = (beta4 * mRNA_csgBCEFG) - (gamma1 * csgE_cyt * (secB ** 4)) + (gamma1d * secBcsgE)
    dsecBcsgE = (gamma1 * csgE_cyt * (secB **4)) - (gamma1d * secBcsgE) - (gamma2 * secBcsgE * (secA ** 2) * (secYEG ** 2))
    dsecABYEGcsgE = (gamma2 * secBcsgE * (secA ** 2) * (secYEG ** 2)) - (gamma3 * secABYEGcsgE)
    dcsgF_cyt = (beta5 * mRNA_csgBCEFG) - (gamma1 * csgF_cyt * (secB ** 4)) + (gamma1d * secBcsgF)
    dsecBcsgF = (gamma1 * csgF_cyt * (secB **4)) - (gamma1d * secBcsgF) - (gamma2 * secBcsgF * (secA ** 2) * (secYEG ** 2))
    dsecABYEGcsgF = (gamma2 * secBcsgF * (secA ** 2) * (secYEG ** 2)) - (gamma3 * secABYEGcsgF)
    dcsgG_cyt = (beta6 * mRNA_csgBCEFG) - (gamma1 * csgG_cyt * (secB ** 4)) + (gamma1d * secBcsgB)
    dsecBcsgG = (gamma1 * csgG_cyt * (secB **4)) - (gamma1d * secBcsgG) - (gamma2 * secBcsgG * (secA ** 2) * (secYEG ** 2))
    dsecABYEGcsgG = (gamma2 * secBcsgG * (secA ** 2) * (secYEG ** 2)) - (gamma3 * secABYEGcsgG)
    dsecB = (gamma3 * (secABYEGcsgA + secABYEGcsgB + secABYEGcsgC + secABYEGcsgE + secABYEGcsgF + secABYEGcsgG)) - (gamma1 * (secB ** 4) * (csgA_cyt + csgB_cyt + csgC_cyt + csgE_cyt + csgF_cyt + csgG_cyt)) + (gamma1d * (secBcsgA + secBcsgB + secBcsgC + secBcsgE + secBcsgF + secBcsgG))
    dsecA = (gamma3 * (secABYEGcsgA + secABYEGcsgB + secABYEGcsgC + secABYEGcsgE + secABYEGcsgF + secABYEGcsgG)) - (gamma2 * (secA ** 2) * (secYEG ** 2) * (secBcsgA + secBcsgB + secBcsgC + secBcsgE + secBcsgF + secBcsgG))
    dsecYEG = (gamma3 * (secABYEGcsgA + secABYEGcsgB + secABYEGcsgC + secABYEGcsgE + secABYEGcsgF + secABYEGcsgG)) - (gamma2 * (secA ** 2) * (secYEG ** 2) * (secBcsgA + secBcsgB + secBcsgC + secBcsgE + secBcsgF + secBcsgG))
    
    dcsgE_9 = (delta1 * (csgE_per ** 9)) - (delta1d * csgE_9)
    dcsgG_9 = (delta3 * (csgG_per ** 9)) - (delta3d * csgG_9) - (delta4 * csgE_9 * (csgF_ECM ** 2) * csgG_9) + (delta4d * csgGEF)
    dcsgGEF = (delta4 * csgE_9 * (csgF_ECM ** 2) * csgG_9) - (delta4d * csgGEF)
    dcsgA_per = (gamma3 * secABYEGcsgA) - (delta5 * csgA_per * csgC_per) - (epsilon1 * (csgA_per ** 2)) + (epsilon1d * F_per) - (epsilon2 * F_per * csgA_per) + (epsilon2d * F_per)
    dF_per = (epsilon1 * (csgA_per ** 2)) - (epsilon1d * F_per)
    dcsgB_per = (gamma3 * secABYEGcsgB) - (delta5 * csgB_per * csgC_per)
    dcsgC_per = (gamma3 * secABYEGcsgC) - (delta5 * csgA_per * csgC_per) - (delta5 * csgB_per * csgC_per) + (D * ((csgCcsgA - csgA_ECM) * 1e-3 / omega) * SA * N_A) + (D * ((csgCcsgB - csgB_ECM) * 1e-3 / omega) * SA * N_A)
    dcsgCcsgA = (delta5 * csgA_per * csgC_per) - (D * ((csgCcsgA - csgA_ECM) * 1e-3 / omega))
    dcsgCcsgB = (delta5 * csgB_per * csgC_per) - (D * ((csgCcsgB - csgB_ECM) * 1e-3 / omega))
    dcsgE_per = (gamma3 * secABYEGcsgE)
    dcsgF_per = (gamma3 * secABYEGcsgF) - D * ((csgF_per - csgF_ECM) * 1e-3 / omega) * SA * N_A
    dcsgG_per = (gamma3 * secABYEGcsgF) - (delta3 * (csgG_per ** 9))
    
    dcsgF_ECM = D * ((csgF_per - csgF_ECM) * 1e-3 / omega) * SA * N_A
    dcsgA_ECM = (D * ((csgCcsgA - csgA_ECM) * 1e-3 / omega) * SA * N_A) - (epsilon1 * (csgA_cyt ** 2)) + (epsilon1d * F_cyt) - (epsilon2 * F_cyt * csgA_cyt) + (epsilon2d * F_cyt)
    dF_ECM = (epsilon1 * (csgA_ECM ** 2)) - (epsilon1d * F_ECM)
    dcsgB_ECM = D * ((csgCcsgB - csgB_ECM) * 1e-3 / omega) * SA * N_A
    
    # return derivatives
    return [g_csgA, g_csgBCEFG, dmRNA_csgA, dmRNA_csgBCEFG, dcsgA_cyt, dsecBcsgA, dsecABYEGcsgA, dF_cyt, dcsgB_cyt,
            dsecBcsgB, dsecABYEGcsgB, dcsgC_cyt, dsecBcsgC, dsecABYEGcsgC, dcsgE_cyt,
            dsecBcsgE, dsecABYEGcsgE, dcsgF_cyt, dsecBcsgF, dsecABYEGcsgF, dcsgG_cyt, dsecBcsgG,
            dsecABYEGcsgG, dsecB, dsecA, dsecYEG, dcsgE_9, dcsgG_9, dcsgGEF, dcsgA_per, dF_per,
            dcsgB_per, dcsgC_per, dcsgCcsgA, dcsgCcsgB, dcsgE_per, dcsgF_per, dcsgG_per, 
            dcsgF_ECM, dcsgA_ECM, dF_ECM, dcsgB_ECM]
    
### INITIAL CONDITIONS
g_csgA_0 = 71.4e-8 # uM
g_csgBCEFG_0 = 71.4e-8 # uM
mRNA_csgA_0 = 0.0
mRNA_csgBCEFG_0 = 0.0

csgA_cyt_0 = 0.0
secBcsgA_0 = 0.0
secABYEGcsgA_0 = 0.0
F_cyt_0 = 0.0
csgB_cyt_0 = 0.0
secBcsgB_0 = 0.0
secABYEGcsgB_0 = 0.0
csgC_cyt_0 = 0.0
secBcsgC_0 = 0.0
secABYEGcsgC_0 = 0.0
csgE_cyt_0 = 0.0
secBcsgE_0 = 0.0
secABYEGcsgE_0 = 0.0
csgF_cyt_0 = 0.0
secBcsgF_0 = 0.0
secABYEGcsgF_0 = 0.0
csgG_cyt_0 = 0.0
secBcsgG_0 = 0.0
secABYEGcsgG_0 = 0.0
SecB_0 = 0.4 # uM
SecA_0 = 0.2
SecYEG_0 = 0.2

csgG_9_0 = 0.0
csgE_9_0 = 0.0
csgGEF_0 = 0.0
csgA_per_0 = 0.0
F_per_0 = 0.0
csgB_per_0 = 0.0
csgC_per_0 = 0.0
csgCcsgA_0 = 0.0
csgCcsgB_0 = 0.0
csgE_per_0 = 0.0
csgF_per_0 = 0.0
csgG_per_0 = 0.0

csgF_ECM_0 = 0.0
csgA_ECM_0 = 0.0
F_ECM_0 = 0.0
csgB_ECM_0 = 0.0
    
state_0 = [g_csgA_0, g_csgBCEFG_0, mRNA_csgA_0, mRNA_csgBCEFG_0, csgA_cyt_0, secBcsgA_0, 
           secABYEGcsgA_0, F_cyt_0, csgB_cyt_0, secBcsgB_0, secABYEGcsgB_0, csgC_cyt_0, 
           secBcsgC_0, secABYEGcsgC_0, csgE_cyt_0, secBcsgE_0, secABYEGcsgE_0, csgF_cyt_0, 
           secBcsgF_0, secABYEGcsgF_0, csgG_cyt_0, secBcsgG_0, secABYEGcsgG_0, SecB_0, SecA_0, 
           SecYEG_0, csgG_9_0, csgE_9_0, csgGEF_0, csgA_per_0, F_per_0, csgB_per_0, csgC_per_0,
           csgCcsgA_0, csgCcsgB_0, csgE_per_0, csgF_per_0, csgG_per_0, csgF_ECM_0, csgA_ECM_0,
           F_ECM_0, csgB_ECM_0]

### RUN SIMULATION
state = odeint(CsgPathway, state_0, t, args=(beta[10], 1, 1, 1, 1, 1), mxstep=5000000)

# Maximize A
state = np.zeros((3600, 21))
state_df = pd.DataFrame(state, index=t, columns=np.arange(0, 21, 1))
for i in range(21):
    tmp = odeint(CsgPathway, state_0, t, args=(beta[i], 1, 1, 1, 1, 1))
    tmp_df = pd.DataFrame(tmp, index=t, columns=np.arange(0, 42, 1))
    state_df.iloc[:, i] = tmp_df.iloc[:, 40]

### CHECK NUMBERS
test_range = np.arange(0, 60, 10)
state_df = pd.DataFrame(state, index=t, columns=np.arange(0, 42, 1))
print(state_df.shape)
#for i in test_range:
 #   print(state_df.iloc[i, 40])

### VISUALIZATION
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')

# Plot the surface.
surf = ax.plot_surface(t, beta, state_df, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
