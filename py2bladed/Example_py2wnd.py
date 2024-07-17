#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:13:57 2024

@author: danielamoreno
"""
#from pyforwind import SWF

import sys
import scipy.io as si

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc
from scipy.special import gamma, factorial
import seaborn as sns
from matplotlib import rc
import scipy.stats as stats
import scipy  
import math
import scipy.special as sc
import os
from scipy.special import gamma, factorial
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
plt.rc('font', family='serif')

run 'pyforwind.py'

# --- set the wind field box ---
N_T = 1024//2 # grid points in "temporal" direction
N_rotor = 21 # rotor plane
N_xi = 30#30 # number of realization in Gaussian scale mixture 


# --- Model Parameters ---
H = 1./3. # Hurst exponent, determines power law of Kaimal spectrum                                                                                                                                                 
mu = 0.22 # intermittency coefficient                                                                                                                                                                               
L = 340.2 # integral length scale in Kaimal spectrum                                                                                                                                                                
L_c = 340.2 # correlation length in Kaimal coherences                                                                                                                                                               
T = 600. # length of time series 
diam = 200. # diameter
tilde_T = 600.
V_hub = 8.
tilde_L = tilde_T*V_hub

# --- Saving --
path = '/WindFields'

# --- Generate Fields --- 
realizations = 1 # number of realizations
n_comp = 3 # number of components of the field (u, v, w)

swf_gauss = SWF(L, mu, V_hub, (T, diam), (N_T, N_rotor), kind='gauss', full_vector=True)
u_gauss = np.zeros((realizations, n_comp, N_rotor, N_rotor, N_T))
#u_temporal = np.zeros((realizations,n_comp, N_rotor, N_rotor, N_T))
#u_spatiotemporal = np.zeros((realizations,n_comp, N_rotor, N_rotor, N_T))

for nn in range(realizations):
    seed = 31 
    u_gauss[nn] = swf_gauss.field(seed)
    #u_temporal[nn] = swf_temporal.field(seed)
    #u_spatiotemporal[nn] = swf_spatiotemporal.field(seed)
    file_name_base_gauss = os.path.join(path,'Gauss_' + str(seed))
    #np.save(file_name_base_gauss+'_u.npy', u_gauss)
    # Generating log file
    logfilename= file_name_base_gauss+'.log'
    ts = time.time()
    str_time = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H:%M:%S')            
    f = open(logfilename, "a")
    f.write('time stamp: ' + str_time + '\n' )
    f.write('arguments of wind field:' + '\n')
    f.write(str(list(swf_gauss.__dict__)[:17]) + '\n')
    f.write(str(list(swf_gauss.__dict__.values())[:17]) + '\n'+ '\n')
    f.close()
    # Generating mat file
    si.savemat(file_name_base_gauss+'.mat',{'u_gauss':u_gauss,'u_param_names':list(swf_gauss.__dict__)[:17], 'u_param':list(swf_gauss.__dict__.values())[:17] } )
    seed += 1

    

    
    
