#!/usr/bin/env python3

"""
An open source package to generate synthetic IEC wind 
fields with extended turbulence characteristics 
url="https://github.com/fiddir/pyforwind"
license="MIT license"
"""

import matplotlib                                                                                                                                                                                                                                                                                                                                               
import numpy as np
import math                                                                                                                                                                                                                                                                                                                                                   
import scipy                                                                                                                                                                   
import scipy.special as sc
from matplotlib import pyplot as plt
from scipy.misc import derivative
import pandas as pd
import logging

class SWF:
    def __init__(self, H, L, L_c, mu, dy, dz, y_i, y_e, z_i, z_e, T, tilde_T, N_x, N_y, N_z, N_xi, V_hub, sigma=None): 
        self.H = H
        self.L = L
        self.L_c = L_c
        self.mu = mu
        self.dy = dy
        self.dz = dz
        self.y_i = y_i
        self.y_e = y_e
        self.z_i = z_i
        self.z_e = z_e
        self.T = T
        self.tilde_T = tilde_T
        self.df = 1./self.T
        self.N_x = N_x
        self.N_y = N_y
        self.N_z = N_z
        self.N_xi = N_xi
        self.V_hub = V_hub
        if sigma is None:
             self.sigma = 1. 
        else:
            self.sigma = sigma
        self.dt = self.T/self.N_x
        self.t = np.linspace(self.dt, self.dt*(self.N_x/2), self.N_x//2+1)
        self.y = np.linspace(y_i, y_e, self.N_y)
        self.z = np.linspace(z_i, z_e, self.N_z)      
        self.Fs = self.N_x/self.T                                                                                                                                                                                                                                                                                  
        self.f = np.linspace(self.df, 0.5*self.Fs, self.N_x//2+1)
        self.Y, self.Z = np.meshgrid(self.y, self.z)   
        positions = np.vstack([self.Y.ravel(), self.Z.ravel()]).T
        self.R_ij = scipy.spatial.distance.cdist(positions, positions, 'euclidean')
    
    def kaimal_spec(self, F):
        return 4.*self.sigma**2*self.L/(1.+6.*F*self.L/self.V_hub)**(1.+2.*self.H)/self.V_hub
    
    def kaimal_coh(self, F):
        return np.exp(-12.*self.R_ij*np.sqrt((F/self.V_hub)**2+(0.12/self.L_c)**2))
    
    def cov(self):
        kaimal_spec_vec = np.vectorize(self.kaimal_spec)
        return np.fft.irfft(kaimal_spec_vec(self.f))[:self.N_x//2]

    def rescaled_spec(self, xi):
        T_rescaled = np.array(xi**np.sqrt(self.mu*np.log(self.tilde_T/(self.t)))*(self.t/self.tilde_T)**(self.mu/2)*self.t)
        rescaled_int = np.array(T_rescaled/self.dt, dtype='int')    
        if rescaled_int.max() < self.N_x//2-1:
            df = pd.DataFrame({'corr' : self.cov()[rescaled_int]})  
            df[df.duplicated()]  = None
            cov_rescaled = df.interpolate(method='linear')
        else:
            cov = self.cov()
            rescaled_tilde = rescaled_int[rescaled_int<self.N_x//2]
            cov_rescaled = np.append(cov[rescaled_tilde], cov[rescaled_tilde][-1]*np.ones(self.N_x//2-rescaled_tilde.size+1))
        rescaled_spec = np.fft.rfft(np.append(cov_rescaled, cov_rescaled[::-1][1:]))  
        return rescaled_spec  
    
    def gauss_field(self, seed):
        np.random.seed(seed)
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        u_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
        for ff in range(1, self.N_x//2+1):
            coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff]))
            u_hat[:,ff] = coh_decomp*np.sqrt(self.kaimal_spec(self.f[ff]))@random_phases[:, ff]
        u = np.fft.irfft(u_hat, axis=1)
        return u.reshape(self.N_y, self.N_y, self.N_x)
    
    def mask_field(self):
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        u_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
        for ff in range(1, self.N_x//2+1):
            coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff]))
            u_hat[:,ff] = coh_decomp*np.sqrt(self.kaimal_spec(self.f[ff]))@random_phases[:, ff]
        u = np.fft.irfft(u_hat, axis=1)
        return u.reshape(self.N_y, self.N_y, self.N_x)

    def field(self, seed):
        np.random.seed(seed)
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        u_field = np.zeros((self.N_xi, self.N_y, self.N_y, self.N_x))
        u = np.zeros((self.N_y, self.N_y, self.N_x))
        xi_array = np.sort(np.random.lognormal(0, 1, self.N_xi))
        for xx in range(self.N_xi):
            xi = xi_array[xx]
            print('xi_values', xi)
            u_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
            spec = self.rescaled_spec(xi)
            for ff in range(1,self.N_x//2+1):
                coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff]))
                u_hat[:, ff] = coh_decomp*np.sqrt(spec[ff])@random_phases[:, ff]
            u_field[xx] = np.fft.irfft(u_hat, axis=1).reshape(self.N_y, self.N_y, self.N_x)
        mask = self.mask_field()
        mask /= np.std(mask)
        mask = 1./2*(1.+sc.erf(mask/np.sqrt(2)))                                                                                           
        mask -= mask.min()
        mask /= mask.max()
        mask *= (self.N_xi-1)
        mask = mask.astype(int)  
        for xx in range(self.N_y):
            for yy in range(self.N_y):
                for zz in range(self.N_x):
                    u[xx, yy, zz] = u_field[mask[xx, yy, zz], xx, yy, zz]
        f = open('field_seed_'+str(seed)+'.log', "x")
        f.write(str(list(self.__dict__)[:14]))
        f.write(str(list(self.__dict__.values())[:14]))
        f.close()
        return u, mask
    
