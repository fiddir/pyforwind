"""
An open source package to generate synthetic IEC wind 
fields with extended turbulence characteristics 
url="https://github.com/fiddir/pyforwind"
license="LGPL license"
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
import tensorflow as tf

class SWF:
    """
    Generate an IEC-conform Kaimal wind field with extended turbulence characteristics.

    .. legacy:: class

    This class returns a function which generates a ``superstatistical`` Kaimal wind field
    with extended turbulence parametrization.

    Parameters
    ----------

    L : float
        Specifies the integral length scale. 
    mu : float
        Intermittency coefficient that determines deviations from Gaussianity of wind field
        fluctuations.
    V_hub: Wind speed at rotor hub.
    range : (float, float)
        Length scales of wind field in the form (T, diam) where T is the temporal length of wind field
        and diam the extends in the rotor plane.
    dim : (int, int)
        Dimensions of the wind field in the form (N_T, N_rotor), where N_T is the number of grid points in time, and 
        N_rotor^2 the number of grid points in the rotor plane.
    kind : str, optional
        Specifies non-Gaussian features: 'gauss' generates original Kaimal wind field, 'temporal' generates 
        temporal non-Gaussian wind field fluctuations, 'spatial' generates fluctuations in the rotor plane, 
        and 'spatiotemporal' both.
    corrected_spec:
        If True 
    H : float, optional
        Hurst exponent determining the monofractal scaling properties of the wind field, if None it 
        defaults to H=1/3, i.e., the classical ``-5/3``-law of the energy spectrum.
    L_c : float, optional
        The coherence length that characterizes the exponential decay of the coherence in the rotor
        plane.
    tilde_L : float, optional
        The ``internal`` turbulence length scale characterizing non-Gaussian coherence in the rotor 
        plane.
    tilde_T : float, optional
        The ``internal`` turbulence time scale characterizing ``temporal`` non-Gaussianity.
    sigma: float, optional
        Standard deviation at rotor hub, default corresponds to turbulence intensity of 10 per cent.
    N_xi : int, optional
        Number of realizations in Gaussian scale mixture, default is 30 realizations. 
        (convergence check might be appropriate)
    """
    def __init__(self, L, mu,  V_hub, range, dim, kind='spatiotemporal', 
                 H=None, L_c=None, tilde_L=None,  tilde_T=None, sigma=None, N_xi=None): 
        """ Initializes the wind field class."""
        if H is None:
            self.H = 1./3.
        else:
            self.H = H
        
        self.L = L

        if L_c is None:
            self.L_c = L
        else:
            self.L_c = L_c

        self.mu = mu
        self.V_hub = V_hub
        self.T = range[0]
        if tilde_T is None:
            self.tilde_T = self.T
        else:
            self.tilde_T = tilde_T
        
        if tilde_L is None:
            self.tilde_L = V_hub*self.tilde_T
        else:
            self.tilde_L = tilde_L

        self.df = 1./self.T
        self.N_x = dim[0]
        self.N_y = self.N_z = dim[1]
        if N_xi is None:
            self.N_xi = 30
        else:
            self.N_xi = N_xi

        self.y = self.z = np.linspace(-range[1]/2., range[1]/2., self.N_y)
        self.N_hub = (self.N_y-1)//2
        if self.y[self.N_hub] != 0.:
            raise ValueError("the grid in the rotor plane must contain the rotor hub (0,0).") 
    
        if sigma is None:
            self.sigma = 0.1*self.V_hub
        else:
            self.sigma = sigma
        
        if kind == 'gauss':
            self.field = self.gauss_field
        else:
            self.field = self.wind_field

        self.kind = kind
        self.dt = self.T/self.N_x
        self.t = np.linspace(self.dt, self.dt*(self.N_x/2), self.N_x//2+1)
        self.Fs = self.N_x/self.T                                                                                                                                                                                                                                                                                  
        self.f = np.linspace(self.df, 0.5*self.Fs, self.N_x//2+1)
        self.Y, self.Z = np.meshgrid(self.y, self.z)   
        self.positions = np.vstack([self.Y.ravel(), self.Z.ravel()]).T
        self.R_ij = scipy.spatial.distance.cdist(self.positions, self.positions, 'euclidean')
    
    def get_positions(self):
        """ Returns the positions of grid points in the rotor plane. """
        return self.positions
    
    def kaimal_spec(self, F):
        """ Frequency spectrum of Kaimal wind field as specified by IEC 61400-1"""
        return 4.*self.sigma**2*self.L/(1.+6.*F*self.L/self.V_hub)**(1.+2.*self.H)/self.V_hub

    def kaimal_coh(self, F, R):
        """ Coherences of Kaimal wind field as specified by IEC 61400-1"""
        return np.exp(-12.*R*np.sqrt((F/self.V_hub)**2+(0.12/self.L_c)**2))

    def cov(self):
        kaimal_spec_vec = np.vectorize(self.kaimal_spec)
        return np.fft.irfft(kaimal_spec_vec(self.f))[:self.N_x//2]

    def rescaled_spec(self, xi):
        """ Modify/re-scale the Kaimal frequency spectrum 

        Parameters
        ----------
        xi : float
            Scale parameter characterizing long- or short-range correlations compared
            to original temporal Kaimal correlations; distributed according to lognormal.

        Returns
        -------
        rescaled_spec : array
            The Fourier transform of the re-scaled temporal Kaimal correlations, i.e., a
            re-scaled spectrum.
        """
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
        """ 

        Parameters
        ----------
        seed : int
            Random seed for reproducibility of wind fields.

        Returns
        -------
        u : array
            Returns an original Kaimal wind field of size (N_rotor, N_rotor, N_T). 
        """
        np.random.seed(seed)
        print(seed)
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        u_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
        for ff in range(1, self.N_x//2+1):
            coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff], self.R_ij))
            u_hat[:,ff] = coh_decomp*np.sqrt(self.kaimal_spec(self.f[ff]))@random_phases[:, ff]
        u = np.fft.irfft(u_hat, axis=1).reshape(self.N_y, self.N_y, self.N_x)
        u *= self.sigma/np.std(u[self.N_hub, self.N_hub, :])
        u += self.V_hub - np.mean(u[self.N_hub, self.N_hub, :])
        return u
    
    def mask_field(self):
        """ Integer mask field for generation of Gaussian scale mixture. """
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        mask_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
        for ff in range(1, self.N_x//2+1):
            coh_decomp = scipy.linalg.cholesky(self.kaimal_coh(self.f[ff], self.R_ij), overwrite_a=True, check_finite=False)
            mask_hat[:,ff] = coh_decomp*np.sqrt(self.kaimal_spec(self.f[ff]))@random_phases[:, ff]
        mask = np.fft.irfft(mask_hat, axis=1).reshape(self.N_y, self.N_y, self.N_x)
        mask /= np.std(mask)
        mask = 1./2*(1.+sc.erf(mask/np.sqrt(2)))                                                                                           
        mask -= mask.min()
        mask /= mask.max()
        mask *= (self.N_xi-1)
        mask = mask.astype(int)  
        return mask

    def wind_field(self, seed):
        """ 

        Parameters
        ----------
        seed : int
            Random seed for reproducibility of wind fields.

        Returns
        -------
        u : array
            Returns a superstatistical Kaimal wind field of size (N_rotor, N_rotor, N_T) of kind
            'temporal', 'spatial', or 'spatiotemporal'.
        """
        np.random.seed(seed)
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        u_field = np.zeros((self.N_xi, self.N_y, self.N_y, self.N_x))
        u = np.zeros((self.N_y, self.N_y, self.N_x))
        xi_array = np.sort(np.random.lognormal(0, 1, self.N_xi))
        if self.kind in ['spatial']:
            spec = self.kaimal_spec(self.f)
        
        if self.kind in ['temporal']:
            R = self.R_ij
        
        for xx in range(self.N_xi):
            xi = xi_array[xx]
            u_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
            if self.kind in ['spatial', 'spatiotemporal']:
                with np.errstate(divide='ignore'):
                    R = np.where(self.R_ij==0., 0., np.array(xi**np.sqrt(self.mu*np.log(self.tilde_L/(self.R_ij)))*(self.R_ij/self.tilde_L)**(self.mu/2)*self.R_ij))

            if self.kind in ['temporal', 'spatiotemporal']:
                spec = self.rescaled_spec(xi)
            
            for ff in range(1,self.N_x//2+1):
                coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff], R))
                u_hat[:, ff] = coh_decomp*np.sqrt(spec[ff])@random_phases[:, ff]
            u_field[xx] = np.fft.irfft(u_hat, axis=1).reshape(self.N_y, self.N_y, self.N_x)
        mask = self.mask_field()
        for zz in range(self.N_y):
            for yy in range(self.N_y):
                for xx in range(self.N_x):
                    u[zz, yy, xx] = u_field[mask[zz, yy, xx], zz, yy, xx]
        u *= self.sigma/np.std(u[self.N_hub, self.N_hub, :])
        u += self.V_hub - np.mean(u[self.N_hub, self.N_hub, :])
        return u
    
   