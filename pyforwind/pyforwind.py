"""
An open source package to generate synthetic IEC wind 
fields with extended turbulence characteristics 

"""

#Copyright (C) 2024 GNU LGPLv3, pyforwind 2024, J. Friedrich, D. Moreno, and J. Lengyel.
#All rights reserved.
#Contact: jan.friedrich@uol.de, aura.daniela.moreno.mora@uol.de, janka.lengyel@uol.de
                                                                                                                                                                                                                                                                                                                                             
import numpy as np                                                                                                                                                                                                                                                                                                                                                 
import scipy                                                                                                                                                                   
import scipy.special as sc
import pandas as pd
import logging
import datetime
import time
import timeit

class SWF:
    """
    Generate an IEC-conform Kaimal wind field with extended turbulence characteristics.

    .. legacy:: class

    This class returns a function which generates a ``superstatistical`` Kaimal wind field
    with extended turbulence parametrization.

    Parameters
    ----------

    L : float
        Specifies the integral length scale of longitudinal velocity component.
    mu : float
        Intermittency coefficient that determines deviations from Gaussianity of wind field
        fluctuations.
    V_hub: float
        Wind speed at rotor hub.
    h_hub: float
        Rotor hub height.
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
    full_vector: 
        If True, longitudinal, lateral, and upward Kaimal wind field components will be generated. 
    corr_spec:
        If True, the kinetic energy spectrum will be corrected to exactly reproduce the Gaussian one, i.e.,
        without intermittency corrections.
    log:
        Write parameter list to .info file and logging.
    save: 
        Save wind field.
    """

    def __init__(self, L_int, mu,  V_hub, h_hub, range, dim, kind='spatiotemporal', 
                 H=None, L_c=None, tilde_L=None,  tilde_T=None, sigma=None, N_xi=None, 
                 full_vector=False, corr_spec=False, log=False, save=False): 
        """ Initializes the wind field class."""
        if H is None:
            self.H = 1./3.
        else:
            self.H = H

        if sigma is None:
            sigma = 0.1*V_hub
        
        if full_vector is True:
            self.n_comp = 3
            self.L_int = np.array([L_int, 2.7*L_int/8.1, 0.66*L_int/8.1])
            self.sigma = np.array([sigma, 0.8*sigma, 0.5*sigma])
        else:
            self.n_comp = 1
            self.L_int = np.array([L_int, np.nan, np.nan])
            self.sigma = np.array([sigma, np.nan, np.nan])

        if L_c is None:
            self.L_c = L_int
        else:
            self.L_c = L_c

        self.mu = mu
        self.V_hub = V_hub
        self.h_hub = h_hub
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
    
        self.y = np.linspace(-range[1]/2., range[1]/2., self.N_y)
        self.z = np.linspace(-range[1]/2.+self.h_hub, range[1]/2.+self.h_hub, self.N_z)
        self.N_hub = (self.N_z-1)//2
        if self.y[self.N_hub] != 0.:
            raise ValueError("the grid in the rotor plane must contain the rotor hub.") 
    
        if kind == 'gauss':
            self.field_type = self.gauss_field
        else:
            self.field_type = self.wind_field

        self.kind = kind
        self.full_vector = full_vector
        self.save = save
        self.range = range
        self.dim = dim
        self.corr_spec = corr_spec
        self.log = log
        self.dt = self.T/self.N_x
        self.t = np.linspace(0, self.dt*(self.N_x/2), self.N_x//2+1)
        self.Fs = self.N_x/self.T                                                                                                                                                                                                                                                                                  
        self.f = np.linspace(self.df, 0.5*self.Fs, self.N_x//2+1)
        self.Y, self.Z = np.meshgrid(self.y, self.z)   
        self.positions = np.vstack([self.Y.ravel(), self.Z.ravel()]).T
        self.R_ij = scipy.spatial.distance.cdist(self.positions, self.positions, 'euclidean')

        if log is True:
            self.logger()
        else:
            if save in ['npy', 'bin']:
                raise ValueError("No seed given for log-file.") 
    
    def get_positions(self):
        """ Returns the positions of grid points in the rotor plane. """
        return self.positions
    
    def kaimal_spec(self, F, L, sigma):
        """ Frequency spectrum of Kaimal wind field as specified by IEC 61400-1. """
        return 4.*sigma**2*L/(1.+6.*F*L/self.V_hub)**(1.+2.*self.H)/self.V_hub

    def kaimal_coh(self, F, R):
        """ Coherences of Kaimal wind field as specified by IEC 61400-1. """
        return np.exp(-12.*R*np.sqrt((F/self.V_hub)**2+(0.12/self.L_c)**2))

    def cov(self, L, sigma):
        """ Temporal correlations calculated from inverse Fourier transform of Kaimal spectrum. """
        kaimal_spec_vec = np.vectorize(self.kaimal_spec)
        return np.fft.irfft(kaimal_spec_vec(self.f, L, sigma))[:self.N_x//2]

    def rescaled_spec(self, xi, L, sigma):
        """ Modify/re-scale the Kaimal frequency spectrum 

        Parameters
        ----------
        L : float
            Integral length scale for Kaimal spectrum.
        sigma: float
            Standard deviations of velocity field entering Kaimal spectrum.
        xi : float
            Scale parameter characterizing long- or short-range correlations compared
            to original temporal Kaimal correlations; distributed according to lognormal.

        Returns
        -------
        rescaled_spec : array
            The Fourier transform of the re-scaled temporal Kaimal correlations, i.e., a
            re-scaled spectrum.
        """
        
        with np.errstate(divide='ignore', invalid='ignore'):
            T_rescaled = np.where(self.t==0., 0.,  np.array(xi**np.sqrt(self.mu*np.log(self.tilde_T/(self.t)))
                                                        *(self.t/self.tilde_T)**(self.mu/2)*self.t))
        rescaled_int = np.array(T_rescaled/self.dt, dtype='int')
        if rescaled_int.max() < self.N_x//2-1:
            df = pd.DataFrame({'corr' : self.cov(L, sigma)[rescaled_int]})
            df[df.duplicated()]  = None
            cov_rescaled = df.interpolate(method='linear')
        else:
            cov = self.cov(L, sigma)
            rescaled_tilde = rescaled_int[rescaled_int<self.N_x//2]
            cov_rescaled = np.append(cov[rescaled_tilde], cov[rescaled_tilde][-1]*np.ones(self.N_x//2-rescaled_tilde.size+1))
        
        rescaled_spec = np.fft.rfft(np.append(cov_rescaled, cov_rescaled[::-1][:-1]))
        return rescaled_spec

    def gauss_field(self, L, sigma):
        """ Kaimal wind field components 

        Parameters
        ----------
        L : float
            Integral length scale for Kaimal spectrum (differs for longitudinal, lateral, and vertical
            components).
        sigma: float
            Standard deviations of velocity field entering Kaimal spectrum (differs for 
            longitudinal, lateral, and vertical components).

        Returns
        -------
        u : array
            Returns an original Kaimal wind field of size (N_rotor, N_rotor, N_T). 
        """
       
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        u_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
        for ff in range(1, self.N_x//2+1):
            coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff], self.R_ij))
            u_hat[:,ff] = coh_decomp*np.sqrt(self.kaimal_spec(self.f[ff], L, sigma))@random_phases[:, ff]
        u = np.fft.irfft(u_hat, axis=1).reshape(self.N_y, self.N_y, self.N_x)
        u *= sigma/np.std(u[self.N_hub, self.N_hub, :])
        u -= np.mean(u[self.N_hub, self.N_hub, :])
        return u
    
    def mask_field(self, L, sigma):
        """ Integer mask field for generation of Gaussian scale mixture. """
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        mask_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
        for ff in range(1, self.N_x//2+1):
            coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff], self.R_ij))
            mask_hat[:,ff] = coh_decomp*np.sqrt(self.kaimal_spec(self.f[ff], L, sigma))@random_phases[:, ff]
        mask = np.fft.irfft(mask_hat, axis=1).reshape(self.N_y, self.N_y, self.N_x)
        mask /= np.std(mask)
        mask = 1./2*(1.+sc.erf(mask/np.sqrt(2)))                                                                                           
        mask -= mask.min()
        mask /= mask.max()
        mask *= (self.N_xi-1)
        mask = mask.astype(int)  
        return mask

    def wind_field(self, L, sigma):
        """  Superstatistical Kaimal wind field components 
        
        Parameters
        ----------
        L : float
            Integral length scale for Kaimal spectrum (differs for longitudinal, lateral, and vertical
            components).
        sigma: float
            Standard deviations of velocity field entering Kaimal spectrum (differs for 
            longitudinal, lateral, and vertical components).


        Returns
        -------
        u : array
            Returns a superstatistical Kaimal wind field of size (N_rotor, N_rotor, N_T) of kind
        'temporal', 'spatial', or 'spatiotemporal'.
        """
       
        random_phases = np.exp(1j*np.random.random_sample((self.N_y*self.N_y, self.N_x//2+1))*2*np.pi)
        u_field = np.zeros((self.N_xi, self.N_y, self.N_y, self.N_x))
        u = np.zeros((self.N_y, self.N_y, self.N_x))
        xi_array = np.sort(np.random.lognormal(0, 1, self.N_xi))
        if self.kind in ['spatial']:
            spec = self.kaimal_spec(self.f, L, sigma)
        
        if self.kind in ['temporal']:
            R = self.R_ij
        u_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
        for xx in range(self.N_xi):
            xi = xi_array[xx]
            if self.kind in ['spatial', 'spatiotemporal']:
                with np.errstate(divide='ignore', invalid='ignore'):
                    R = np.where(self.R_ij==0., 0., np.array(xi**np.sqrt(self.mu*np.log(self.tilde_L/(self.R_ij)))
                                                             *(self.R_ij/self.tilde_L)**(self.mu/2)*self.R_ij))
            if self.kind in ['temporal', 'spatiotemporal']:
                spec = self.rescaled_spec(xi, L, sigma)

            for ff in range(1,self.N_x//2+1):
                coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff], R))
                u_hat[:, ff] = coh_decomp*np.sqrt(spec[ff])@random_phases[:, ff]
            u_field[xx] = np.fft.irfft(u_hat, axis=1).reshape(self.N_y, self.N_y, self.N_x)
        mask = self.mask_field(L, sigma)
        for zz in range(self.N_y):
            for yy in range(self.N_y):
                for xx in range(self.N_x):
                    u[zz, yy, xx] = u_field[mask[zz, yy, xx], zz, yy, xx]
        if self.corr_spec is True:
            u = self.corrected_spec(u, random_phases, L, sigma)

        u *= sigma/np.std(u[self.N_hub, self.N_hub, :])
        u -= np.mean(u[self.N_hub, self.N_hub, :])
        return u
    
    def corrected_spec(self, u, random_phases, L, sigma):
        """ Correction of energy spectrum of non-Gaussian wind field.

        Parameters
        ----------
        u : array
            Initial non-Gaussian wind field.
        random_phases: array
            Random-phases of wind field. 
        seed : int
            Initial seed for calculation of corresponding Gauss field.
        Returns
        -------
        u_corrected : array
            Returns a wind field with spectra identical to original Kaimal of size (N_rotor, N_rotor, N_T) 
        """

        u_gauss_hat = np.zeros((self.N_y*self.N_y, self.N_x//2+1), dtype='complex')
        for ff in range(1, self.N_x//2+1):
            coh_decomp = np.linalg.cholesky(self.kaimal_coh(self.f[ff], self.R_ij))
            u_gauss_hat[:,ff] = coh_decomp*np.sqrt(self.kaimal_spec(self.f[ff], L, sigma))@random_phases[:, ff]
        u_gauss = np.fft.irfft(u_gauss_hat, axis=1).reshape(self.N_y, self.N_y, self.N_x)
        u_int_hat = np.fft.fftshift(np.fft.rfftn(u),axes=(0,1))
        u_gauss_hat = np.fft.fftshift(np.fft.rfftn(u_gauss),axes=(0,1))
        u_int_hat *= np.sqrt(np.abs(u_gauss_hat)**2/np.abs(u_int_hat)**2)
        u_corrected = np.fft.irfftn(np.fft.ifftshift(u_int_hat,axes=(0,1)))
        return u_corrected

    def field(self, seed=None):
        """ Final call for wind field.

        Parameters
        ----------
        seed : int, optional
            Random seed for reproducibility of wind fields.

        Returns
        -------
        u : array
            Returns a wind field of size (N_rotor, N_rotor, N_T) of kind
            'gaussian', 'temporal', 'spatial', or 'spatiotemporal'. If full_vector is specified, 
            all three components of the wind field are returned in an array (N_dim=3, N_rotor, N_rotor, N_T) 
            in accordance with the IEC 61400-1.
        """

        if seed is None:
            np.random.seed()
            if self.log is True:
                raise ValueError("No seed given for log-file.") 
            if self.save is True:
                raise ValueError("No seed given for save-file.")
        else:
            np.random.seed(seed)
            if self.log is True:
                logging.info(f"Generating wind field with seed: {seed}")
            u = self.field_type(self.L_int[0], self.sigma[0])+self.V_hub
        if self.full_vector is True:
            v = self.field_type(self.L_int[1], self.sigma[1])
            w = self.field_type(self.L_int[2], self.sigma[2])
            u = np.array([u, v, w])
        if self.save in ['npy', 'bin']:
            self.saving(u, seed)
        return u   
    
    def logger(self):
        """ Initialize a log-file and write wind field model parameters to .info-file."""
        logging.basicConfig(filename=str(self.kind)+'_V_hub_'+str(self.V_hub)+'.log', format='%(asctime)s %(levelname)s %(message)s',
                            filemode='w', level=logging.INFO, force=True)
        logging.info('Instantiating SWF with parameters saved in: '+ str(self.kind)+'_V_hub_'+str(self.V_hub)+'.info') 

        params = np.round(np.array([self.L_int[0], self.L_int[1], self.L_int[2], self.sigma[0], self.sigma[1], 
                                     self.sigma[2], self.mu, self.V_hub, self.h_hub, self.H, self.T, self.range[1], 
                                     self.dim[0], self.dim[1], self.n_comp, self.tilde_L, self.tilde_T, self.N_xi]),3)
        print(params.shape)
        params_str = ['L_int_x', 'L_int_y', 'L_int_z', 'sigma_x', 'sigma_y', 'sigma_z', 'mu', 'V_hub', 'h_hub',
                      'H', 'T', 'diameter', 'N_T', 'N_rotor', 'N_dim', 'tilde_L', 'tilde_T', 'N_xi']
        df = pd.DataFrame([params], columns=params_str)
        df.to_csv(str(self.kind)+'_V_hub_'+str(self.V_hub)+'.info', index=False)  

    def saving(self, u, seed):
        """ Saving of wind field to .npy or .bin file."""
        if self.save in ['npy']:
            np.save('u_'+self.kind+'_V_hub_'+str(self.V_hub)+'_seed_'+str(seed)+'.'+str(self.save), u)
        else:
            np.array(u, dtype='float32').tofile('u_'+self.kind+'_V_hub_'+str(self.V_hub)+'_seed_'+str(seed)+'.'+str(self.save))
        logging.info('Saving wind field to file: u_'+self.kind+'_V_hub_'+str(self.V_hub)+'_seed_'+str(seed)+'.'+str(self.save))
        

   
