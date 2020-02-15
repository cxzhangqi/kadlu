# ================================================================================ #
#   Authors: Casey Hillard and Oliver Kirsebom                                     #
#   Contact: oliver.kirsebom@dal.ca                                                #
#   Organization: MERIDIAN (https://meridian.cs.dal.ca/)                           #
#   Team: Data Analytics                                                           #
#   Project: kadlu                                                                 #
#   Project goal: The kadlu library provides functionalities for modeling          #
#   underwater noise due to environmental source such as waves.                    #
#                                                                                  #
#   License: GNU GPLv3                                                             #
#                                                                                  #
#       This program is free software: you can redistribute it and/or modify       #
#       it under the terms of the GNU General Public License as published by       #
#       the Free Software Foundation, either version 3 of the License, or          #
#       (at your option) any later version.                                        #
#                                                                                  #
#       This program is distributed in the hope that it will be useful,            #
#       but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#       GNU General Public License for more details.                               # 
#                                                                                  #
#       You should have received a copy of the GNU General Public License          #
#       along with this program.  If not, see <https://www.gnu.org/licenses/>.     #
# ================================================================================ #

""" Parabolic Equation propagator module within the kadlu library

    This module provides an implementation of the Parabolic Equation 
    propagation scheme for numerically solving the wave equation 
    for sound pressure.

    Contents:
        Propagator class
"""
import numpy as np
import math
from tqdm import tqdm
from numpy.lib import scimath
from kadlu.utils import create_boolean_array
from kadlu.geospatial.data_sources.source_map import default_val


class Propagator():
    """ Propagates the sound pressure field from zero range to the boundary 
        of the computational domain using a parabolic-equation numerical scheme.
        
        Args:
            ocean: Ocean
                Ocean data
            sound_speed: SoundSpeed
                Sound speed data
            grid: Grid
                Computational grid
            k0: float
                Reference wavenumber in inverse meters
            bathy_step: int
                How often the bathymetry data is updated. If for example bathy_step=3, the 
                bathymetry is updated at every 3rd step.
            steps_btw_sound_speed_updates: int
                How often the sound-speed data is updated. If for example steps_btw_sound_speed_updates=3, 
                the sound speed is updated at every 3rd step. By default steps_btw_sound_speed_updates is 
                set to infinity, corresponding to a range-independent sound speed profile.
            verbose: bool
                Print information during execution
            progress_bar: bool
                Show progress bar. Only shown if verbose if False.            

        Attributes:

        Example:
    """
    def __init__(self, ocean, seafloor, c, c0, grid, k0,\
                smooth_len_den, smooth_len_c,\
                absorption_layer, ignore_bathy_gradient,\
                bathy_step, c_step,\
                verbose=False, progress_bar=True):

        self.ocean = ocean
        self.seafloor = seafloor
        self.c = c
        self.c0 = c0
        self.grid = grid
        self.k0 = k0

        self.verbose = verbose
        if self.verbose:
            self.progress_bar = False
        else:
            self.progress_bar = progress_bar

        self.smooth_len_den = smooth_len_den
        self.smooth_len_c = smooth_len_c

        self.ignore_bathy_gradient = ignore_bathy_gradient

        self.update_bathy = create_boolean_array(self.grid.Nr, bathy_step)
        self.update_c = create_boolean_array(self.grid.Nr, c_step)

        # compute cos(theta) and sin(theta) for all angular bins
        self.costheta = np.cos(grid.q)
        self.sintheta = np.sin(grid.q)

        # refractive index squared of seafloor
        self.n2_b = seafloor.nsq() 
        self.n2_w = np.zeros(self.grid.Z.shape)

        # attenuation at lower and upper surface of computational domain
        absorp_coeff =  1. / np.log10(np.e) / np.pi
        absorp_thick = absorption_layer * np.max(np.abs(grid.z))
        D = absorp_thick / 3
        self.attenuation = 1j * absorp_coeff * np.exp(-(np.abs(grid.Z) - np.max(np.abs(grid.z)))**2 / D**2)


    def run(self, psi, receiver_depth=[.1], vertical_slice=True):
        """ Propagate the pressure field to the boundary of the computional domain.

            The sound pressure is computed at every grid point on 
            one or several horizontal planes at user-specified depth(s).

            Optionally, the sound pressure may also be computed at every 
            every grid point on a vertical plane intersecting the source position.
            
            Args:
                psi: 2d numpy array
                    Starting sound pressure field computed with the PEStarter.
                    Has shape (Nz,Nq) where Nz and Nq are the number of 
                    vertical and angular grid points, respectively.
                receiver_depth: list of floats
                    Depths at which horizontal slices of the sound pressure 
                    field are computed.
                vertical_slice: bool
                    Compute the sound pressure at all grid points on 
                    a vertical plane intersecting the source position. 
                    Note: This will slow down the computation.

            Returns:
                output: OutputCollector
                    Result of the computation
        """
        # initialize output containers
        self._init_output(psi.shape[0], receiver_depth, vertical_slice) 

        # free propagator
        U_fr = self._free_propagation_matrix()

        # save output field at 0
        self._save_output(step_no=0, dist=0, psi=psi)

        # PE marching starts here
        dist = 0
        Nr = self.grid.Nr
        dr = self.grid.dr
        for i in tqdm(range(Nr-1), disable = not self.progress_bar):

            # (1) r --> r + dr/2 free propagation
            psi = U_fr * psi

            # update bathy
            if self.update_bathy[i]:
                self._update_bathy(dist + dr/2)

            # update sound speed
            if self.update_c[i]:
                self._update_c(dist + dr/2)

            # compute propagation matrix
            if self.update_bathy[i] or self.update_c[i]:
                U = self._propagation_matrix()            

            # (2) phase adjustment at r + dr/2            
            psi = np.fft.fft(U * np.fft.ifft(psi, axis=1), axis=1)    

            # (3) r + dr/2 --> r + dr free propagation
            psi = U_fr * psi

            # increment distance
            dist = dist + dr

            # collect output
            self._save_output(step_no=i+1, dist=dist, psi=psi)


    def _free_propagation_matrix(self):
        """ Compure matrices used to propagate the sound pressure field 
            in the radial direction by half a grid spacing.

            The matrix has shape (Nz,Nq) where Nz and Nq are the number 
            of vertical and angular grid points, respectively.

            Returns:
                U_fr: 2d numpy array
                    Half-step free propagation matrix. 
        """
        k0 = self.k0
        dr = self.grid.dr
        kz = self.grid.kz
        Nq = self.grid.Nq

        U_fr = np.exp(1j * dr / 2 * (scimath.sqrt(k0**2 - kz**2) - k0))
        U_fr = U_fr[:,np.newaxis]
        U_fr = U_fr * np.ones(shape=(1,Nq))
        U_fr = U_fr[np.newaxis,:,:]

        return U_fr


    def _propagation_matrix(self):

        h = self.height_above_seafloor
        lc = self.smooth_len_c
        ld = self.smooth_len_den
        #den_w = self.ocean.water_density
        den_w = default_val['water_density']
        #den_b = self.seafloor.density
        # in the future, fetch/load sediment values here ?
        den_b = default_val['seafloor_density']

        # smooth sound speed
        f = 0.5 * (1 + np.tanh(0.5 * h / lc))
        n2 = self.n2_w + (self.n2_b - self.n2_w) * f
        itmp = (self.depth[0,:] == 0) 
        if np.any(itmp):
            n2[0,itmp] = self.n2_b

        # smooth density
        th = np.tanh(0.5 * h / ld)
        f = 0.5 * (1 + th)
        den = den_w + (den_b - den_w) * f
        sech2 = 1. / np.cosh(0.5 * h / ld)
        sech2 = sech2**2
        dden =  0.5 * sech2 / ld * np.sqrt(1 + self.gradient**2)
        dden =  0.5 * (den_b - den_w) * dden
        d2den = -0.5 * sech2 / ld * (th / ld * (1 + self.gradient**2))
        d2den = 0.5 * (den_b - den_w) * d2den

        # calculate propagation matrix, U
        U = np.exp(1j * self.grid.dr * self.k0 * (-1 + scimath.sqrt(n2 + self.attenuation +\
            0.5 / self.k0**2 * (d2den / den - 1.5 * (dden / den)**2))))

        U = U[np.newaxis,:,:]

        # square root of density matrix
        self.sqrt_den = np.sqrt(den)
        
        return U


    def _update_bathy(self, dist):
        """ Load bathymetry at the specified distance from the source.

            This updates the attributes 
            
                * depth
                * gradient
                * height_above_seafloor

            Args:
                dist: float
                    Distance from source in meters
        """
        if self.verbose:
            print('Updating bathymetry at {0:.2f} m'.format(dist))

        Nz = self.grid.Nz

        x = self.costheta * dist
        y = self.sintheta * dist

        depth = self.ocean.bathy(lon=x, lat=y)
        depth *= (-1.)

        if self.ignore_bathy_gradient:
            gradient = np.zeros(x.shape)

        else:
            dfdx = self.ocean.bathy_deriv(lon=x, lat=y, axis='lon')
            dfdy = self.ocean.bathy_deriv(lon=x, lat=y, axis='lat')
            gradient = self.costheta * dfdx + self.sintheta * dfdy
            gradient *= (-1.)

        depth = depth[np.newaxis,:]
        gradient = gradient[np.newaxis,:]

        self.depth = np.ones((Nz,1)) * depth
        self.gradient = np.ones((Nz,1)) * gradient

        self.height_above_seafloor = np.abs(self.grid.Z) - self.depth


    def _update_c(self, dist):
        """ Load sound speed at the specified distance from the source.

            This updates the attributes
            
                * n2w

            Args:
                dist: float
                    Distance from source in meters

        """
        if self.verbose:
            print('Updating sound speed at {0:.2f} m'.format(dist))

        # sound speed
        x = self.costheta * dist
        y = self.sintheta * dist
        z = self.grid.z[self.grid.z <= 0]
        x, _ = np.meshgrid(x,z)
        y, z = np.meshgrid(y,z) 
        x = x.flatten()
        y = y.flatten()
        z = z.flatten()
        c = self.c.eval(x=x, y=y, z=-z)        

        # refractive index squared
        self.n2_w[self.grid.Z_below] = (self.c0 / c)**2

        # mirror sound-speed profile above/below sea surface
        self.n2_w = self.grid.mirror(self.n2_w)


    def _save_output(self, step_no, dist, psi):
        """ Post-processe and collect output data at specified distance
            from the source.
            
            Args:
                step_no: int
                    Step number
                dist: float
                    Radial distance from the source in meters
                psi: 2d numpy array
                    Sound pressure field at the specified radial distance. 
                    Has shape (Nz,Nq) where Nz and Nq are the number of 
                    vertical and angular grid points, respectively.
        """
        if step_no > 0:
            dz = self.grid.dz
            idx = np.squeeze(np.round(self.receiver_depth/dz).astype(int))
            if np.ndim(idx) == 0:
                idx = np.array([idx])

            A = np.matmul(self.ifft_kernel, psi) 
            B = self.sqrt_den[idx][np.newaxis,:,:]

            self.field_horiz[:,:,:,step_no] = A * B * np.exp(1j * self.k0 * dist) / np.sqrt(dist)

            if self.vertical_slice:
                psi = np.fft.ifft(psi, axis=1) * np.exp(1j * self.k0 * dist) / np.sqrt(dist) * self.sqrt_den

        else:
            if self.vertical_slice:
                psi = np.fft.ifft(psi, axis=1)

        if self.vertical_slice:
            n = int(self.grid.Nz / 2)
            self.field_vert[:,:,step_no,:] = psi[:, :n, :]


    def _init_output(self, num_source_depth, receiver_depth, vertical_slice):

        Nr = self.grid.Nr
        Nq = self.grid.Nq
        Nz = self.grid.Nz
        Ns = num_source_depth
        Nd = len(receiver_depth)

        self.vertical_slice = vertical_slice

        self.receiver_depth = receiver_depth[:, np.newaxis]

        # ifft kernel
        kz = self.grid.kz
        self.ifft_kernel = np.exp(1j * np.matmul(self.receiver_depth, kz[np.newaxis,:])) / len(kz)
        self.ifft_kernel = self.ifft_kernel[np.newaxis,:,:]

        # output arrays
        self.field_horiz = np.empty(shape=(Ns, Nd, Nq, Nr), dtype=complex)  # sound intensity values
        self.field_vert = np.empty(shape=(Ns, int(Nz/2), Nr, Nq), dtype=complex)
