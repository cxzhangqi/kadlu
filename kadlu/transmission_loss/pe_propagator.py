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
        PEPropagator class: 
"""

import numpy as np
from tqdm import tqdm
from numpy.lib import scimath
from kadlu.transmission_loss.model_output import OutputCollector


class PEPropagator():
    """ Propagates the sound pressure field from zero range to the boundary 
        of the computational domain using a parabolic-equation numerical scheme.
        
        Args:
            ref_wavenumber: float
                Reference wavenumber in inverse meters
            grid: PEGrid
                Computational grid
            env_input: EnviromentInput
                Environmental data
            verbose: bool
                Print information during execution
            progress_bar: bool
                Show progress bar. Only shown if verbose if False.            

        Attributes:
            k0: float
                Reference wavenumber in inverse meters
            grid: PEGrid
                Computational grid
            env_input: EnviromentInput
                Environmental data
            verbose: bool
                Print information during execution
            progress_bar: bool
                Show progress bar. Only shown if verbose if False.            

        Example:
    """
    def __init__(self, ref_wavenumber, grid, env_input, verbose=False, progress_bar=True):

        self.k0 = ref_wavenumber
        self.grid = grid
        self.env_input = env_input
        self.verbose = verbose
        if self.verbose:
            self.progress_bar = False
        else:
            self.progress_bar = progress_bar


    def run(self, psi, depths=[.1], vertical_slice=True):
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
                depths: list of floats
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
        # output collector
        output = OutputCollector(ref_wavenumber=self.k0, grid=self.grid, env_input=self.env_input,\
            vertical_slice=vertical_slice, depths=depths)

        # initial Nx2D free propagator
        fr_half, fr_full = self.__free_propagator__()

        # output field at 0
        output.collect(dist=0, psi=psi)

        # PE marching starts here
        dist = 0
        Nr = self.grid.Nr
        dr = self.grid.dr
        for _ in tqdm(range(1, Nr+1), disable = not self.progress_bar):
    
            # (1) r --> r + dr/2 free propagation
            psi = fr_half * psi

            # (2) phase adjustment at r + dr/2
            U = self.env_input.update(dist + dr/2)            
            psi = np.fft.fft(U * np.fft.ifft(psi, axis=0), axis=0)    

            # (3) r + dr/2 --> r + dr free propagation
            psi = fr_half * psi

            # increment distance
            dist = dist + dr

            # collect output
            output.collect(dist=dist, psi=psi)

        return output


    def __free_propagator__(self):
        """ Compure matrices used to propagate the sound pressure field 
            in the radial direction by half or the entire grid spacing.

            The matrices have shape (Nz,Nz) where Nz and Nq are the number 
            of vertical and angular grid points, respectively.

            Returns:
                fr_half: 2d numpy array
                    Half-step free propagation matrix. 
                fr_full: 2d numpy array
                    Full-step free propagation matrix.
        """
        k0 = self.k0
        dr = self.grid.dr
        kz = self.grid.kz
        ny = self.grid.Nq

        fr_half = np.exp(1j * dr / 2 * (scimath.sqrt(k0**2 - kz**2) - k0))
        fr_half = fr_half[:,np.newaxis]
        fr_half = fr_half * np.ones(shape=(1,ny))

        fr_full = np.exp(1j * dr * (scimath.sqrt(k0**2 - kz**2) - k0)) 
        fr_full = fr_full[:,np.newaxis]        
        fr_full = fr_full * np.ones(shape=(1,ny))

        return fr_half, fr_full


