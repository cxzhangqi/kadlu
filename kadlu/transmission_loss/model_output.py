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

""" Model output module within the kadlu library

    This module handles post-processing and collects 
    output data from the transmission-loss calculation.

    Contents:
        OutputCollector class
"""

import numpy as np

class OutputCollector():
    """ Post-processes and collects output data from the 
        transmission loss calculation.
        
        Args:
            ref_wavenumber: float
                Reference wavenumber in inverse meters
            grid: PEGrid
                Computational grid
            env_input: EnviromentInput
                Environmental data
            depths: list of floats
                Depths at which horizontal slices of the sound pressure 
                field are computed.
            vertical_slice: bool
                Compute the sound pressure at all grid points on 
                a vertical plane intersecting the source position. 
                Note: This will slow down the computation.

        Attributes:
            k0: float
                Reference wavenumber in inverse meters
            grid: PEGrid
                Computational grid
            env_input: EnviromentInput
                Environmental data

        Example:
    """

    def __init__(self, ref_wavenumber, grid, env_input, vertical_slice=True, depths=[.1]):

        self.k0 = ref_wavenumber
        self.grid = grid
        # sound intensity on the vertical x-z plane crossing the source
        self.vertical_slice = vertical_slice
        self.env_input = env_input

        Z = self.grid.Z
        Y = self.grid.Q
        kz = self.grid.kz
        x = self.grid.r 
        ny = self.grid.Nq
        nz = self.grid.Nz

        self.counter = 0
        self.Af = np.empty(shape=(int(nz/2), len(self.grid.r)), dtype=complex)


        # which points to be output?
        self.yso = np.argwhere(Y[0,:] >= 0)[0]
        self.nzhalf = int(Z.shape[0] / 2)

        # this are used for storing the calculated 3d field values
        self.Ez_z = np.array(depths) 
        self.Ez_z = self.Ez_z[:, np.newaxis] 
        self.Ez_y = Y[0,:]  # y values (azimuthal)
        self.Ez = np.empty(shape=(len(self.Ez_z), ny, len(x)), dtype=complex)  # sound intensity values

        self.Ez_ifft_kernel = np.exp(1j * np.matmul(self.Ez_z, kz[np.newaxis,:])) / len(kz)

        self.iout_Ez = 0


    def collect(self, dist, psi):

        self.iout_Ez += 1

        if dist != 0:
            dz = self.grid.dz

            idx = np.squeeze(np.round(self.Ez_z/dz).astype(int))
            
            if np.ndim(idx) == 0:
                idx = np.array([idx])

            sqrt_denin = self.env_input.sqrt_denin

            A = np.matmul(self.Ez_ifft_kernel, psi)
            B = sqrt_denin[idx]

            self.Ez[:,:,self.iout_Ez-1] = A * B * np.exp(1j * self.k0 * dist) / np.sqrt(dist)

            if self.vertical_slice:
                psi = np.fft.ifft(psi, axis=0) * np.exp(1j * self.k0 * dist) / np.sqrt(dist) * sqrt_denin

        else:
            if self.vertical_slice:
                psi = np.fft.ifft(psi, axis=0)

        if self.vertical_slice:
            self.Af[:, self.counter] = np.squeeze(psi[:self.nzhalf, self.yso])
            self.counter += 1
