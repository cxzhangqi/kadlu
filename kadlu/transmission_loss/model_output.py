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
                For all angular bins, compute the sound pressure on a
                vertical plane intersecting the source position. 
                Note: This will slow down the computation.

        Attributes:
            k0: float
                Reference wavenumber in inverse meters
            grid: PEGrid
                Computational grid
            env_input: EnviromentInput
                Environmental data
            vertical_slice: bool
                For all angular bins, compute the sound pressure on a
                vertical plane intersecting the source position. 
                Note: This will slow down the computation.
            counter: int
                Counter that keeps track of the radial step number. 
                Increments by 1 unit for every call to the collect method.
            field_vert: 3d numpy array
                Vertical slices of the sound pressure field. Has shape 
                (int(Nz/2), Nr+1, Nq) where Nz,Nr,Nq are the number of 
                grid points in the vertical, radial and azimuthal direction.
            depths: numpy array
                Receiver depths in meters
            field_horiz: 3d numpy array
                Horizontal slices of the sound pressure field at the specified 
                depths. Has shape (Nd, Nq, Nr+1) where Nd is the number of depths, 
                while Nq and Nr are the number of azimuthal and radial grid points, 
                respectively.
            ifft_kernel: numpy array
                Matrix used to compute the horizontal field.

        Example:
    """
    def __init__(self, ref_wavenumber, grid, env_input, depths=[.1], vertical_slice=True):

        self.k0 = ref_wavenumber
        self.grid = grid
        self.vertical_slice = vertical_slice
        self.env_input = env_input

        Nr = self.grid.Nr
        Nq = self.grid.Nq
        Nz = self.grid.Nz
        Nd = len(depths)

        self.counter = 0

        # vertical slices
        self.field_vert = np.empty(shape=(int(Nz/2), Nr+1, Nq), dtype=complex)

        # horizontal slices
        self.depths = np.array(depths)
        self.depths = self.depths[:, np.newaxis]
        self.field_horiz = np.empty(shape=(Nd, Nq, Nr+1), dtype=complex)  # sound intensity values

        # ifft kernel
        kz = self.grid.kz
        self.ifft_kernel = np.exp(1j * np.matmul(self.depths, kz[np.newaxis,:])) / len(kz)


    def collect(self, dist, psi):
        """ Post-processe and collect output data at specified distance
            from the source.
            
            Args:
                dist: float
                    Radial distance from the source in meters
                psi: 2d numpy array
                    Sound pressure field at the specified radial distance. 
                    Has shape (Nz,Nq) where Nz and Nq are the number of 
                    vertical and angular grid points, respectively.
        """
        self.counter += 1

        if dist != 0:
            dz = self.grid.dz
            idx = np.squeeze(np.round(self.depths/dz).astype(int))
            if np.ndim(idx) == 0:
                idx = np.array([idx])

            sqrt_denin = self.env_input.sqrt_denin

            A = np.matmul(self.ifft_kernel, psi)
            B = sqrt_denin[idx]

            self.field_horiz[:,:,self.counter-1] = A * B * np.exp(1j * self.k0 * dist) / np.sqrt(dist)

            if self.vertical_slice:
                psi = np.fft.ifft(psi, axis=0) * np.exp(1j * self.k0 * dist) / np.sqrt(dist) * sqrt_denin

        else:
            if self.vertical_slice:
                psi = np.fft.ifft(psi, axis=0)

        if self.vertical_slice:
            n = int(self.grid.Nz / 2)
            self.field_vert[:, self.counter-1, :] = psi[:n, :]
