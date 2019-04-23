# ================================================================================ #
#   Authors: Fabio Frazao and Oliver Kirsebom                                      #
#   Contact: fsfrazao@dal.ca, oliver.kirsebom@dal.ca                               #
#   Organization: MERIDIAN (https://meridian.cs.dal.ca/)                           #
#   Team: Data Analytics                                                           #
#   Project: ketos                                                                 #
#   Project goal: The ketos library provides functionalities for handling          #
#   and processing acoustic data and applying deep neural networks to sound        #
#   detection and classification tasks.                                            #
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

""" Parabolic Equation starter module within the kadlu library

    This module provides various starting fields for the 
    Parabolic Equation propagation scheme.

    Contents:
        PEStarter class: 
"""

import numpy as np
from numpy.lib import scimath
from enum import Enum
from kadlu.utils import get_member


class PEStarterMethod(Enum):
    """ Enum class for PE starter methods
    """
    GAUSSIAN = 1
    GREENE = 2
    THOMSON = 3


class PEStarter():
    """ Computes starting field for Parabolic-Equation propagator.
        
        Args:
            ref_wavenumber: float or numpy array
                Reference wavenumber in inverse meters
            grid: PEGrid
                Computational grid
            method: str
                Options are: GAUSSIAN, GREENE, THOMSON
            aperture: float
                Aperture in degrees

        Example:
            >>> from kadlu.transmission_loss.transmission_loss_calculator import PEGrid
            >>> from kadlu.transmission_loss.pe_starter import PEStarter
            >>>
            >>> # Create a regular azimuthal grid with depth and range of 1 km
            >>> # and grid spacing of 100 meters and 10 degrees
            >>> grid = PEGrid(100., 1000., 10.*np.pi/180., 2.*np.pi, 100., 1000.)
            >>>
            >>> # Initialize a Thomson PE starter with an aperture of 86 degrees
            >>> starter = PEStarter(ref_wavenumber=0.04, grid=grid, aperture=86)
            >>>
            >>> # Compute the initial field for a source depth of 9 meters
            >>> psi = starter.eval(zs=9)
            >>> psi = np.round(psi, 4) # round to 4 decimals
            >>> print(psi)
            [[ 0.    +0.j    ]
             [ 0.0101-0.0101j]
             [ 0.0205-0.0205j]
             [ 0.0319-0.0319j]
             [ 0.0451-0.0451j]
             [ 0.    +0.j    ]
             [-0.0451+0.0451j]
             [-0.0319+0.0319j]
             [-0.0205+0.0205j]
             [-0.0101+0.0101j]]
    """
    def __init__(self, ref_wavenumber, grid, method='THOMSON', aperture=90):

        self.method = get_member(PEStarterMethod, method)
        self.aperture = aperture
        self.k0 = ref_wavenumber
        self.grid = grid

    def eval(self, zs):
        """ Evaluate starting field for PE propagator at specified source depth
            
            Args:
                zs: float
                    Source depth in meters

            Returns:
                psi: numpy array
                    Initial sound pressure field along the vertical 
                    axis at zero range. Has shape (Nz,1) where Nz is 
                    the number of vertical grid points.
        """
        if self.method is PEStarterMethod.GAUSSIAN:
            psi = self.__gaussian__(zs)
        elif self.method is PEStarterMethod.GREENE:
            psi = self.__greene__(zs)
        elif self.method is PEStarterMethod.THOMSON:
            psi = self.__thomson__(zs)
            
        return psi

    def __gaussian__(self, zs):
        """ Gaussian starter
            
            Args:
                zs: float
                    Source depth in meters

            Returns:
                psi: 2d numpy array
                    Initial sound pressure field
        """
        k0 = self.k0
        Z = self.grid.Z

        # compute psi
        psi = np.sqrt(k0) * np.exp( -0.5*k0**2 *(Z-zs)**2 )
        psi = psi - ( np.sqrt(k0) * np.exp( -0.5*k0**2 *(Z+zs)**2 ))
        psi = np.fft.fft(psi)

        return psi

    def __greene__(self, zs):
        """ Greene's starter
            
            Args:
                zs: float
                    Source depth in meters

            Returns:
                psi: 2d numpy array
                    Initial sound pressure field
        """
        a = 1.4467
        b = .04201
        c = 3.0512
        k0 = self.k0
        Z = self.grid.Z

        # compute psi        
        psi = np.sqrt(k0) * (a - b * k0**2 * (Z - zs)**2) * np.exp(-(k0**2 * (Z - zs)**2) / c )
        psi = psi - (np.sqrt(k0) * (a - b * k0**2 * (Z + zs)**2) * np.exp(-(k0**2 * (Z + zs)**2) / c ))
        psi = np.fft.fft(psi)

        return psi

    def __thomson__(self, zs):
        """ Thomson's starter
            
            Args:
                zs: float
                    Source depth in meters

            Returns:
                psi: 2d numpy array
                    Initial sound pressure field
        """
        k0 = self.k0
        kz = self.grid.kz
        dz = self.grid.dz
        Nz = self.grid.Nz

        # compute psi
        psi = np.exp(-1j * np.pi / 4.) * 2 * scimath.sqrt(2 * np.pi) * np.sin(kz * zs) / scimath.sqrt(scimath.sqrt(k0**2 - kz**2))

        # normalize the starter
        psi = psi / dz
        psi[int(Nz/2)] = 0 

        # taper the spectrum to obtain desired angle using Turkey window
        kcut1 = k0 * np.sin(self.aperture / 180 * np.pi) 
        kcut0 = k0 * np.sin((self.aperture - 1.5) / 180 * np.pi)
        W = 0.5 * (1 + np.cos(np.pi / (kcut1 - kcut0) * (np.abs(kz) - kcut0)))
        W[np.abs(kz) >= kcut1] = 0
        W[np.abs(kz) <= kcut0] = 1
        psi = psi * W
        psi[np.abs(kz) >= kcut1] = 0
        psi = psi[:, np.newaxis]

        return psi

