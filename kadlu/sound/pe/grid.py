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

""" Parabolic Equation grid module within the kadlu library

    This module generates the computational grid for the 
    Parabolic Equation propagation scheme.

    Contents:
        Grid class: 
"""
import numpy as np
from numpy.lib import scimath
from enum import Enum
from kadlu.utils import get_member


class Grid():
    """ Grid for Parabolic Equation solver.

        Creates a regular grid using a cylindrical 
        r,q,z coordinate system, where

            r: radial distance
            q: azimuthal angle
            z: vertical depth

        The radial (r) and vertical coordinates are in meters
        while the angular coordinate (q) is in radians.

        Args:
            dr: float
                Radial step size
            rmax: float
                Radial range
            dq: float
                Angular bin size
            qmax: float
                Angular domain
            dz: float
                Vertical bin size
            zmax: float
                Vertical range

        Attributes:
            r: 1d numpy array
                Radial values
            dr: float
                Radial step size
            Nr: int
                Number of radial steps
            q: 1d numpy array
                Angular values
            dq: float
                Angular bin size
            Nq: int
                Number of angular bins
            z: 1d numpy array
                Vertical values
            dz: float
                Vertical bin size
            Nz: int
                Number of vertical bins
            kz: 1d numpy array
                Wavenumber grid values
            Q: 2d numpy array
                Angular values of q,z meshgrid; has shape (Nz,Nq).
            Z: 2d numpy array
                Vertical values of q,z meshgrid; has shape (Nz,Nq).
    """
    def __init__(self, dr, rmax, dq, qmax, dz, zmax):

        # radial
        self.r, self.dr, self.Nr = self._radial_coordinates(dr, rmax)

        # azimuthal
        self.q, self.dq, self.Nq = self._azimuthal_coordinates(dq, qmax)

        # vertical
        self.z, self.dz, self.Nz, self._indices_below, self._indices_above, self._indices_mirror = self._vertical_coordinates(dz, zmax)

        # wavenumber
        L = self.Nz * self.dz  
        self.kz = self.z * 2 * np.pi / (L * self.dz)

        # meshgrids
        self.Q, self.Z = np.meshgrid(self.q, self.z)

        # z-q grid points with negative z (i.e. below sea surface)
        self.Z_below = np.nonzero(self.Z <= 0)


    def _radial_coordinates(self, dr, rmax):
        N = int(round(rmax / dr))
        r = np.arange(N+1, dtype=float)
        r *= dr
        return r, dr, len(r)


    def _azimuthal_coordinates(self, dq, qmax):        
        N = int(np.ceil(qmax / dq))
        if N%2 == 1: 
            N = N + 1 # ensure even number of angular bins

        q_pos = np.arange(start=0, stop=N/2, step=1, dtype=float)
        q_neg = np.arange(start=-N/2, stop=0, step=1, dtype=float)
        q = np.concatenate((q_pos, q_neg)) 
        q *= dq
        return q, dq, N


    def _vertical_coordinates(self, dz, zmax):
        N = 2 * zmax / dz
        N = int(round(N / 2) * 2)  # ensure even number of vertical bins
        z_pos = np.arange(start=0, stop=N/2, step=1, dtype=float)
        z_neg = np.arange(start=-N/2, stop=0, step=1, dtype=float)
        z = np.concatenate((z_pos, z_neg))
        z *= dz
        # z indices below and above surface
        one = np.array([0], dtype=int)
        below = np.arange(start=N-1, step=-1, stop=N/2-1, dtype=int)
        below = np.concatenate([one, below])
        above = np.arange(start=1, step=1, stop=N/2, dtype=int)
        mirror = below[1:-1]
        return z, dz, N, below, above, mirror


    def mirror(self, a, axis=0):

        if axis != 0:
            a = np.swapaxes(a, 0, axis)

        a[self._indices_above] = a[self._indices_mirror]

        if axis != 0:
            a = np.swapaxes(a, 0, axis)

        return a

