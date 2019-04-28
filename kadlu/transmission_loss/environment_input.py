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

""" Environment input module within the kadlu library

    This module handles loading and pre-processing of environment data for the 
    the transmission-loss calculation.

    Contents:
        EnvironmentInput class
"""

import numpy as np
from numpy.lib import scimath

class EnvironmentInput():

    def __init__(self, ref_wavenumber, grid, xs, ys, freq,\
            steps_btw_bathy_updates, steps_btw_sound_speed_updates,\
            c0, cb, bottom_loss, bottom_density, water_density,\
            smoothing_length_sound_speed, smoothing_length_density,\
            absorption_layer, bathymetry=None, ignore_bathy_gradient=False,\
            flat_seafloor_depth=None, sound_speed=None, verbose=False):

        assert bathymetry or flat_seafloor_depth, 'bathymetry or flat_seafloor_depth must be provided'

        self.bathymetry = bathymetry
        self.ignore_bathy_gradient = ignore_bathy_gradient
        self.flat_seafloor_depth = flat_seafloor_depth
        self.sound_speed = sound_speed
        self.verbose = verbose

        self.k0 = ref_wavenumber

        self.xs = xs
        self.ys = ys

        self.smoothing_length_sound_speed = smoothing_length_sound_speed
        self.smoothing_length_density = smoothing_length_density

        self.Z = grid.Z
        self.dr = grid.dr

        self._range_independent = (grid.Nr == 1) and (grid.Nq == 1) # range independent environment

        # allocate memory
        m = grid.Z.shape
        n = grid.Nq

        self.U = np.zeros(shape=m, dtype=complex)
        self._unchanged = np.zeros(shape=n, dtype=bool)

        self.wd_old = np.empty(shape=(1,n))
        self.wd_mask = np.empty(m)
        self.DwdDy = np.empty(m)
        self.Z_sub_wd = np.empty(m)

        self.n2w = np.zeros(m)
        self.n2w_new = np.zeros(m)
        
        # updated in update_phase_screen()
        self.n2in = np.empty(m, dtype=complex)
        self.H_c = np.empty(m, dtype=float)
        self.H_rho = np.empty(m, dtype=float)
        self.ddenin = np.empty(m, dtype=float)   
        self.d2denin = np.empty(m, dtype=float)
        self.denin = np.empty(m, dtype=float)
        self.sqrt_denin = np.empty(m, dtype=float)

        # compute cos(theta) and sin(theta) for all angular bins
        self.costheta = np.cos(grid.q)
        self.sintheta = np.sin(grid.q)

        # how often is bathymetry and sound-speed data updated
        self.dn_bathy = max(1, steps_btw_bathy_updates)
        self.dn_sound_speed = max(1, steps_btw_sound_speed_updates)
        self.dist_next_bathy_update = 0
        self.dist_next_sound_speed_update = 0
        
        # complex bottom sound speed
        self.bottom_density = bottom_density
        kbi = bottom_loss / (cb / freq) / 20. / np.log10(np.e) 
        betab = kbi / 2 / np.pi / freq
        cbi = np.roots([betab, -1, betab*cb**2])  # roots of polynomial p[0]*x^n+...p[n]
        cbi = cbi[np.imag(cbi) == 0] 
        cbi = cbi[np.logical_and(cbi >= 0, cbi < cb)]
        cb = cb - 1j * cbi
        self.n2b = (c0 / cb)**2

        self.water_density = water_density

        # boudary conditions in vertical dimension (z)
        absorption_coefficient =  1. / np.log10(np.e) / np.pi
        absorption_layer_thickness = absorption_layer * np.max(np.abs(grid.z))
        D = absorption_layer_thickness / 3
        self.attenuation = 1j * absorption_coefficient * np.exp(-(np.abs(grid.Z) - np.max(np.abs(grid.z)))**2 / D**2)


    def update(self, dist):        
        """ Update the environment model and compute the propagation matrix, U, 
            at the specified distance.
            
            Args:
                dist: float
                    Radial distance from the source in meters

            Returns:
                self.U: 2d numpy array
                    Propagation matrix. Has shape (Nz,Nq) where Nz and Nq are the 
                    number of vertical and angular grid points, respectively.
        """
        do_update, indices = self.__update_env_model__(dist)        
        if do_update:
            self.U[:, indices] = np.exp(1j * self.dr * self.k0 * (-1 + scimath.sqrt( self.n2in[:,indices] + self.attenuation[:,indices] +\
                1/2 / self.k0**2 * (self.d2denin[:,indices] / self.denin[:,indices] - 3/2 * (self.ddenin[:,indices] / self.denin[:,indices])**2))))

            print('Computing U matrix {0:.2f} m'.format(dist))

        return self.U


    def __update_env_model__(self, dist):

        new_env = self.__update_env_input__(dist)[0]
        new_env_any = np.any(new_env)

        if new_env_any:

            # smooth ssp
            self.H_c[:,new_env] = (1 + np.tanh(self.Z_sub_wd[:,new_env] / self.smoothing_length_sound_speed / 2)) / 2
            self.n2in[:,new_env] = self.n2w[:,new_env] + (self.n2b - self.n2w[:,new_env]) * self.H_c[:,new_env]
            itmp = (self.wd_mask[0,:] == 0) 
            if np.any(itmp):
                self.n2in[0,itmp] = self.n2b
            
            # smooth density
            TANH = np.tanh(self.Z_sub_wd[:,new_env] / self.smoothing_length_density / 2)
            self.H_rho[:,new_env] = (1 + TANH) / 2
            self.denin[:,new_env] = self.water_density + (self.bottom_density - self.water_density) * self.H_rho[:,new_env]
            self.sqrt_denin[:,new_env] = np.sqrt(self.denin[:,new_env])
            
            SECH2 = 1 / np.cosh(self.Z_sub_wd[:,new_env] / self.smoothing_length_density / 2)
            SECH2 = SECH2 * SECH2
            self.ddenin[:,new_env] =  SECH2 / self.smoothing_length_density / 2 * np.sqrt(1 + self.DwdDy[:,new_env]**2)
            self.ddenin[:,new_env] =  (self.bottom_density - self.water_density) / 2 * self.ddenin[:,new_env]

            self.d2denin[:,new_env] = -SECH2 / self.smoothing_length_density / 2 * (TANH / self.smoothing_length_density * (1 + self.DwdDy[:,new_env]**2))
            self.d2denin[:,new_env] =  (self.bottom_density - self.water_density) / 2 * self.d2denin[:,new_env]

        return new_env_any, new_env


    def __update_env_input__(self, dist):
        """ Update the environment data (i.e. bathymetry and sound speed) at 
            the specified distance from the source.

            Args:
                dist: float
                    Distance from source in meters

            Returns:
                new: 1d numpy array
                    Indices of those entries that experience a change in 
                    bathymetry and/or sound speed. 
        """
        print(self.__update_bathy__(dist))
        print(self.__update_refractive_index__(dist))


        new = np.logical_or(self.__update_bathy__(dist), self.__update_refractive_index__(dist))
        return new


    def __update_bathy__(self, dist):

        new_bathy = self._unchanged

        if (dist == self.dr/2) or (dist >= self.dist_next_bathy_update):
            
            # next distance at which to update bathymetry
            self.dist_next_bathy_update = dist + self.dn_bathy * self.dr
            
            # get bathymetry
            wd_new, DwdDy_new = self.__seafloor_depth__(dist)
            wd_new = wd_new[np.newaxis,:]
            DwdDy_new = DwdDy_new[np.newaxis,:]

            new_bathy = np.logical_and(self.wd_old != wd_new, np.logical_not(np.isnan(wd_new)))

            self.wd_old = wd_new

            k = self.Z.shape

            self.wd_mask[:,new_bathy[0]] = np.ones((k[0],1)) * wd_new[new_bathy]  # water depth mask
            self.DwdDy[:,new_bathy[0]] = np.ones((k[0],1)) * DwdDy_new[new_bathy]
            self.Z_sub_wd[:,new_bathy[0]] = np.abs(self.Z[:,new_bathy[0]]) - self.wd_mask[:,new_bathy[0]]
            
            print('Updating bathymetry at {0:.2f} m'.format(dist))

        return new_bathy


    def __update_refractive_index__(self, dist):

        new_refr = self._unchanged

        if (dist == self.dr/2) or ((dist >= self.dist_next_sound_speed_update) and not self._range_independent):  # update water column

            # next distance at which to update water column
            self.dist_next_sound_speed_update = dist + self.dn_sound_speed * self.dr

            # interpolate water nsq
            self.n2w_new = np.copy(self.n2w)
            IDZ = (self.Z <= 0)
            idy = np.nonzero(IDZ)[1]
            IDZ = np.nonzero(IDZ)
            self.n2w_new[IDZ] = self.__refractive_index__(dist, IDZ)   #sub_NSQ(x(idy).',y(idy).',Z(IDZ));

            k = self.Z.shape
            one = np.array([0], dtype=int)
            indeces = np.arange(start=k[0]-1, step=-1, stop=k[0]/2, dtype=int)
            indeces = np.concatenate([one, indeces])

            new_refr = np.any(self.n2w[indeces,:] - self.n2w_new[indeces,:] != 0, axis=0)

            indeces2 = np.arange(start=1, step=1, stop=k[0]/2, dtype=int)

            self.n2w_new[np.ix_(indeces2, new_refr)] = self.n2w_new[np.ix_(indeces[1:], new_refr)]

            self.n2w[:,new_refr] = self.n2w_new[:,new_refr]       

        return new_refr


    def __seafloor_depth__(self, dist):
        """ Compute seafloor depth and gradient.

            The computation is performed at equally spaced points on a circle, 
            at the specified distance from the source.

            The depth is positive below the sea surface, positive above.

            The gradient is computed in the direction perpendicular to the circle.

            Args:
                dist: float
                    Distance from source in meters

            Returns:
                depth: 1d numpy array
                    Depth at each point
                gradient: 1d numpy array
                    Gradient at each point
        """
        x = self.xs + self.costheta * dist
        y = self.ys + self.sintheta * dist

        if self.flat_seafloor_depth:
            n = len(x)
            depth = self.flat_seafloor_depth * np.ones(n)
            gradient = np.zeros(n)

        else:
            depth = self.bathymetry.eval_xy(x=x, y=y)
            depth *= (-1.)

            if self.ignore_bathy_gradient:
                gradient = np.zeros(n)
            else:            
                dfdx = self.bathymetry.eval_xy(x=x, y=y, x_deriv_order=1)
                dfdy = self.bathymetry.eval_xy(x=x, y=y, y_deriv_order=1)
                gradient = self.costheta * dfdx + self.sintheta * dfdy
                gradient *= (-1.)
            
        return depth, gradient


    def __refractive_index__(self, dist, IDZ):
        """ Compute refractive index squared. 

            The computation is performed at equally spaced points on a circle, 
            at the specified distance from the source.

            Args:
                dist: float
                    Distance from source in meters

            Returns:
                nsq: 2d numpy array
                    Refractive index squared. Has shape (Nz,Nq) where Nz and Nq 
                    are the number of vertical and angular bins, respectively.
        """
        print(' WARNING: Adopting uniform sound-speed profile')

        n = self.Z[IDZ].shape

        nsq = np.ones(n)

        return nsq