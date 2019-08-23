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
import gsw
import numpy as np
from numpy.lib import scimath
from kadlu.utils import XYtoLL

class EnvironmentInput():
    """ Compute the reduction in intensity (transmission loss) of 

        Args:
            ref_wavenumber: float
                Reference wavenumber in inverse meters
            grid: PEGrid
                Computational grid
            xs: float
                Source x coordinate in meters
            ys: float
                Source y coordinate in meters
            freq: float
                Frequency in Hz
            steps_btw_bathy_updates: int
                How often the bathymetry data is updated. If for example steps_btw_bathy_updates=3, the 
                bathymetry is updated at every 3rd step.
            steps_btw_sound_speed_updates: int
                How often the sound-speed data is updated. If for example steps_btw_sound_speed_updates=3, 
                the sound speed is updated at every 3rd step. By default steps_btw_sound_speed_updates is 
                set to infinity, corresponding to a range-independent sound speed profile.
            

        Attributes:

        Example:
    """

    def __init__(self, ref_wavenumber, grid, xs, ys, freq,\
            steps_btw_bathy_updates, steps_btw_sound_speed_updates,\
            c0, cb, bottom_loss, bottom_density, water_density,\
            smoothing_length_sound_speed, smoothing_length_density,\
            absorption_layer, env_data=None, ignore_bathy_gradient=False,\
            flat_seafloor_depth=None, sound_speed=None, verbose=False):

        assert env_data or flat_seafloor_depth, 'env_data or flat_seafloor_depth must be specified'

        self.env_data = env_data
        self.ignore_bathy_gradient = ignore_bathy_gradient
        self.flat_seafloor_depth = flat_seafloor_depth
        self.sound_speed = sound_speed
        self.verbose = verbose

        self.k0 = ref_wavenumber
        self.c0 = c0

        self.xs = xs
        self.ys = ys

        self.smoothing_length_sound_speed = smoothing_length_sound_speed
        self.smoothing_length_density = smoothing_length_density

        self.grid = grid
        self.Z = grid.Z
        self.dr = grid.dr

        self._range_independent = (grid.Nr == 1) and (grid.Nq == 1) # range independent environment

        m = grid.Z.shape
        n = grid.Nq

        # allocate memory

        # bathymetry
        self._unchanged = np.zeros(shape=n, dtype=bool)
        self.depth_old = np.empty(shape=(1,n))
        self.depth = np.empty(m)
        self.gradient = np.empty(m)
        self.height_above_seafloor = np.empty(m)

        # refractive index squared of water, n_w^2
        self.n2w = np.zeros(m)
        self.n2w_new = np.zeros(m)
        
        # derived real quantities
        self.H_c = np.empty(m, dtype=float)
        self.H_rho = np.empty(m, dtype=float)
        self.ddenin = np.empty(m, dtype=float)   
        self.d2denin = np.empty(m, dtype=float)
        self.denin = np.empty(m, dtype=float)
        self.sqrt_denin = np.empty(m, dtype=float)

        # derived complex quantities
        self.n2in = np.empty(m, dtype=complex)
        self.U = np.zeros(m, dtype=complex)

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

        # water density
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

            if self.verbose:
                print('Computing U matrix {0:.2f} m'.format(dist))

        return self.U


    def __update_env_model__(self, dist):
        """ Load the environment data at the specified distance from the source.

            Then update the environment model, by computing the attributes

                * H_c
                * n2in
                * H_rho
                * denin
                * sqrt_denin
                * ddenin
                * d2denin

            Args:
                dist: float
                    Distance from source in meters

            Returns:
                new_any: bool
                    True, if one or more angular bins bins experience a change in the environment attributes.
                new: 1d numpy array
                    Indices of those angular bins that experience a change in the environment attributes.
        """
        new = self.__update_env_input__(dist)[0]
        new_any = np.any(new)

        if new_any:

            # smooth ssp
            self.H_c[:,new] = (1 + np.tanh(self.height_above_seafloor[:,new] / self.smoothing_length_sound_speed / 2)) / 2
            self.n2in[:,new] = self.n2w[:,new] + (self.n2b - self.n2w[:,new]) * self.H_c[:,new]
            itmp = (self.depth[0,:] == 0) 
            if np.any(itmp):
                self.n2in[0,itmp] = self.n2b
            
            # smooth density
            _tanh = np.tanh(self.height_above_seafloor[:,new] / self.smoothing_length_density / 2)
            self.H_rho[:,new] = (1 + _tanh) / 2
            self.denin[:,new] = self.water_density + (self.bottom_density - self.water_density) * self.H_rho[:,new]
            self.sqrt_denin[:,new] = np.sqrt(self.denin[:,new])
            
            _sech2 = 1 / np.cosh(self.height_above_seafloor[:,new] / self.smoothing_length_density / 2)
            _sech2 = _sech2 * _sech2
            self.ddenin[:,new] =  _sech2 / self.smoothing_length_density / 2 * np.sqrt(1 + self.gradient[:,new]**2)
            self.ddenin[:,new] =  (self.bottom_density - self.water_density) / 2 * self.ddenin[:,new]

            self.d2denin[:,new] = -_sech2 / self.smoothing_length_density / 2 * (_tanh / self.smoothing_length_density * (1 + self.gradient[:,new]**2))
            self.d2denin[:,new] =  (self.bottom_density - self.water_density) / 2 * self.d2denin[:,new]

        return new_any, new


    def __update_env_input__(self, dist):
        """ Load the environment data at the specified distance from the source.

            Args:
                dist: float
                    Distance from source in meters

            Returns:
                new: 1d numpy array
                    Indices of those angular bins that experience a change in 
                    bathymetry and/or sound speed.
        """
        new = np.logical_or(self.__update_bathy__(dist), self.__update_sound_speed__(dist))
        return new


    def __update_bathy__(self, dist):
        """ Load bathymetry at the specified distance from the source.

            This updates the attributes 
            
                * depth
                * gradient
                * height_above_seafloor
                * dist_next_bathy_update
                * depth_old

            Args:
                dist: float
                    Distance from source in meters

            Returns:
                new: 1d numpy array
                    Indices of those angular bins that experience a change in 
                    bathymetry.
        """
        new = self._unchanged

        if (dist == self.dr/2) or (dist >= self.dist_next_bathy_update):
            
            if self.verbose:
                print('Updating bathymetry at {0:.2f} m'.format(dist))

            # next distance at which to update bathymetry
            self.dist_next_bathy_update = dist + self.dn_bathy * self.dr
            
            # get bathymetry
            depth_new, gradient_new = self.__seafloor_depth__(dist)

            # which bins have changed bathymetry?
            new = np.logical_and(self.depth_old != depth_new, np.logical_not(np.isnan(depth_new)))

            # update 'old' attribute
            self.depth_old = depth_new

            # number of vertical bins
            Nz = self.grid.Nz

            # mask those entries that haven't changed
            self.depth[:,new[0]] = np.ones((Nz,1)) * depth_new[new]
            self.gradient[:,new[0]] = np.ones((Nz,1)) * gradient_new[new]

            # update height above seafloor matrix
            self.height_above_seafloor[:,new[0]] = np.abs(self.Z[:,new[0]]) - self.depth[:,new[0]]
            
        new = np.squeeze(new)

        return new


    def __update_sound_speed__(self, dist):
        """ Load sound speed at the specified distance from the source.

            This updates the attributes
            
                * n2w
                * n2w_new
                * dist_next_sound_speed_update

            Args:
                dist: float
                    Distance from source in meters

            Returns:
                new: 1d numpy array
                    Indices of those angular bins that experience a change in 
                    sound speed.
        """
        new = self._unchanged

        if (dist == self.dr/2) or ((dist >= self.dist_next_sound_speed_update) and not self._range_independent):

            if self.verbose:
                print('Updating sound speed at {0:.2f} m'.format(dist))

            # next distance at which to update sound speed
            self.dist_next_sound_speed_update = dist + self.dn_sound_speed * self.dr

            self.n2w_new = np.copy(self.n2w)

            # z-q grid points with negative z (i.e. below sea surface)
            below_sea_surf = np.nonzero(self.Z <= 0)

            # load refractive index squared
            self.n2w_new[below_sea_surf] = self.__refractive_index__(dist, self.grid.z[self.grid.z <= 0])

            # number of vertical bins
            Nz = self.grid.Nz

            # z indeces below sea surface
            one = np.array([0], dtype=int)
            indeces = np.arange(start=Nz-1, step=-1, stop=Nz/2, dtype=int)
            indeces = np.concatenate([one, indeces])

            # which z-q grid points have changed refractive index?
            new = np.any(self.n2w[indeces,:] - self.n2w_new[indeces,:] != 0, axis=0)

            # mirror sound-speed profile above/below sea surface
            indeces2 = np.arange(start=1, step=1, stop=Nz/2, dtype=int)
            self.n2w_new[np.ix_(indeces2, new)] = self.n2w_new[np.ix_(indeces[1:], new)]

            # update 'current' attribute
            self.n2w[:,new] = self.n2w_new[:,new]       

        return new


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
                depth: 2d numpy array
                    Depth at each point. Has shape (1,Nq) where Nq is the number of angular bins.
                gradient: 2d numpy array
                    Gradient at each point. Has shape (1,Nq) where Nq is the number of angular bins.
        """
        x = self.xs + self.costheta * dist
        y = self.ys + self.sintheta * dist

        if self.flat_seafloor_depth:
            n = len(x)
            depth = self.flat_seafloor_depth * np.ones(n)
            gradient = np.zeros(n)

        else:
            depth = self.env_data.bathy(x=x, y=y)
            depth *= (-1.)

            if self.ignore_bathy_gradient:
                gradient = np.zeros(n)
            else:            
                dfdx = self.env_data.bathy_gradient(x=x, y=y, axis='x')
                dfdy = self.env_data.bathy_gradient(x=x, y=y, axis='y')
                gradient = self.costheta * dfdx + self.sintheta * dfdy
                gradient *= (-1.)
            
        depth = depth[np.newaxis,:]
        gradient = gradient[np.newaxis,:]

        return depth, gradient


    def __refractive_index__(self, dist, z):
        """ Compute refractive index squared. 

            The computation is performed at equally spaced points on a circle, 
            at the specified distance from the source.

            Args:
                dist: float
                    Distance from source in meters
                z: numpy array
                    Depths

            Returns:
                nsq: 2d numpy array
                    Refractive index squared. Has shape (Nz,Nq) where Nz and Nq 
                    are the number of vertical and angular bins, respectively.
        """
        x = self.xs + self.costheta * dist
        y = self.ys + self.sintheta * dist

        x, _ = np.meshgrid(x,z)
        y, z = np.meshgrid(y,z) 

        x = x.flatten()
        y = y.flatten()
        z = z.flatten()

        if self.env_data is None:
            c = np.ones(z.shape) * self.sound_speed

        else:
            t = self.env_data.temp(x=x, y=y, z=-z)  # in-situ temperature
            SP = self.env_data.salinity(x=x, y=y, z=-z)  # practical salinity
            lat_ref = self.env_data.origin.latitude
            lon_ref = self.env_data.origin.longitude
            lats, lons = XYtoLL(x=x, y=y, lat_ref=lat_ref, lon_ref=lon_ref)
            c = self._sound_speed(lats=lats, lons=lons, z=z, t=t, SP=SP)
        
        nsq = (self.c0 / c)**2

        return nsq


    def _sound_speed(self, lats, lons, z, t, SP):
        """ Compute sound speed.

            Args:
                lats: numpy array
                    Latitudes (-90 to 90 degrees)
                lons: numpy array
                    Longitudes (-180 to 180 degrees)
                z: numpy array
                    Depths (meters)
                t: numpy array
                    In-situ temperature (Celsius)
                SP: numpy array
                    Practical Salinity (psu)

            Returns:
                c: numpy array
                    Sound speed (m/s) 
        """
        p = gsw.p_from_z(z=z, lat=lats)  # sea pressure
        SA = gsw.SA_from_SP(SP, p, lons, lats)  # absolute salinity
        CT = gsw.CT_from_t(SA, t, p)  # conservative temperature
        c = gsw.density.sound_speed(SA=SA, CT=CT, p=p)
        return c


    def seafloor_depth_transect(self, dist, angle):
        """ Compute seafloor depth along a straight line from the source position.

            Args:
                dist: numpy array
                    Distance from source in meters
                angle: float
                    Angle in degrees

            Returns:
                depth: numpy array
                    Depth at each point. 
        """
        angle_rad = angle * np.pi / 180.
        x = self.xs + np.cos(angle_rad) * dist
        y = self.ys + np.sin(angle_rad) * dist

        if self.flat_seafloor_depth:
            n = len(x)
            depth = self.flat_seafloor_depth * np.ones(n)

        else:
            depth = self.env_data.bathy(x=x, y=y)
            depth *= (-1.)
            
        return depth
