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

""" Sound propagation module within the kadlu library

    The sound pressure wave equation is solved numerically using a
    Parabolic Equation propagation scheme.

    Contents:
        TLCalculator class:
        TLGrid class 
"""
import math
import numpy as np
from numpy.lib import scimath
from kadlu.sound.sound_speed import SoundSpeed
from kadlu.sound.pe.grid import Grid
from kadlu.sound.pe.starter import Starter
from kadlu.sound.pe.propagator import Propagator
from kadlu.utils import XYtoLL

from sys import platform as sys_pf
if sys_pf == 'darwin' or sys_pf == 'win32':
    import matplotlib
    matplotlib.use("TkAgg")
from matplotlib import pyplot as plt


class Seafloor():
    """ Properties of the sea floor.

        Args:
            c: float
                Uniform and isotropic sound speed in m/s
            density: float
                Uniform density in g/cm^3
            thickness: float
                Thickness in meters
            loss: float
                Attenuation in dB/lambda, where lambda is the reference 
                wave length computed as c0 / frequency

        Attributes:
            c: float
                Uniform and isotropic sound speed in m/s
            density: float
                Uniform density in g/cm^3
            thickness: float
                Thickness in meters
            loss: float
                Attenuation in dB/lambda, where lambda is the reference 
                wave length computed as c0 / frequency
            c0: float
                Reference sound speed in m/s
            frequency: float
                Sound frequency in Hz
    """
    def __init__(self, c=1700, density=1.5, thickness=2000, loss=0.5):

        self.c = c
        self.density = density
        self.thickness = thickness
        self.loss = loss

        self.c0 = None
        self.frequency = None

    def nsq(self):
        """ Compute the refractive index squared.

            Produces an AssertionError if the reference sound speed (c0)
            and the sound frequency (frequency) have not been specified.

            Returns:
                n2: complex
                    Refractive index squared
        """
        assert self.frequency is not None, 'Frequency must be specified to allow calculation of the refractive index'
        assert self.c0 is not None, 'Reference sound speed must be specified to allow calculation of the refractive index'

        f = self.frequency
        ki = self.loss / (self.c / f) / 20. / np.log10(np.e) 
        beta = ki / 2 / np.pi / f
        ci = np.roots([beta, -1, beta * self.c**2])  # roots of polynomial p[0]*x^n+...p[n]
        ci = ci[np.imag(ci) == 0] 
        ci = ci[np.logical_and(ci >= 0, ci < self.c)]
        c = self.c - 1j * ci
        n2 = (self.c0 / c[0])**2

        return n2


class TLCalculator():
    """ Compute the reduction in intensity (transmission loss) of 
        sound waves propagating in an underwater environment.

        The wave equation is solved numerically using the 
        parabolic equation method.

        The computation is performed on a regular grid in a 
        cylindrical r,q,z coordinate system, where:

            r: radial distance in meters
            q: azimuthal angle in radians
            z: vertical depth in meters

        Args:
            env_data: DataProvider
                Environment data provider
            flat_seafloor_depth: float
                Depth of flat seafloor in meters. Useful for testing purposes. 
                Overwrites the bathymetry from `env_data`. 
            sound_speed: SoundSpeedInterpolator
                Interpolated sound-speed profile. If None is specified, a uniform 
                sound-speed profile equal to the reference sound speed will be assumed.
            ref_sound_speed: float
                Reference sound speed in meters/second.
            water_density: float 
                Water density in grams/cm^3
            bottom_sound_speed: float
                Homogenous bottom sound speed in meters/seconds
            bottom_loss: float
                Homogenous bottom attenuation in dB/lambda, where lambda is the reference 
                wave length given by lambda = ref_sound_speed / frequency
            bottom_density: float
                Homogenous bottom density in grams/cm^3
            step_size: float
                Radial step size in meters. If None is specified (default), the step size 
                is computed as lambda/2, where lambda = ref_sound_speed / frequency is 
                the reference wave length.
            range: float
                Radial range in meters
            angular_bin_size: float
                Angular bin size in degrees
            vertical_bin_size: float
                Vertical bin size in meters
            max_depth: float
                Maximum depth in meters. The vertical range used for the computation is  
                [-z_max,z_max], where z_max = max_depth * (1 + absorption_layer). 
            absorption_layer: float
                Thickness of the artificial absorption layer expressed as a fraction of the vertical range.
                For example, if the vertical range is 1.2 km (max_depth=1200) and absorption_layer=1./6. 
                (the default value), the thickness of the artificial absorption layer will be 200 meters.
            starter_method: str
                PE starter method. Options are: GAUSSIAN, GREENE, THOMSON
            starter_aperture: float
                Aperture of PE starter in degrees
            steps_btw_bathy_updates: int
                How often the bathymetry data is updated. If for example steps_btw_bathy_updates=3, the 
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
            grid: PEGrid
                Computational grid
            TL: numpy array
                Transmission loss in dB on a horizontal plane at the specified depths. 
                Has shape (Nd, Nq, Nr) where Nd is the number of depths, Nq is number of 
                angular bins, and Nr is the number of radial steps. If Nd=1, the first dimension
                is discarded. 
            TL_vertical: numpy array
                Transmission loss in dB on a vertical plane for all angular bins. 
                Has shape (Nz/2, Nr+1, Nq) where Nz is the number of vertical grid points, Nr is 
                the number of radial steps, and Nq is number of angular bins. 
            receiver_depth: list of floats
                Depths of receivers. 
            env_input: EnvironmentInput
                Module holding environment data

        Example:
    """
    def __init__(self, ocean, seafloor, sound_speed=None, ref_sound_speed=1500,\
            radial_bin=None, radial_range=50e3,\
            angular_bin=10, angular_range=2*np.pi,\
            vertical_bin=10, vertical_range=None,\
            absorption_layer=1/6.,\
            starter_method='THOMSON', starter_aperture=88,\
            steps_btw_bathy_updates=1, steps_btw_sound_speed_updates=1,\
            verbose=False, progress_bar=True):

        self.ocean = ocean
        self.seafloor = seafloor
        self.c0 = ref_sound_speed

        self.seafloor.c0 = self.c0

        self.steps_btw_bathy_updates = steps_btw_bathy_updates
        self.steps_btw_c_updates = steps_btw_sound_speed_updates

        if sound_speed is None:
            self._compute_sound_speed = True
            self.c = None
        else:
            self._compute_sound_speed = False
            self.c = SoundSpeed(ssp=sound_speed)
            self.steps_btw_c_updates = math.inf

        self.bin_size = {'r':radial_bin, 'q':angular_bin, 'z':vertical_bin}
        self.range = {'r':radial_range, 'q':angular_range, 'z':vertical_range}

        self.absorption_layer = absorption_layer

        self.starter_method = starter_method
        self.starter_aperture = starter_aperture

        self.TL = list()
        self.TLv = list()

        self.verbose = verbose
        self.progress_bar = progress_bar

        if self.verbose:
            print('\nTransmission loss calculator successfully initialized')
            print('Bathymetry will be updated every {0} steps'.format(self.steps_btw_bathy_updates))
            if self.steps_btw_c_updates is math.inf:
                print('Adopting range-independent sound-speed profile')
            else:
                print('Sound speed will be updated every {0} steps'.format(self.steps_btw_c_updates))


    def _update_source_location(self, lat, lon):

        r = self.range['r'] + 10e3
        S, N, W, E = self._bounding_box(lat=lat, lon=lon, r=r)

        self.ocean.load(south=S, north=N, west=W, east=E)

        if self._compute_sound_speed:
            self.c = SoundSpeed(self.ocean)


    def _bounding_box(self, lat, lon, r):

        angles = np.linspace(start=0, stop=2*np.pi, num=361)

        x = r * np.cos(angles)
        y = r * np.sin(angles)

        lat, lon = XYtoLL(x, y, lat_ref=lat, lon_ref=lon)

        S = np.min(lat)
        N = np.max(lat)
        W = np.min(lon)
        E = np.max(lon)

        return S, N, W, E


    def _create_grid(self, frequency):

        dz = self.bin_size['z']
        rmax = self.range['r']
        qmax = self.range['q']

        # automatic determination of radial step size
        if self.bin_size['r'] is None:
            dr = 0.5 * self.c0 / frequency
        else:
            dr = self.bin_size['r']

        # convert angular bin size to radians
        dq = self.bin_size['q'] / 180 * np.pi

        # automatic determination of vertical range
        max_depth = -np.min(self.ocean.bathy())
        zmax = (max_depth + self.seafloor.thickness) * (1. + self.absorption_layer)

        # create grid
        grid = Grid(dr=dr, rmax=rmax, dq=dq, qmax=qmax, dz=dz, zmax=zmax)

        if self.verbose:
            print('Created computational grid:')
            print('range (r) x angle (q) x depth (z)')
            print('Dimensions: {0} x {1} x {2}'.format(grid.Nr, grid.Nq, grid.Nz))
            print('min(r) = {0:.1f}, max(r) = {1:.1f}, delta(r) = {2:.1f}'.format(np.min(grid.r), np.max(grid.r), grid.dr))
            print('min(q) = {0:.1f}, max(q) = {1:.1f}, delta(q) = {2:.1f}'.format(np.min(grid.q)*180./np.pi, np.max(grid.q)*180./np.pi, grid.dq*180./np.pi))
            print('min(z) = {0:.1f}, max(z) = {1:.1f}, delta(z) = {2:.1f}'.format(np.min(grid.z), np.max(grid.z), grid.dz))

        return grid


    def run(self, frequency, source_lat, source_lon, source_depth,\
            receiver_depth=[.1], vertical_slice=False,\
            ignore_bathy_gradient=False, ignore_below_seafloor=True):
        """ Compute the transmission loss at the specified frequency, source depth, 
            and receiver depths.
            
            The transmission loss is computed at every grid point on 
            one or several horizontal planes at the specified receiver depth(s), 
            and optionally also on vertical planes for all angular bins.

            The results are stored in the attributes TL_dB and TL_dB_vert.

            Args:
                frequency: float
                    Sound frequency in Hz
                source_depth: float or list
                    Source depth in meters
                receiver_depth: list of floats
                    Depths of receivers. 
                vertical_slice: bool
                    For all angular bins, compute the sound pressure on a
                    vertical plane intersecting the source position. 
                    Note: This will slow down the computation.
                ignore_bathy_gradient: bool
                    Set the bathymetry gradient to zero everywhere.
                    This can be used to speed up the interpolation of the bathymetry 
                    data if the depth only changes gradually, implying that the gradient 
                    can be ignored.
        """
        source_depth = self._toarray(source_depth)
        receiver_depth = self._toarray(receiver_depth)

        if self.verbose:
            import time
            start = time.time()
            print('Begin transmission-loss calculation')
            print('Source depths:', source_depth)
            print('Receiver depths:', receiver_depth)
            if ignore_bathy_gradient:
                print('Ignoring bathymetry gradient')
            if vertical_slice:
                print('Computing the transmission loss on a vertical plane')
        
        self.seafloor.frequency = frequency

        # load data and initialize grid
        self._update_source_location(lat=source_lat, lon=source_lon)

        seafloor_depth = np.abs(self.ocean.bathy(x=0, y=0))

        # create grid
        grid = self._create_grid(frequency)

        # PE starter
        k0 = 2 * np.pi * frequency / self.c0
        starter = Starter(grid=grid, k0=k0, method=self.starter_method, aperture=self.starter_aperture)

        # loop over source depths
        # TODO: Vectorize the loop over source depths
        for zs in source_depth:

            # only compute if source is above seafloor
            if ignore_below_seafloor and (np.abs(zs) > seafloor_depth):
                self.TL.append(None)
                self.TLv.append(None)

            # compute initial field
            psi = starter.eval(zs=zs) * np.ones(shape=(1,grid.Nq))
            if self.verbose:
                print('Initial field computed')

            # smoothing lengths
            smooth_len_den = self.c0 / frequency / 4
            smooth_len_c = np.finfo(float).eps
 
            # Configure the PE propagator
            propagator = Propagator(ocean=self.ocean, seafloor=self.seafloor,\
                c=self.c, c0=self.c0, k0=k0, grid=grid,\
                smooth_len_den=smooth_len_den, smooth_len_c=smooth_len_c,\
                absorption_layer=self.absorption_layer,\
                ignore_bathy_gradient=ignore_bathy_gradient,\
                bathy_step=self.steps_btw_bathy_updates,\
                c_step=self.steps_btw_c_updates,\
                verbose=self.verbose, progress_bar=self.progress_bar)

            # propagate
            propagator.run(psi=psi, receiver_depth=receiver_depth, vertical_slice=vertical_slice)

            # output
            fh = propagator.field_horiz
            fv = propagator.field_vert

            # transmission loss in dB (horizontal plane)
            TL = np.fft.fftshift(fh[:,:,1:], axes=1)
            TL = 20 * np.log10(np.abs(TL))
            TL = np.squeeze(TL)

            # transmission loss in dB (vertical plane)  
            if vertical_slice:
                TLv = 20 * np.ma.log10(np.abs(fv))  # OBS: this computation is rather slow
                TLv = np.squeeze(TLv)
            else:
                TLv = None

            self.TL.append(TL)
            self.TLv.append(TLv)


        if self.verbose:
            end = time.time()
            print('Calculation completed in {0:.2f} seconds'.format(end - start))

        # store relevant data
        self.grid = grid
        self.source_depth = source_depth
        self.receiver_depth = receiver_depth


    def plot(self, source_depth_index=0, receiver_depth_index=0):
        """ Plot the transmission loss on a horizontal plane at fixed depth.

            The argument `depth_index` referes to the array `receiver_depth` 
            provided has input argument for the `run` method when the 
            transmission loss was calculated.

            Returns None if the transmission loss has not been computed.

            Args:
                depth_index: int
                    Depth index.

            Returns:
                fig: matplotlib.figure.Figure
                    A figure object.
        """
        assert source_depth_index < len(self.source_depth), 'Invalid source depth index. The index must be < {0}'.format(len(self.source_depth))

        assert receiver_depth_index < len(self.receiver_depth), 'Invalid receiver depth index. The index must be < {0}'.format(len(self.receiver_depth))

        TL = self.TL[source_depth_index]

        # check that transmission loss has been calculated
        if self.TL is None:
            print("Use the run method to calculate the transmission loss before plotting")
            return None

        if np.ndim(self.TL) == 2:
            tl = self.TL

        elif np.ndim(self.TL) == 3:
            tl = self.TL[receiver_depth_index,:,:]

        r = self.grid.r[1:]
        q = np.fft.fftshift(np.squeeze(self.grid.q))

        r, q = np.meshgrid(r, q)
        x = r * np.cos(q)
        y = r * np.sin(q)

        fig = plt.contourf(x, y, tl, 100)

        ax = plt.gca()
        ax.set_xlabel('x(m)')
        ax.set_ylabel('y (m)')
        plt.title('Transmission loss, source at {0:.1f} m, receiver at {1:.1f} m'.format(self.receiver_depth[receiver_depth_index], self.source_depth[source_depth_index]))

        plt.colorbar(fig, format='%+2.0f dB')

        return fig


    def plot_vertical(self, angle=0, source_depth_index=0, show_bathy=False):
        """ Plot the transmission loss on a vertical plane for a selected angular bin.

            Returns None if the transmission loss has not been computed.

            Args:
                angle: float
                    Angle in degrees with respect to the longitudinal axis
                show_bathy: bool
                    Display bathymetry superimposed on transmission loss

            Returns:
                fig: matplotlib.figure.Figure
                    A figure object.
        """
        assert source_depth_index < len(self.source_depth), 'Invalid source depth index. The index must be < {0}'.format(len(self.source_depth))

        TLv = self.TLv[source_depth_index]

        # check that transmission loss has been calculated
        if self.TLv is None:
            print("Use the run method with option vertical_slice set to True to calculate the transmission loss before plotting")
            return None

        z = self.grid.z[:int(self.grid.Nz/2)]
        x, y = np.meshgrid(self.grid.r, z)

        # find nearest angular bin
        if angle <= 180:
            q0 = angle * np.pi / 180
        else:
            q0 = (angle - 360) * np.pi / 180

        idx = np.abs(self.grid.q - q0).argmin()

        # transmission loss
        tl = np.squeeze(self.TL_vertical[:,:,idx])

        # bathy
        angle_rad = angle * np.pi / 180.
        x = np.cos(angle_rad) * self.grid.r
        y = np.sin(angle_rad) * self.grid.r
        bathy = self.ocean.bathy(x=x, y=y)
        bathy *= (-1.)

        # min and max transmission loss (excluding sea surface bin)
        tl_min = np.min(tl[1:,:])
        tl_max = np.max(tl[1:,:])

        fig = plt.figure()
        fig = plt.contourf(x, y, tl, 100, vmin=tl_min, vmax=tl_max)

        plt.colorbar(fig, format='%+2.0f dB')

        ax = plt.gca()
        ax.invert_yaxis()
        ax.set_xlabel('Range (m)')
        ax.set_ylabel('Depth (m)')
        angle = self.grid.q[idx] * 180 / np.pi
        plt.title('Transmission loss at {0:.2f} degrees'.format(angle))

        if show_bathy:
            plt.plot(self.grid.r, bathy, 'w')
            
        return fig


    def _toarray(self, v):
        if isinstance(v, float) or isinstance(v, int):
            v = [v]
        
        v = np.array(v)
        return v