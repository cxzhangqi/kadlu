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

""" Parabolic Equation propagator module within the kadlu library

    This module provides an implementation of the Parabolic Equation 
    propagation scheme for numerically solving the wave equation 
    for sound pressure.

    Contents:
        TransmissionLossCalculator class:
        PEGrid class 
"""

import numpy as np
from numpy.lib import scimath
from kadlu.transmission_loss.pe_starter import PEStarter
from kadlu.transmission_loss.pe_propagator import PEPropagator
from kadlu.transmission_loss.environment_input import EnvironmentInput
import math


class TransmissionLossCalculator():
    """ Compute the reduction in intensity (transmission loss) of 
        sound waves propagating in an underwater environment.

        The wave equation is solved numerically using the 
        parabolic equation method.

        The computation is performed on a regular, cylindrical 
        grid with axes r, q, z:

            r: radial distance in meters
            q: azimuthal angle in radians
            z: vertical depth in meters

        Args:
            bathymetry: BathyInterpolator
                Interpolated bathymetry data
            flat_seafloor_depth: float
                Depth of flat seafloor in meters. Useful for testing purposes. 
                If flat_seafloor_depth is specified, the bathymetry input 
                argument will be ignored. 
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
            bathy_update: int
                How often the bathymetry data is updated. 
                If for example bathy_update=3, the bathymetry is updated at every 3rd step 
                of the propagation algorithm.
            sound_speed_update: int
                How often the sound-speed data is updated. 
                If for example sound_speed_update=3, the sound speed is updated at every 3rd step 
                of the propagation algorithm. By default sound_speed_update is set to infinity, 
                corresponding to a range-independent sound speed profile.
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
    def __init__(self, bathymetry=None, flat_seafloor_depth=None, sound_speed=None, ref_sound_speed=1500,\
            water_density=1.0, bottom_sound_speed=1700, bottom_loss=0.5, bottom_density=1.5,\
            step_size=None, range=50e3, angular_bin_size=1, vertical_bin_size=10, max_depth=12e3,\
            absorption_layer=1./6., starter_method='THOMSON', starter_aperture=88,\
            bathy_update=3, sound_speed_update=math.inf, verbose=False, progress_bar=True):

        self.bathymetry = bathymetry
        self.flat_seafloor_depth = flat_seafloor_depth
        self.sound_speed = sound_speed

        self.c0 = ref_sound_speed
        self.rhow = water_density
        self.cb = bottom_sound_speed
        self.bloss = bottom_loss
        self.rhob = bottom_density

        self.step_size = step_size
        self.range = range
        self.angular_bin_size = angular_bin_size
        self.vertical_bin_size = vertical_bin_size
        self.max_depth = max_depth
        self.absorption_layer = absorption_layer

        self.starter_method = starter_method
        self.starter_aperture = starter_aperture

        self.ndx_ChangeWD = max(1, bathy_update)
        self.ndx_ChangeNSQ = max(1, sound_speed_update)
        if sound_speed is None:
            self.ndx_ChangeNSQ = math.inf

        self.verbose = verbose
        self.progress_bar = progress_bar

        self.nsq = 1 # = (self.c0 / self.sound_speed)^2  refractive index squared

        if self.verbose:
            print('\nTransmission loss calculator successfully initialized')
            print('Bathymetry will be updated every {0} steps'.format(self.ndx_ChangeWD))
            if self.ndx_ChangeNSQ is math.inf:
                print('Adopting range-independent sound-speed profile')
            else:
                print('Sound speed will be updated every {0} steps'.format(self.ndx_ChangeNSQ))
                

    def run(self, frequency, source_depth, receiver_depths=[.1], vertical_slice=True,\
            ignore_bathy_gradient=False):
        """ Compute the transmission loss at the specified frequency, source depth, 
            and receiver depths.
            
            The transmission loss is computed at every grid point on 
            one or several horizontal planes at the specified receiver depth(s).

            Args:
                frequency: float
                    Sound frequency in Hz
                source_depth: float
                    Source depth in meters
                receiver_depths: list of floats
                    Depths of receivers. 
                vertical_slice: bool
                    Compute the transmission loss at all grid points on 
                    a vertical plane intersecting the source position. 
                    Note: This will slow down the computation.
                ignore_bathy_gradient: bool
                    Set the bathymetry gradient to zero everywhere.
                    This can be used to speed up the interpolation of the bathymetry 
                    data if the depth only changes gradually, implying that the gradient 
                    can be ignored.
        """
        if self.verbose:
            import time
            start = time.time()

        # frequency in Hz
        freq = frequency  
        
        # source position in meters
        xs = 0
        ys = 0
        zs = source_depth

        # reference wavelength and wavenumber
        lambda0 = self.c0 / freq        
        k0 = 2 * np.pi * freq / self.c0  

        smoothing_length_rho = self.c0 / freq / 4
        smoothing_length_ssp = np.finfo(float).eps  # machine epsilon

        # radial step size
        if self.step_size is None:
            dr = lambda0 / 2
        else:
            dr = self.step_size

        # azimuthal step size
        azimuthal_step = self.angular_bin_size / 180 * np.pi

        # vertical range
        vertical_range = 2 * self.max_depth * (1. + self.absorption_layer)
        vertical_step = self.vertical_bin_size

        # construct the computational grid
        grid = PEGrid(radial_step=dr, radial_range=self.range,\
                azimuthal_step=azimuthal_step, azimuthal_range=2*np.pi,\
                vertical_step=vertical_step, vertical_range=vertical_range)

        if self.verbose:
            print('Computational grid:')
            print('range (r) x angle (q) x depth (z)')
            print('Dimensions: {0} x {1} x {2}'.format(grid.Nr, grid.Nq, grid.Nz))
            print('min(r) = {0:.1f}, max(r) = {1:.1f}, delta(r) = {2:.1f}'.format(np.min(grid.r), np.max(grid.r), grid.dr))
            print('min(q) = {0:.1f}, max(q) = {1:.1f}, delta(q) = {2:.1f}'.format(np.min(grid.q)*180./np.pi, np.max(grid.q)*180./np.pi, grid.dq*180./np.pi))
            print('min(z) = {0:.1f}, max(z) = {1:.1f}, delta(z) = {2:.1f}'.format(np.min(grid.z), np.max(grid.z), grid.dz))
 
        # PE starter
        starter = PEStarter(ref_wavenumber=k0, grid=grid, method=self.starter_method, aperture=self.starter_aperture)

        # compute initial field
        psi = starter.eval(zs) * np.ones(shape=(1,grid.Nq))

        if self.verbose:
            print('Initial field computed')

        # module handling updates of environmental input
        env_input = EnvironmentInput(ref_wavenumber=k0, grid=grid, xs=xs, ys=ys,\
            freq=freq, ndx_ChangeWD=self.ndx_ChangeWD, ndx_ChangeNSQ=self.ndx_ChangeNSQ, c0=self.c0,\
            cb=self.cb, bloss=self.bloss, rhob=self.rhob, rhow=self.rhow,\
            smoothing_length_ssp=smoothing_length_ssp, smoothing_length_rho=smoothing_length_rho,
            absorption_layer=self.absorption_layer,\
            bathymetry=self.bathymetry, flat_seafloor_depth=self.flat_seafloor_depth, ignore_bathy_gradient=ignore_bathy_gradient)

        # Configure the PE propagator
        propagator = PEPropagator(ref_wavenumber=k0, grid=grid, env_input=env_input,\
                                verbose=self.verbose, progress_bar=self.progress_bar)

        # propagate
        output = propagator.run(psi=psi, depths=receiver_depths, vertical_slice=vertical_slice)

        # sound presure in dB
        sound_pressure = np.fft.fftshift(output.Ez[:,:,1:], axes=1)
        sound_pressure = 20 * np.log10(np.abs(sound_pressure))
        sound_pressure = np.squeeze(sound_pressure)

        if self.verbose:
            end = time.time()
            print('Calculation completed in {0:.2f} seconds'.format(end - start))

        return sound_pressure


        # ------- plotting ------- #

        # take only first 1/2 of z axis (?)
        z = grid.z[:int(grid.Nz/2)]

        # rearrange y axis (azimuthal) so values are increasing order
        y = np.fft.fftshift(grid.q)

        import matplotlib.pyplot as plt
        plot_r = grid.r[1:]
        plot_theta = np.fft.fftshift(np.squeeze(output.Ez_y))


        if False:
            R, TH = np.meshgrid(plot_r, plot_theta)
            XX = R * np.cos(TH)
            YY = R * np.sin(TH)
            fig = plt.contourf(XX, YY, SPfield, 100)
            plt.colorbar(fig)
            plt.show()

        if False:
            fig=plt.figure()
            ZZ = 20 * np.log10(np.abs(output.Af[:,:]))
            XX, YY = np.meshgrid(grid.r, z)
            fig = plt.contourf(XX, YY, ZZ, 100)
            plt.colorbar(fig)
            ax = plt.gca()
            ax.invert_yaxis()
            plt.show()



class PEGrid():
    """ Grid for Parabolic Equation solver.

        Creates a regular cylindrical grid, where

            r: radial distance
            q: azimuthal angle
            z: vertical depth

        The radial (r) and vertical coordinates are in meters
        while the angular coordinate (q) is in radians.

        Args:
            radial_step: float
                Radial step size
            radial_range: float
                Radial range
            azimuthal_step: float
                Angular bin size
            azimuthal_range: float
                Angular domain
            vertical_step: float
                Vertical bin size
            vertical_range: float
                Vertical range

        Attributes:
            r: 1d numpy array
                Radial grid values
            dr: float
                Radial step size
            Nr: int
                Number of radial steps
            q: 1d numpy array
                Angular grid values
            dq: float
                Angular bin size
            Nq: int
                Number of angular bins
            z: 1d numpy array
                Vertical grid values
            dz: float
                Vertical bin size
            Nz: int
                Number of vertical bins
            z: 1d numpy array
                Wavenumber grid values
            Q: 2d numpy array
                Angular values of azimuthal-vertical matrix; has shape (Nz,Nq).
            Z: 2d numpy array
                Vertical values of azimuthal-vertical matrix; has shape (Nz,Nq).
    """
    def __init__(self, radial_step, radial_range,\
            azimuthal_step, azimuthal_range, vertical_step, vertical_range):

        # radial
        self.r, self.dr, self.Nr = self.__make_radial_grid__(radial_step, radial_range)

        # azimuthal
        self.q, self.dq, self.Nq = self.__make_azimuthal_grid__(azimuthal_step, azimuthal_range)

        # vertical
        self.z, self.dz, self.Nz = self.__make_vertical_grid__(vertical_step, vertical_range)

        # wavenumber
        L = self.Nz * self.dz  
        self.kz = self.z * 2 * np.pi / (L * self.dz)

        # mesh-grid
        self.Q, self.Z = np.meshgrid(self.q, self.z)

    def __make_radial_grid__(self, dr, rmax):
        N = round(rmax / dr)
        r = np.arange(N+1, dtype=float)
        r *= dr
        return r, dr, N

    def __make_azimuthal_grid__(self, dq, qmax):
        N = int(np.ceil(qmax / dq))
        if N%2 == 1: 
            N = N + 1 # ensure even number of angular bins

        q_pos = np.arange(start=0, stop=N/2, step=1, dtype=float)
        q_neg = np.arange(start=-N/2, stop=0, step=1, dtype=float)
        q = np.concatenate((q_pos, q_neg)) 
        q *= dq
        return q, dq, N

    def __make_vertical_grid__(self, dz, zmax):
        N = zmax / dz
        N = round(N / 2) * 2  # ensure even number of vertical bins
        z_pos = np.arange(start=0, stop=N/2, step=1, dtype=float)
        z_neg = np.arange(start=-N/2, stop=0, step=1, dtype=float)
        z = np.concatenate((z_pos, z_neg))
        z *= dz
        return z, dz, N