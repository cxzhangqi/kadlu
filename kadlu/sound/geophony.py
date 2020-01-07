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

""" Geophony module within the kadlu library
"""
import numpy as np
from kadlu.geospatial.ocean import Ocean
from kadlu.sound.sound_propagation import TLCalculator, Seafloor
from kadlu.utils import xdist, ydist, XYtoLL
from tqdm import tqdm
from scipy.interpolate import interp1d
from datetime import datetime


class Geophony():
    """ Geophony modeling on a regular 3D grid.

        TODO: Check that the specified region is within the coverage of 
              of the ocean data. 

        Args:
            tl_calculator: instance of TLCalculator
                Transmission loss calculator
            depth: float or 1d array
                Depth(s) at which the noise level is computed.
            xy_res: float
                Horizontal spacing (in meters) between points at which the 
                noise level is computed. If None is specified, the spacing 
                will be set equal to sqrt(2) times the range of the transmission
                loss calculator.
            progress_bar: bool
                Display calculation progress bar. Default is True.

        Attributes:
            tl: instance of TLCalculator
                Transmission loss calculator
            south, north: float
                ymin, ymax coordinate boundaries to fetch bathymetry. range: -90, 90
            west, east: float
                xmin, xmax coordinate boundaries to fetch bathymetry. range: -180, 180
            depth: float or 1d array
                Depth(s) at which the noise level is computed.
            xy_res: float
                Horizontal spacing (in meters) between points at which the 
                noise level is computed. If None is specified, the spacing 
                will be set equal to sqrt(2) times the range of the transmission
                loss calculator.
            progress_bar: bool
                Display calculation progress bar. Default is True.
    """
    def __init__(self, tl_calculator, south, north, west, east, depth, xy_res=None, progress_bar=True):

        self.tl = tl_calculator

        if isinstance(depth, list):
            depth = np.array(depth)
        elif isinstance(depth, float) or isinstance(depth, int):
            depth = np.array([depth])

        self.depth = np.sort(depth)

        if xy_res is None:
            self.xy_res = np.sqrt(2) * self.tl.range['r']
        else:
            self.xy_res = xy_res

        self.tl.progress_bar = False
        self.progress_bar = progress_bar

        # prepare grid
        self.lats, self.lons, self.x, self.y = self._create_grid(south, north, west, east)
        self.bathy = self.tl.ocean.bathy(x=self.lons, y=self.lats, geometry='spherical')

        # wind source level interpolation table
        self._kewley1990 = interp1d(x=[2.57, 5.14, 10.29, 15.23, 20.58], y=[34, 39, 48, 53, 58], kind='linear', fill_value="extrapolate")

    def compute(self, frequency, below_seafloor=False, start=None, end=None):
        """ Compute the noise level within a specified geographic 
            region at a specified date and time.

            If below_seafloor is False (default), the noise level is only computed 
            at grid points above the seafloor, and is set to NaN below.

            TODO: Allow user to specify time, as an alternative to 
            start and end.

            Args:
                frequency: float
                    Sound frequency in Hz.
                below_seafloor: bool
                    Whether to compute the noise below the seafloor. Default is False.

            Returns:
                SPL: 3d numpy array
                    Sound pressure levels, has shape (Nx,Ny,Nz) where Nx is the number 
                    of west-east (longitude) grid points, Ny is the number of south-north 
                    (latitude) grid points, and Nz is the number of depths.
        """
        N = len(self.lats)
        SPL = None

        for i in tqdm(range(N), disable = not self.progress_bar):

            lat = self.lats[i]
            lon = self.lons[i]

            if below_seafloor: # include depths below seafloor
                depth = self.depth

            else: # ignore depths below seafloor
                seafloor_depth = -self.bathy[i]
                depth = self.depth[self.depth <= seafloor_depth] 

            if len(depth) == 0:
                dB = np.empty((1,len(self.depth)), dtype=float)
                dB[:,:] = np.nan

            else:
                # set receiver depth to 1/4 of the characteristic wave length
                receiver_depth = 0.25 * self.tl.c0 / frequency

                # load data and compute transmission loss
                self.tl.run(frequency=frequency, source_lat=lat, source_lon=lon,
                        source_depth=depth, receiver_depth=receiver_depth)

                TL = self.tl.TL[:,0,:,:]

                # source level
                SL = self._source_level(freq=frequency, grid=self.tl.grid, start=start, end=end)

                # integrate SL-TL to obtain sound pressure level
                p = np.power(10, (SL + TL) / 20)
                p = np.squeeze(np.apply_over_axes(np.sum, p, range(1, p.ndim))) # sum over all but the first axis
                dB = 20 * np.log10(p)

                if np.ndim(dB) == 0:
                    dB = np.array([dB])

                # pad, if necessary
                n = len(self.depth) - len(dB)
                if n > 0:
                    pad = np.empty(n)
                    pad[:] = np.nan
                    dB = np.concatenate((dB, pad))

                dB = dB[np.newaxis, :]

            if SPL is None:
                SPL = dB
            else:
                SPL = np.concatenate((SPL, dB), axis=0)


        SPL = np.reshape(SPL, newshape=(len(self.x), len(self.y), SPL.shape[1]))
            
        return SPL


    def _source_level(self, freq, grid, start, end, method='wind'):

        # x,y coordinates
        r = grid.r[1:]
        q = grid.q
        r, q = np.meshgrid(r, q)
        x = r * np.cos(q)
        y = r * np.sin(q)

        # area elements (m^2)
        a = grid.dr * grid.dq * r

        # flatten arrays
        x = x.flatten()
        y = y.flatten()
        a = a.flatten()

        # wind source level impirical parametrization
        assert method == 'wind', 'The only allowed method is wind'

        if method == 'wind':
            SL = self._wind_source_level_per_area(freq=freq, x=x, y=y, start=start, end=end)
            SL += 20 * np.log10(a)

        SL = np.reshape(SL, newshape=r.shape)
        SL = SL[np.newaxis, :, :]

        return SL


    def _wind_source_level_per_area(self, freq, x, y, start, end):

        # interpolate wind speed
        wind_speed = self.tl.ocean.wave(x=x, y=y) 

        # use parametrization of Kewley 1990
        # (Ocean Ambient Noise p. 114)
        WSL = self._kewley1990(x=wind_speed)

        return WSL


    def _create_grid(self, south, north, west, east):

        # select latitude closest to the equator
        if np.abs(south) < np.abs(north):
            lat = south
        else:
            lat = north

        # compute x and y range
        xd = xdist(lon2=east, lon1=west, lat=lat)
        xd += 2 * self.tl.range['r']
        yd = ydist(lat2=north, lat1=south) 
        yd += 2 * self.tl.range['r']

        # number of bins
        nx = int(xd / self.xy_res)
        ny = int(yd / self.xy_res)
        nx += nx%2 
        ny += ny%2 

        # create x and y arrays
        x = np.arange(start=-nx/2, stop=nx/2+1)
        x *= self.xy_res
        y = np.arange(start=-ny/2, stop=ny/2+1)
        y *= self.xy_res

        # convert to lat-lon
        lat_ref = 0.5 * (north + south)
        lon_ref = 0.5 * (east + west)
        lats, lons = XYtoLL(x=x, y=y, lat_ref=lat_ref, lon_ref=lon_ref, grid=True)

        lats = lats.flatten()
        lons = lons.flatten()

        return lats, lons, x, y
