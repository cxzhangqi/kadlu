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
from scipy.interpolate import interp1d, interp2d
from datetime import datetime
from scipy.interpolate import RectBivariateSpline


# wind source level tabulation of Kewley et al. 1990
_kewley_table = np.array([[40.0, 44.0, 48.0, 53.0, 58.0],\
                 [37.5, 42.5, 48.0, 53.0, 58.0],\
                 [34.0, 39.0, 48.0, 53.0, 58.0]])
_kewley_table = np.swapaxes(_kewley_table, 0, 1)
_kewley_interp = interp2d(x=[40, 100, 300], y=[2.57, 5.14, 10.29, 15.23, 20.58], z=_kewley_table, kind='linear')


def source_level_kewley(freq, wind_speed):
    """ Compute the wind source level according to the 
        tabulation of Kewley et al. 1990. (Ocean Ambient Noise p. 114).

        Values outside the tabulation domain are extrapolated via 
        nearest-neighbor extrapolation.

        Args:
            freq: float
                Frequency in Hz
            wind_speed: float or array
                Wind speed in m/s

        Returns:
            SL: float or array
                Source level in units of dB re 1 uPa^2 / Hz @ 1m / m^2
    """    
    sl = _kewley_interp(x=freq, y=wind_speed)
    sl = np.squeeze(sl)
    return sl

def source_level(freq, x, y, area, ocean, method, grid=False, geometry='planar'):
    """ Compute the source levels at the specified frequency and coordinates.
    
        Args:
            frequency: float
                Sound frequency in Hz.
            x: float or array
                x-coordinate(s) or longitude(s)
            y: float or array
                y-coordinate(s) or latitude(s)
            area: float or array
                Source area in units of meters squared. Must have same shape as x and y.
            ocean: instance of the Ocean class
                Ocean environmental data.
            method: str
                Method used to compute the source levels.
            grid: bool
                Specify how to combine elements of x and y. If x and y have different
                lengths, specifying grid has no effect as it is automatically set to True.
            geometry: str
                Can be either 'planar' (default) or 'spherical'.

        Returns:
            SL: float or 1d numpy array
                Source levels in units of dB re 1 uPa^2 / Hz @ 1m.
    """
    if method == 'Kewley':
        wind_speed = ocean.wave(x=x, y=y) 
        sl = source_level_kewley(freq=freq, wind_speed=wind_speed) # source level per unit area
        sl += 20 * np.log10(area) # scale by area

    return sl

def create_geophony_xy_grid(south, north, west, east, x_res, y_res):
    """ Create 2D surface grid for the geophony computation.

        Args:
            south, north: float
                Latitude range
            west, east: float
                Longitude range
            x_res: float
                Longitude resolution in meters
            y_res: float
                Latitude resolution in meters
    """
    # select latitude closest to the equator
    if np.abs(south) < np.abs(north):
        lat = south
    else:
        lat = north

    # compute x and y range
    xd = xdist(lon2=east, lon1=west, lat=lat)
    yd = ydist(lat2=north, lat1=south) 

    # number of bins
    nx = int(xd / x_res)
    ny = int(yd / y_res)
    nx += nx%2 
    ny += ny%2 

    # create x and y arrays
    x = np.arange(start=-nx/2, stop=nx/2+1)
    x *= x_res
    y = np.arange(start=-ny/2, stop=ny/2+1)
    y *= y_res

    # convert to lat-lon
    lat_ref = 0.5 * (north + south)
    lon_ref = 0.5 * (east + west)
    lats, lons = XYtoLL(x=x, y=y, lat_ref=lat_ref, lon_ref=lon_ref, grid=True)

    lats = lats.flatten()
    lons = lons.flatten()

    return lats, lons, x, y

class Geophony():
    """ Geophony modeling on a regular 3D grid.

        TODO: Check that the specified region is within the coverage of 
              of the ocean data. 

        Args:
            tl_calculator: instance of TLCalculator
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
            source_level_method: str
                Method used to compute the source levels.
            progress_bar: bool
                Display calculation progress bar. Default is True.

        Attributes:
            tl: instance of TLCalculator
                Transmission loss calculator
            depth: float, list, or 1d numpy array
                Depth(s) at which the noise level is computed.
            lats, lons: 1d numpy array
                Latitude and longitudes of the computational grid
            x, y: 1d numpy array
                Planar coordinates of the computational grid
            bathy: 1d numpy array
                Bathymetry at each grid point
            source_level_method: str
                Method used to compute the source levels.
            progress_bar: bool
                Display calculation progress bar. Default is True.
    """
    def __init__(self, tl_calculator, south, north, west, east, depth, xy_res=None, source_level_method='Kewley', progress_bar=True):

        self.tl = tl_calculator

        if isinstance(depth, list):
            depth = np.array(depth)
        elif isinstance(depth, float) or isinstance(depth, int):
            depth = np.array([depth])

        self.depth = np.sort(depth)

        if xy_res is None:
            xy_res = np.sqrt(2) * tl_calculator.range['r']
        else:
            xy_res = xy_res

        self.tl.progress_bar = False
        self.progress_bar = progress_bar

        assert source_level_method == 'Kewley', 'Invalid method for computing source levels'
        self.source_level_method = source_level_method

        # prepare grid
        self.lats, self.lons, self.x, self.y = create_geophony_xy_grid(south, north, west, east, x_res=xy_res, y_res=xy_res)

        # interpolate bathymetry
        self.bathy = self.tl.ocean.bathy(x=self.lons, y=self.lats, geometry='spherical')

    def compute(self, frequency, below_seafloor=False, start=None, end=None):
        """ Compute the noise level within a specified geographic 
            region at a specified date and time.

            If below_seafloor is False (default), the noise level is only computed 
            at grid points above the seafloor, and is set to NaN below.

            TODO: Allow user to specify a single date-time value, as an alternative to 
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
                SL = self._compute_source_level(freq=frequency)

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

        # transform output array to desired shape
        SPL = np.reshape(SPL, newshape=(len(self.y), len(self.x), SPL.shape[1]))
        SPL = np.swapaxes(SPL, 0, 1)
            
        return SPL

    def _compute_source_level(self, freq):
        """ Compute the source levels at the specified frequency at each surface 
            grid point.

            The source levels are computed according to the method prescribed 
            by the source_level_method attribute.

            Args:
                frequency: float
                    Sound frequency in Hz.

            Returns:
                SL: 3d numpy array
                    Source levels in units of dB re 1 uPa^2 / Hz @ 1m.
                    Has shape (1,Nq,Nr) where Nq is the number of angular bins 
                    and Nr is the number of radial bins in grid.
        """
        grid = self.tl.grid

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

        # compute source levels
        sl = source_level(freq=freq, x=x, y=y, area=a, ocean=self.tl.ocean, method=self.source_level_method)
        sl = np.reshape(sl, newshape=r.shape) # transform to desired shape
        sl = sl[np.newaxis, :, :]

        return sl