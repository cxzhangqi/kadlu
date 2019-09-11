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


class Geophony():

    def __init__(self, tl_calculator, depth, xy_res=None, time_dependent_tl=False, progress_bar=True):

        self.tl = tl_calculator
        self.tl.progress_bar = False
        self.progress_bar = progress_bar

        self.depth = depth
        self.time_dependent_tl = time_dependent_tl

        if xy_res is None:
            self.xy_res = np.sqrt(2) * self.tl.range['r']
        else:
            self.xy_res = xy_res


    def model(self, frequency, south, north, west, east, start=None, stop=None):

        lats, lons, x, y = self._create_grid(south, north, west, east)

        SPL = None

        N = len(lats)
        for i in tqdm(range(N), disable = not self.progress_bar):

            lat = lats[i]
            lon = lons[i]

            # loop over times
            # compute transmission loss once at the center time, or for ever time step
            time = None

            # load data and compute transmission loss
            self.tl.run(frequency=frequency, source_lat=lat, source_lon=lon, source_depth=self.depth, time=time)          
            TL = self.tl.TL

            # source level
            SL = self._source_level(lat, lon, time)

            # integrate SL-TL to obtain sound pressure level
            p = np.power(10, (SL + TL) / 20)
            p = np.squeeze(np.apply_over_axes(np.sum, p, range(1, p.ndim))) # sum over all but first axis
            dB = 20 * np.log10(p)
            dB = dB[np.newaxis, :]

            if SPL is None:
                SPL = dB
            else:
                SPL = np.concatenate((SPL, dB), axis=0)

        SPL = np.reshape(SPL, newshape=(len(x), len(y), SPL.shape[1]))
            
        return SPL, x, y


    def _source_level(self, lat, lon, time):

        # TODO: complete implementation of this method

        # load wave data for the specified time step
        # (if not already loaded via the TL calculation)

        r = self.tl.grid.r[1:]
        q = self.tl.grid.q
        r, q = np.meshgrid(r, q)
        x = r * np.cos(q)
        y = r * np.sin(q)
        x = x.flatten()
        y = y.flatten()
        SL = 200 * self.tl.ocean.wave(x=x, y=y)  # 200 dB times wave height
        SL = np.reshape(SL, newshape=r.shape)

        return SL


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