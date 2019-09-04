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


class Geophony():

    def __init__(self, tl_calculator, depth, xy_res=None):

        self.tl = tl_calculator
        self.depth = depth

        if xy_res is None:
            self.xy_res = 2 * self.tl.range['r']
        else:
            self.xy_res = xy_res


    def model(self, frequency, south, north, west, east, time=None):

        lats, lons = self._create_grid(south, north, west, east)

        for lat,lon in zip(lats, lons):

            # transmission loss
            self.tl.run(frequency=frequency, source_lat=lat, source_lon=lon, source_depth=self.depth, time=time)          
            TL = self.tl.TL

            # source levels
            r = self.tl.grid.r[1:]
            q = self.tl.grid.q
            r, q = np.meshgrid(r, q)
            x = r * np.cos(q)
            y = r * np.sin(q)
            x = x.flatten()
            y = y.flatten()
            SL = 200 * self.tl.ocean.wave(x=x, y=y)  # 200 dB times wave height
            SL = np.reshape(SL, newshape=r.shape)

            # integrate SL-TL to obtain sound pressure level
            p = np.power(10, (SL + TL) / 20)
            p = np.squeeze(np.apply_over_axes(np.sum, p, range(1, p.ndim))) # sum over all but first axis
            SPL = 20 * np.log10(p)
            
        return SPL


    def _create_grid(self, south, north, west, east):

        # TODO: implement this method

        lats = np.array([45])
        lons = np.array([45])

        return lats, lons