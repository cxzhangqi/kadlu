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

""" Unit tests for the sound speed module within the kadlu library
"""

import pytest
import os
import numpy as np
from kadlu.sound.sound_speed import SoundSpeed
from kadlu.geospatial.ocean import Ocean

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)),"assets")


def test_ocean_or_ssp_must_be_specified():
    with pytest.raises(AssertionError):
        SoundSpeed(num_depths=50, rel_err=None)

def test_sound_speed_from_uniform_data():
    # environment data provider
    o = Ocean(default=False, cache=False,
              load_bathymetry=-1000, load_temp=4, load_salinity=3,
              south=44, north=45, west=60, east=61)
    # instance of sound speed class 
    _ = SoundSpeed(o, num_depths=50, rel_err=None)

def test_sound_speed_from_uniform_ssp():
    # instance of sound speed class 
    ss = SoundSpeed(ssp=1499, num_depths=50, rel_err=None)
    # evaluate
    c = ss.eval(x=0, y=0, z=0)
    assert c == 1499

def test_sound_speed_from_ssp():
    # sound speed profile
    z0 = np.array([0, -10, -20, -30, -60])
    c0 = np.array([1500, 1510, 1512, 1599, 1489])
    # instance of sound speed class 
    ss = SoundSpeed(ssp=(c0,z0), num_depths=50, rel_err=None)
    # evaluate
    c = ss.eval(x=0, y=0, z=z0, grid=True)
    assert np.all(np.abs(c-c0) < 1E-6)

def test_query_sound_speed_interpolation_data():
    # sound speed profile
    z0 = np.array([0, -10, -20, -30, -60])
    c0 = np.array([1500, 1510, 1512, 1599, 1489])
    # instance of sound speed class 
    ss = SoundSpeed(ssp=(c0,z0), num_depths=50, rel_err=None)
    # query underlying data
    d = ss.eval()
    assert np.all(np.abs(d[1]-z0) < 1E-6)
    assert np.all(np.abs(d[0]-c0) < 1E-6)
