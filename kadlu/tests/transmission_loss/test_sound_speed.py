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
from kadlu.transmission_loss.sound_speed import SoundSpeed
from kadlu.geospatial.data_provider import DataProvider

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)),"assets")


def test_can_create_sound_speed_instance():
    # environment data provider
    env = DataProvider(bathy_source="CHS", south=43, west=-60, north=44, east=-59)
    # instance of sound speed class 
    c = SoundSpeed(env, num_depths=50, rel_err=None)
    assert len(c.depths) == 50
