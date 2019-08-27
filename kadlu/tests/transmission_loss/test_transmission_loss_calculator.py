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

""" Unit tests for the transmission loss module within the kadlu library
"""

import pytest
import os
import numpy as np
from kadlu.transmission_loss.transmission_loss_calculator import TransmissionLossCalculator
from kadlu.geospatial.data_provider import DataProvider

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)),"assets")


def test_can_initialize_TL_calculator():
    TransmissionLossCalculator(env_data=None, sound_speed=None, flat_seafloor_depth=10000,\
        progress_bar=False)

def test_run_TL_calculator_with_flat_seafloor():
    calc = TransmissionLossCalculator(env_data=None, sound_speed=None, flat_seafloor_depth=10000,\
        step_size=1000, range=10e3, angular_bin_size=45, vertical_bin_size=1000, verbose=True, progress_bar=False)
    calc.run(frequency=10, source_depth=9900)
    field = calc.TL
    expected = np.array([[-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535]])
    np.testing.assert_array_almost_equal(field, expected, decimal=3)

def test_run_TL_calculator_with_realistic_bathymetry():
    # load bathymetric data
    provider = DataProvider(bathy_source="CHS", south=43, west=-60, north=44, east=-59)
    # depth at center
    seafloor_depth = -provider.bathy(x=0, y=0)
    max_depth = -np.min(provider.bathy_data[0]) 
    # initialize calculator
    calc = TransmissionLossCalculator(env_data=provider, sound_speed=None,\
        step_size=1000, range=10e3, angular_bin_size=45, vertical_bin_size=1000,\
        max_depth=1.2*max_depth, steps_btw_sound_speed_updates=1,\
        verbose=True, progress_bar=False)
    # run
    calc.run(frequency=10, source_depth=0.9*seafloor_depth)
