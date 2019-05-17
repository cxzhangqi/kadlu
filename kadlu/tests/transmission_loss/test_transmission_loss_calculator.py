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

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_initialize_TL_calculator():
    TransmissionLossCalculator(bathymetry=None, sound_speed=None, flat_seafloor_depth=10000,\
        progress_bar=False)

def test_run_TL_calculator():
    calc = TransmissionLossCalculator(bathymetry=None, sound_speed=None, flat_seafloor_depth=10000,\
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

    import matplotlib.pyplot as plt
    calc.plot_vertical(angle=0, show_bathy=True)
    plt.show()