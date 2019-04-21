""" Unit tests for the the 'transmission_loss_calculator' module in the 'kadlu' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import pytest
import os
import numpy as np
from kadlu.transmission_loss.transmission_loss_calculator import TransmissionLossCalculator

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_initialize_TL_calculator():
    TransmissionLossCalculator(bathymetry=None, sound_speed=None)

def test_run_TL_calculator():
    calc = TransmissionLossCalculator(bathymetry=None, sound_speed=None, step_size=1000, range=10e3, angular_bin_size=45, vertical_bin_size=1000)
    field = calc.run(frequency=10, source_depth=9900)
    expected = np.array([[-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535],\
        [-164.6453, -170.6553, -176.7944, -172.0352, -182.3293, -176.6379, -176.8878, -183.8019, -177.9633, -181.3535]])
    np.testing.assert_array_almost_equal(field, expected, decimal=3)

#def test2():
#    calc = TransmissionLossCalculator(bathymetry=None, sound_speed=None, step_size=1000, range=10e3, angular_bin_size=45, vertical_bin_size=1000)
#    field = calc.run(frequency=10, source_depth=9900)
