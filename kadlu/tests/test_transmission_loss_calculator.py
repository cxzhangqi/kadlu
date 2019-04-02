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
from kadlu.transmission_loss_calculator import TransmissionLossCalculator

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_initialize_TL_calculator():
    c = TransmissionLossCalculator(bathymetry=None, sound_speed=None)

def test_run_TL_calculator():
    c = TransmissionLossCalculator(bathymetry=None, sound_speed=None)
    c.run(frequency=10, source_depth=9905)
