""" Unit tests for the the 'pe_starter' module in the 'kadlu' package

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
from kadlu.transmission_loss.pe_starter import PEStarter

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_initialize_PE_starter_with_default_args(pe_grid):
    _ = PEStarter(ref_wavenumber=1, grid=pe_grid)

def test_initializing_PE_starter_with_invalid_method_gives_error(pe_grid):
    with pytest.raises(ValueError):
        _ = PEStarter(ref_wavenumber=1, grid=pe_grid, method='abc')

def test_can_initialize_PE_starter_with_user_args(pe_grid):
    _ = PEStarter(ref_wavenumber=1, grid=pe_grid, method='GAUSSIAN', aperture=85)
