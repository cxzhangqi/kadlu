
import pytest
import os
import numpy as np
from kadlu.transmission_loss.transmission_loss_calculator import PEGrid

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

@pytest.fixture
def one():
    return 1

@pytest.fixture
def pe_grid():
    grid = PEGrid(radial_step=10, radial_range=200,\
            azimuthal_step=5, azimuthal_range=2*np.pi,\
            vertical_step=20, vertical_range=400)

    return grid
