
import pytest
import os
import numpy as np
from kadlu.sound.pe.grid import Grid

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

@pytest.fixture
def one():
    return 1

@pytest.fixture
def grid():
    grid = Grid(dr=10, rmax=200, dq=5, qmax=2*np.pi, dz=20, zmax=400)
    return grid
