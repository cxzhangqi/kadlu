
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

@pytest.fixture
def bathy_canyon():
    depth = 2e3 #2km
    sigma = 0.5
    def canyon_axis(x):
        y = 45 + (x - 61) 
        return y

    x = np.linspace(60, 62, num=100) #lons
    y = np.linspace(44, 46, num=100) #lats
    yv, xv = np.meshgrid(y, x)
    bathy = -depth * np.exp(-(yv - canyon_axis(xv))**2 / (2 * sigma**2))
    return (bathy, y, x)