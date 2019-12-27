""" Unit tests for the the 'sound.pe.starter' module in the 'kadlu' package

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
from kadlu.sound.pe.starter import Starter
from kadlu.sound.pe.grid import Grid

def test_can_initialize_starter_with_default_args(grid):
    _ = Starter(k0=1, grid=grid)

def test_initializing_starter_with_invalid_method_gives_error(grid):
    with pytest.raises(ValueError):
        _ = Starter(k0=1, grid=grid, method='abc')

def test_can_initialize_starter_with_user_args(grid):
    _ = Starter(k0=1, grid=grid, method='GAUSSIAN', aperture=85)

def test_thomson_gives_expected_result():
    grid = Grid(100., 1000., 10.*np.pi/180., 2.*np.pi, 100., 500.)
    starter = Starter(k0=0.04, grid=grid, aperture=86)
    psi = starter.eval(zs=9)
    psi = np.round(psi[0,:,:], 4)
    answ = np.array([[0.+0.j],[0.0101-0.0101j],[0.0205-0.0205j],\
        [0.0319-0.0319j],[ 0.0451-0.0451j],[0.+0.j],[-0.0451+0.0451j],\
        [-0.0319+0.0319j],[-0.0205+0.0205j],[-0.0101+0.0101j]])
    assert np.all(psi == answ)

def test_thomson_gives_expected_result_with_multiple_depths():
    grid = Grid(100., 1000., 10.*np.pi/180., 2.*np.pi, 100., 500.)
    starter = Starter(k0=0.04, grid=grid, aperture=86)
    psi = starter.eval(zs=[9, 99])
    psi = np.round(psi, 4)
    answ0 = np.array([[0.+0.j],[0.0101-0.0101j],[0.0205-0.0205j],\
        [0.0319-0.0319j],[ 0.0451-0.0451j],[0.+0.j],[-0.0451+0.0451j],\
        [-0.0319+0.0319j],[-0.0205+0.0205j],[-0.0101+0.0101j]])
    assert np.all(psi[0,:,:] == answ0)
    answ1 = np.array([[0.+0.j],[0.1039-0.1039j],[0.1723-0.1723j],\
        [0.1806-0.1806j],[0.1222-0.1222j],[0.+0.j],[-0.1222+0.1222j],\
        [-0.1806+0.1806j],[-0.1723+0.1723j],[-0.1039+0.1039j]])
    assert np.all(psi[1,:,:] == answ1)