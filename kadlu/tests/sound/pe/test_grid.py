""" Unit tests for the the 'sound.pe.grid' module in the 'kadlu' package

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
from kadlu.sound.pe.grid import Grid

def test_initialize_grid(grid):
    _ = Grid(100, 1000, 10*np.pi/180, 2*np.pi, 100, 500)

def test_grid_axes_are_as_expected(grid):
    g = Grid(100, 1000, 45*np.pi/180, 2*np.pi, 100, 500)
    # radial
    assert g.dr == 100
    assert g.Nr == 11
    answ = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    assert np.all(g.r == answ)
    # azimuthal
    assert g.dq == 45*np.pi/180
    assert g.Nq == 8
    answ = np.array([0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi, -1.0*np.pi, -0.75*np.pi, -0.50*np.pi, -0.25*np.pi])
    assert np.all(np.abs(g.q - answ) < 1E-9)
    # vertical
    assert g.dz == 100
    assert g.Nz == 10
    answ = np.array([0., 100., 200., 300., 400., -500., -400., -300., -200., -100.])
    assert np.all(g.z == answ)
    # below/above indices
    answ = np.array([0, 9, 8, 7, 6, 5])
    assert np.all(g._indices_below == answ)
    answ = np.array([1, 2, 3, 4])
    assert np.all(g._indices_above == answ)
    answ = np.array([9, 8, 7, 6])
    assert np.all(g._indices_mirror == answ)
    answ = np.meshgrid(np.arange(8), np.array([0, 9, 8, 7, 6, 5]))
    # mirror method
    a = np.array([6.0, np.nan, np.nan, np.nan, np.nan, 1.0, 2.0, 3.0, 4.0, 5.0])
    answ = np.array([6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    a = g.mirror(a)
    assert np.all(a == answ)