""" Unit tests for the the 'sound.pe.propagator' module in the 'kadlu' package

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
from kadlu.sound.pe.propagator import Propagator
from kadlu.sound.sound_propagation import Seafloor
from kadlu.sound.pe.grid import Grid
from kadlu.geospatial.ocean import Ocean


def test_can_initialize_propapagator(grid):
    o = Ocean()
    s = Seafloor()
    s.frequency=10
    s.c0 = 1500
    k0 = 2 * np.pi * 10 / 1500
    _ = Propagator(ocean=o, seafloor=s, c=1480, c0=1500, grid=grid, k0=k0,\
                smooth_len_den=1500/10/4, smooth_len_c=0.001,\
                absorption_layer=0.2, ignore_bathy_gradient=False,\
                bathy_step=1, c_step=1)
