""" Unit tests for the the 'sound.geophony' module in the 'kadlu' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import pytest
import math
import numpy as np
from kadlu.sound.geophony import Geophony
from kadlu.sound.sound_propagation import TLCalculator, Seafloor
from kadlu.geospatial.ocean import Ocean
from kadlu.utils import R1_IUGG, deg2rad


def test_initialize_geophony():
    s = Seafloor()
    o = Ocean()
    tl = TLCalculator(ocean=o, seafloor=s)
    geo = Geophony(tl_calculator=tl, depth=[-100, -200, -300])

def test_create_grid():
    s = Seafloor()
    o = Ocean()
    tl = TLCalculator(ocean=o, seafloor=s, radial_range=10E3)
    geo = Geophony(tl_calculator=tl, depth=[-100], xy_res=1000) # 1km xy grid
    lats, lons, x, y = geo._create_grid(south=0, north=2, west=60, east=62)
    num_bins = int((2 * deg2rad * R1_IUGG + 20E3)/1000)
    num_bins += num_bins%2
    num_bins += 1
    assert num_bins*num_bins == lats.shape[0]
    assert num_bins*num_bins == lons.shape[0]

def test_model_geophony():
    s = Seafloor(thickness=2000)
    o = Ocean(bathy=-10000, wave=1.0)
    tl = TLCalculator(ocean=o, seafloor=s, sound_speed=1480, radial_bin=100, radial_range=50e3, angular_bin=45, vertical_bin=100)
    geo = Geophony(tl_calculator=tl, depth=[-100, -2000])
    spl, x, y = geo.model(frequency=1000, south=44, north=46, west=-60, east=-58)
    assert x.shape[0] == 5
    assert y.shape[0] == 5
    assert spl.shape[0] == 5
    assert spl.shape[1] == 5
    assert spl.shape[2] == 2
    assert np.all(np.diff(x) == np.sqrt(2) * 50e3)
    assert np.all(np.diff(y) == np.sqrt(2) * 50e3)
