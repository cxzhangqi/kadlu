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
from kadlu.sound.geophony import Geophony, source_level_kewley
from kadlu.sound.sound_propagation import TLCalculator, Seafloor
from kadlu.geospatial.ocean import Ocean
from kadlu.utils import R1_IUGG, deg2rad


def test_source_level_kewley():
    sl1 = source_level_kewley(freq=10, wind_speed=0)
    sl2 = source_level_kewley(freq=40, wind_speed=2.57)
    assert sl1 == sl2
    assert sl2 == 40.0
    sl3 = source_level_kewley(freq=40, wind_speed=5.14)
    assert sl3 == 44.0
    sl4 = source_level_kewley(freq=100, wind_speed=5.14)
    assert sl4 == 42.5

def test_initialize_geophony():
    s = Seafloor()
    o = Ocean(bathy=-10000, wave=1.0)
    tl = TLCalculator(ocean=o, seafloor=s)
    geo = Geophony(tl_calculator=tl, south=0, north=2, west=60, east=62, depth=[100, 200, 300])

def test_geophony_has_expected_grid():
    s = Seafloor()
    o = Ocean(bathy=-10000, wave=1.0)
    tl = TLCalculator(ocean=o, seafloor=s, radial_range=10E3)
    geo = Geophony(tl_calculator=tl, south=0, north=2, west=60, east=62, depth=[100], xy_res=1000) # 1km xy grid
    lats = geo.lats
    lons = geo.lons
    x = geo.x
    y = geo.y
    num_bins = int((2 * deg2rad * R1_IUGG + 20E3)/1000)
    num_bins += num_bins%2
    num_bins += 1
    assert num_bins*num_bins == lats.shape[0]
    assert num_bins*num_bins == lons.shape[0]

def test_compute_geophony():
    s = Seafloor()
    o = Ocean(bathy=-10000, wave=1.0)
    tl = TLCalculator(ocean=o, seafloor=s, sound_speed=1480, radial_bin=100, radial_range=50e3, angular_bin=45, vertical_bin=100)
    geo = Geophony(tl_calculator=tl, south=44, north=46, west=-60, east=-58, depth=[100, 2000])
    x = geo.x
    y = geo.y
    spl = geo.compute(frequency=1000)
    assert x.shape[0] == 5
    assert y.shape[0] == 5
    assert spl.shape[0] == 5
    assert spl.shape[1] == 5
    assert spl.shape[2] == 2
    assert np.all(np.diff(x) == np.sqrt(2) * 50e3)
    assert np.all(np.diff(y) == np.sqrt(2) * 50e3)

def test_compute_geophony_in_canyon(bathy_canyon):
    s = Seafloor()
    o = Ocean(bathy=bathy_canyon, wave=1.0)
    south = 43
    north = 46
    west = 60
    east = 62
    z = [100, 1500, 3000]
    tl = TLCalculator(ocean=o, seafloor=s, sound_speed=1480, radial_bin=100, radial_range=50e3, angular_bin=45, vertical_bin=100)
    geo = Geophony(tl_calculator=tl, south=south, north=north, west=west, east=east, depth=z)
    x = geo.x
    y = geo.y
    spl = geo.compute(frequency=10)
    assert spl.shape[0] == x.shape[0]
    assert spl.shape[1] == y.shape[0]
    assert spl.shape[2] == len(z)
    assert np.all(np.diff(x) == np.sqrt(2) * tl.range['r'])
    assert np.all(np.diff(y) == np.sqrt(2) * tl.range['r'])    
    # check that noise is NaN below seafloor and non Nan above
    bathy = (-1.) * np.swapaxes(np.reshape(geo.bathy, newshape=(y.shape[0], x.shape[0])), 0, 1)
    bathy = bathy[:,:,np.newaxis]
    xyz = np.ones(shape=bathy.shape) * z
    idx = np.nonzero(xyz >= bathy)
    assert np.all(np.isnan(spl[idx]))
    idx = np.nonzero(xyz < bathy)
    assert np.all(~np.isnan(spl[idx]))

def test_wind_source_level_per_area():
    s = Seafloor()
    o = Ocean(bathy=-10000, wave=5.14)
    tl = TLCalculator(ocean=o, seafloor=s)
    geo = Geophony(tl_calculator=tl, south=44, north=46, west=-60, east=-58, depth=[100, 200, 300])
    SL_f10 = geo._wind_source_level_per_area(freq=10, x=0, y=0, start=None, end=None)
    SL_f20 = geo._wind_source_level_per_area(freq=10, x=0, y=0, start=None, end=None)
    SL_arr = geo._wind_source_level_per_area(freq=10, x=[0,1,2], y=[0,1,2], start=None, end=None)
    assert SL_f10 == SL_f20
    assert SL_f10 == 39.0
    assert len(SL_arr) == 3
    assert np.all(SL_arr == 39.0)

def test_source_level(grid):
    s = Seafloor()
    o = Ocean(bathy=-10000, wave=5.14)
    tl = TLCalculator(ocean=o, seafloor=s)
    geo = Geophony(tl_calculator=tl, south=44, north=46, west=-60, east=-58, depth=[100, 200, 300])
    SL = geo._source_level(freq=10, grid=grid, start=None, end=None, method='wind')
    assert SL.shape[0] == 1
    assert SL.shape[1] == len(grid.q)
    assert SL.shape[2] == len(grid.r) - 1
