""" Unit tests for the the 'geospatial.ocean' module in the 'kadlu' package

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
from kadlu.geospatial.ocean import Ocean

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")

def test_load_bathymetry_from_a_single_chs_file():
    o = Ocean(bathy="CHS")
    storage = os.path.join(path_to_assets, 'tif')
    o.load_bathy(south=43, west=-60, north=44, east=-59, storage=storage)
    bathy_data = o.bathy_data
    bathy = bathy_data[0]
    lats = bathy_data[1]
    lons = bathy_data[2]
    assert np.ma.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.ma.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]

def test_interpolate_chs_bathymetry():
    o = Ocean(bathy="CHS")
    storage = os.path.join(path_to_assets, 'tif')
    o.load_bathy(south=43, west=-60, north=44, east=-59, storage=storage)
    N = 10
    x = y = np.arange(N) + 1
    o.bathy(x,y)
    o.bathy_gradient(x,y,axis='x')
    o.bathy_gradient(x,y,axis='y')

def test_interpolate_uniform_bathymetry():
    o = Ocean(bathy=-2000)
    N = 10
    x = y = np.arange(N) + 1
    b = o.bathy(x,y)
    assert np.all(b == -2000)
    bx = o.bathy_gradient(x,y,axis='x')
    assert np.all(bx == 0)
    by = o.bathy_gradient(x,y,axis='y')
    assert np.all(by == 0)

def test_query_for_bathymetry_uniform_data():
    o = Ocean(bathy=-2000)
    b = o.bathy()
    assert b == -2000

def test_query_for_bathymetry_grid_data():
    lat = np.array([44, 45, 46, 47, 48])
    lon = np.array([60, 61, 62, 63])
    bathy = np.random.rand(len(lat),len(lon))
    o = Ocean(bathy=(bathy,lat,lon))
    b = o.bathy()
    assert isinstance(b, tuple)
    assert np.all(b[0] == bathy)
    assert np.all(b[1] == lat)
    assert np.all(b[2] == lon)

def test_query_for_temperature_uniform_data():
    o = Ocean(temp=4)
    t = o.temp()
    assert t == 4

def test_query_for_temp_grid_data():
    lat = np.array([44, 45, 46, 47, 48])
    lon = np.array([60, 61, 62, 63])
    z = np.array([1000, 2000])
    temp = np.random.rand(len(lat),len(lon),len(z))
    o = Ocean(temp=(temp,lat,lon,z))
    b = o.temp()
    assert isinstance(b, tuple)
    assert np.all(b[0] == temp)
    assert np.all(b[1] == lat)
    assert np.all(b[2] == lon)
    assert np.all(b[3] == z)


def test_interpolate_uniform_temperature():
    o = Ocean(temp=4)
    N = 10
    x = y = z = np.arange(N) + 1
    temp = o.temp(x, y, z)
    assert temp.shape[0] == N
    assert np.all(temp == 4)

def test_interpolate_uniform_temperature_on_grid():
    o = Ocean(temp=8)
    x = np.arange(10) + 1
    y = np.arange(11) + 1
    z = np.arange(12) + 1
    # planar coordinates
    temp = o.temp(x=x, y=y, z=z, grid=True)
    assert temp.shape[0] == 10
    assert temp.shape[1] == 11
    assert temp.shape[2] == 12
    assert np.all(temp == 8)
    # spherical coordinates
    lat = np.arange(10)
    lon = np.arange(11)
    temp = o.temp(x=lon, y=lat, z=z, grid=True, geometry='spherical')
    assert temp.shape[0] == 10
    assert temp.shape[1] == 11
    assert temp.shape[2] == 12
    assert np.all(temp == 8)

def test_interpolate_uniform_salinity():
    o = Ocean(salinity=35)
    N = 10
    x = y = z = np.arange(N) + 1
    x = x * 10000
    y = y * 10000
    z = z * 100
    salinity = o.salinity(x, y, z)
    assert salinity.shape[0] == N
    assert np.all(salinity == 35)

def test_interpolate_uniform_wave():
    o = Ocean(wave=1.5)
    N = 10
    x = y = np.arange(N) + 1
    wave = o.wave(x, y)
    assert wave.shape[0] == N
    assert np.all(wave == 1.5)
