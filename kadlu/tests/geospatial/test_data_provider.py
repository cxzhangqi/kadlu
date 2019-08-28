""" Unit tests for the the 'geospatial.data_provider' module in the 'kadlu' package

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
from kadlu.geospatial.data_provider import DataProvider

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")

def test_load_bathymetry_from_a_single_chs_file():
    provider = DataProvider(bathy="CHS", south=43, west=-60, north=44, east=-59)
    bathy_data = provider.bathy_data
    bathy = bathy_data[0]
    lats = bathy_data[1]
    lons = bathy_data[2]
    assert np.ma.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.ma.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]

def test_interpolate_chs_bathymetry():
    provider = DataProvider(bathy="CHS", south=43, west=-60, north=44, east=-59)
    N = 10
    x = y = np.arange(N) + 1
    provider.bathy(x,y)
    provider.bathy_gradient(x,y,axis='x')
    provider.bathy_gradient(x,y,axis='y')

def test_interpolate_uniform_bathymetry():
    provider = DataProvider(bathy=-2000)
    N = 10
    x = y = np.arange(N) + 1
    b = provider.bathy(x,y)
    assert np.all(b == -2000)
    bx = provider.bathy_gradient(x,y,axis='x')
    assert np.all(bx == 0)
    by = provider.bathy_gradient(x,y,axis='y')
    assert np.all(by == 0)

def test_query_for_bathymetry_uniform_data():
    provider = DataProvider(bathy=-2000)
    b = provider.bathy()
    assert b == -2000

def test_query_for_bathymetry_grid_data():
    lat = np.array([44, 45, 46, 47, 48])
    lon = np.array([60, 61, 62, 63])
    bathy = np.random.rand(len(lat),len(lon))
    provider = DataProvider(bathy=(bathy,lat,lon))
    b = provider.bathy()
    assert isinstance(b, tuple)
    assert np.all(b[0] == bathy)
    assert np.all(b[1] == lat)
    assert np.all(b[2] == lon)

def test_load_uniform_temperature():
    provider = DataProvider(temp=4)
    _ = provider.temp_data

def test_interpolate_uniform_temperature():
    provider = DataProvider(temp=4)
    N = 10
    x = y = z = np.arange(N) + 1
    temp = provider.temp(x, y, z)
    assert temp.shape[0] == N
    assert np.all(temp == 4)

def test_interpolate_uniform_temperature_on_grid():
    provider = DataProvider(temp=8)
    x = np.arange(10) + 1
    y = np.arange(11) + 1
    z = np.arange(12) + 1
    # planar coordinates
    temp = provider.temp(x=x, y=y, z=z, grid=True)
    assert temp.shape[0] == 10
    assert temp.shape[1] == 11
    assert temp.shape[2] == 12
    assert np.all(temp == 8)
    # spherical coordinates
    lat = np.arange(10)
    lon = np.arange(11)
    temp = provider.temp(x=lon, y=lat, z=z, grid=True, geometry='spherical')
    assert temp.shape[0] == 10
    assert temp.shape[1] == 11
    assert temp.shape[2] == 12
    assert np.all(temp == 8)

def test_interpolate_uniform_salinity():
    provider = DataProvider(salinity=35)
    N = 10
    x = y = z = np.arange(N) + 1
    x = x * 10000
    y = y * 10000
    z = z * 100
    salinity = provider.salinity(x, y, z)
    assert salinity.shape[0] == N
    assert np.all(salinity == 35)

def test_interpolate_uniform_wave():
    provider = DataProvider(wave=1.5)
    N = 10
    x = y = np.arange(N) + 1
    wave = provider.wave(x, y)
    assert wave.shape[0] == N
    assert np.all(wave == 1.5)
