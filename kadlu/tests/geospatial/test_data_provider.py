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
    folder = os.path.join(path_to_assets, "tif")
    provider = DataProvider(storage_location=folder, bathy_source="CHS", south=43, west=-60, north=44, east=-59)
    bathy_data = provider.bathy_data
    bathy = bathy_data[0]
    lats = bathy_data[1]
    lons = bathy_data[2]
    assert np.ma.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.ma.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]

def test_interpolate_bathymetry():
    folder = os.path.join(path_to_assets, "tif")
    provider = DataProvider(storage_location=folder, bathy_source="CHS", south=43, west=-60, north=44, east=-59)
    N = 10
    x = y = np.arange(N) + 1
    provider.bathy(x,y)
    provider.bathy_gradient(x,y,axis='x')
    provider.bathy_gradient(x,y,axis='y')
    
def test_load_dummy_temperature():
    provider = DataProvider(storage_location='')
    temp_data = provider.temp_data

def test_interpolate_dummy_temperature():
    provider = DataProvider(storage_location='')
    N = 10
    x = y = z = np.arange(N) + 1
    temp = provider.temp(x, y, z)
    assert temp.shape[0] == N

def test_interpolate_dummy_temperature_on_grid():
    provider = DataProvider(storage_location='')
    x = np.arange(10) + 1
    y = np.arange(11) + 1
    z = np.arange(12) + 1
    # planar coordinates
    temp = provider.temp(x=x, y=y, z=z, grid=True)
    assert temp.shape[0] == 10
    assert temp.shape[1] == 11
    assert temp.shape[2] == 12
    # spherical coordinates
    lat = np.arange(10)
    lon = np.arange(11)
    temp = provider.temp(x=lon, y=lat, z=z, grid=True, geometry='spherical')
    assert temp.shape[0] == 10
    assert temp.shape[1] == 11
    assert temp.shape[2] == 12

def test_interpolate_dummy_salinity():
    provider = DataProvider(storage_location='')
    N = 10
    x = y = z = np.arange(N) + 1
    x = x * 10000
    y = y * 10000
    z = z * 100
    salinity = provider.salinity(x, y, z)
    assert salinity.shape[0] == N

def test_interpolate_dummy_wave():
    provider = DataProvider(storage_location='')
    N = 10
    x = y = np.arange(N) + 1
    wave = provider.wave(x, y)
    assert wave.shape[0] == N