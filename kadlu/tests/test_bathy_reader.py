""" Unit tests for the the 'bathy_reader' module in the 'kadlu' package

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
from kadlu.bathy_reader import BathyReader, LatLon

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_read_bathymetry_from_netcdf_file():
    path = path_to_assets + '/bornholm.nc'
    reader = BathyReader(input=path, bathy_name='bathy')
    lat, lon, bathy = reader.read()
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_can_read_bathymetry_from_matlab_file():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(input=path, bathy_name='bathy')
    lat, lon, bathy = reader.read()
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_can_read_bathymetry_from_gridded_matlab_file():
    path = path_to_assets + '/bathy_grid_test.mat'
    reader = BathyReader(input=path, bathy_name='bathy')
    # read without lat-lon constraints
    lat, lon, bathy = reader.read()
    assert np.min(bathy) == -200
    assert np.max(bathy) == 1200
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]
    # read with lat-lon constraints
    lat, lon, bathy = reader.read(latlon_SW=LatLon(45.5,4.5), latlon_NE=LatLon(47.5,6.5))
    assert np.min(bathy) == -133
    assert np.max(bathy) == 800
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_ensure_lat_and_lon_are_strictly_increasing():
    path = path_to_assets + '/bathy_grid_test.mat'
    reader = BathyReader(input=path, bathy_name='bathy')
    lat, lon, _ = reader.read()
    assert np.all(np.diff(lat) > 0)
    assert np.all(np.diff(lon) > 0)    

def test_ensure_lat_and_lon_are_strictly_increasing_for_dbarclays_data():
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    reader = BathyReader(input=path, lat_name='latgrat', lon_name='longrat', bathy_name='mat', lon_axis=0)
    lat, lon, _ = reader.read()
    assert np.all(np.diff(lat) > 0)
    assert np.all(np.diff(lon) > 0)    

def test_can_read_bathymetry_from_single_geotiff_chs_file():
    path = path_to_assets + '/CA2_4300N06000W.tif'
    reader = BathyReader(input=path, format='GEOTIFF_CHS')
    lat, lon, bathy = reader.read()
    assert np.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[0] == lon.shape[0]
