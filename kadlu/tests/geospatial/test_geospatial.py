""" Unit tests for the the 'geospatial.geospatial' module in the 'kadlu' package

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
from kadlu.geospatial.geospatial import read_geotiff, load_data_from_file, write_data_to_file

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")

def test_can_read_bathymetry_from_matlab_file():
    path = path_to_assets + '/bornholm.mat'
    bathy, lat, lon = load_data_from_file(path=path)
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_can_read_bathymetry_from_geotiff_file():
    path = path_to_assets + '/tif/CA2_4300N06000W.tif'
    bathy = read_geotiff(path=path)
    assert np.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.max(bathy) == pytest.approx(1.645, abs=0.001)

def test_can_read_bathymetry_from_netcdf_file():
    path = path_to_assets + '/bornholm.nc'
    bathy, lat, lon = load_data_from_file(path=path)
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_can_read_bathymetry_from_gridded_matlab_file():
    path = path_to_assets + '/bathy_grid_test.mat'
    # read without lat-lon constraints
    bathy, lat, lon = load_data_from_file(path=path)
    assert np.min(bathy) == -200
    assert np.max(bathy) == 1200
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]
    # read with lat-lon constraints
    bathy, lat, lon = load_data_from_file(path=path, south=45.5, west=4.5, north=47.5, east=6.5)
    assert np.min(bathy) == -133
    assert np.max(bathy) == 800
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_ensure_lat_and_lon_are_strictly_increasing():
    path = path_to_assets + '/bathy_grid_test.mat'
    _, lat, lon = load_data_from_file(path=path)
    assert np.all(np.diff(lat) > 0)
    assert np.all(np.diff(lon) > 0)    

def test_ensure_lat_and_lon_are_strictly_increasing_for_dbarclays_data():
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    _, lat, lon = load_data_from_file(path=path, lat_name='latgrat', lon_name='longrat', val_name='mat', lon_axis=0)
    assert np.all(np.diff(lat) > 0)
    assert np.all(np.diff(lon) > 0)    

def test_can_write_bathymetry_to_matlab_file():
    path = path_to_assets + '/tmp.mat'
    # create dummy data and save to file
    lats = [40, 50, 60]
    lons = [-55, -54]
    bathy = [[-440, -460],\
             [-540, -560],\
             [-640, -660]]
    write_data_to_file(lat=lats, lon=lons, val=bathy, destination=path)
    # read back in 
    bathy1, lats1, lons1 = load_data_from_file(path=path, val_name='data')
    assert lats1.shape[0] == 3
    assert lons1.shape[0] == 2
    assert bathy1.shape[0] == 3
    assert bathy1.shape[1] == 2
    assert np.all(bathy1 == np.array(bathy))
    # clean
    os.remove(path)
