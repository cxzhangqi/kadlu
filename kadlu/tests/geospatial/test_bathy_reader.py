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
from kadlu.geospatial.bathy_reader import BathyReader, LatLon, write_bathy

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")

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

def test_can_read_bathymetry_from_single_geotiff_file():
    path = path_to_assets + '/tif/CA2_4300N06000W.tif'
    reader = BathyReader(input=path)
    lat, lon, bathy = reader.read()
    assert np.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[0] == lon.shape[0]

def test_can_read_bathymetry_from_two_geotiff_files():
    path = path_to_assets + '/tif/CA2_4300N06000W.tif'
    r1 = BathyReader(input=path)
    path = path_to_assets + '/tif/CA2_4400N06000W.tif'
    r2 = BathyReader(input=path)
    path = path_to_assets + '/tif/'
    r12 = BathyReader(input=path)
    lat1, lon1, bathy1 = r1.read()
    lat2, lon2, bathy2 = r2.read()
    lat12, lon12, bathy12 = r12.read()
    assert bathy12.shape[0] == bathy1.shape[0] + bathy2.shape[0]
    assert np.min(lat12) == np.min(lat1)
    assert np.max(lat12) == np.max(lat2)
    assert min(np.min(lon1),np.min(lon2)) == np.min(lon12)
    assert max(np.max(lon1),np.max(lon2)) == np.max(lon12)

def test_can_write_bathymetry_to_matlab_file():
    path = path_to_assets + '/tmp.mat'
    # create dummy data and save to file
    lats = [40, 50, 60]
    lons = [-55, -54]
    bathy = [[-440, -460],\
             [-540, -560],\
             [-640, -660]]
    write_bathy(lat=lats, lon=lons, bathy=bathy, destination=path)
    # read back in 
    reader = BathyReader(input=path, lat_name='lat', lon_name='lon', bathy_name='bathy')
    lats1, lons1, bathy1 = reader.read()
    assert lats1.shape[0] == 3
    assert lons1.shape[0] == 2
    assert bathy1.shape[0] == 3
    assert bathy1.shape[1] == 2
    assert np.all(bathy1 == np.array(bathy))
    # clean
    os.remove(path)
