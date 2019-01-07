""" Unit tests for the the 'bathy_reader' module in the 'pyost' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/pyost
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import pytest
import os
import numpy as np
from pyost.bathy_reader import BathyNetCDFReader, BathyMatReader, LatLon

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_read_bathymetry_from_netcdf_file(one):
    path = path_to_assets + '/bornholm.nc'
    reader = BathyNetCDFReader(path=path, bathy_name='bathy')
    lat, lon, bathy = reader.read()
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_can_read_bathymetry_from_matlab_file(one):
    path = path_to_assets + '/bornholm.mat'
    reader = BathyMatReader(path=path, bathy_name='bathy')
    lat, lon, bathy = reader.read()
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_can_read_bathymetry_from_gridded_matlab_file(one):
    path = path_to_assets + '/bathy_grid_test.mat'
    reader = BathyMatReader(path=path, bathy_name='bathy')
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
