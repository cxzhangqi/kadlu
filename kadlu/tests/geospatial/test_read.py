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
from kadlu.geospatial.read import read_netcdf_2d, read_matlab_2d

path_to_assets = os.path.join(os.path.dirname(__file__),"../assets")

def test_can_read_bathymetry_from_netcdf_file():
    path = path_to_assets + '/bornholm.nc'
    bathy, lat, lon = read_netcdf_2d(path=path, val_name='bathy', lat_name='lat', lon_name='lon')
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_can_read_bathymetry_from_matlab_file():
    path = path_to_assets + '/bornholm.mat'
    bathy, lat, lon = read_matlab_2d(path=path, val_name='bathy', lat_name='lat', lon_name='lon')
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]
