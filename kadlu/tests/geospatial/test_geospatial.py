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
from kadlu.geospatial.geospatial import read_geotiff, read_matlab

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")

def test_can_read_bathymetry_from_matlab_file():
    path = path_to_assets + '/bornholm.mat'
    lat = np.squeeze(read_matlab(path=path, name='lat'))
    lon = np.squeeze(read_matlab(path=path, name='lon'))
    bathy = read_matlab(path=path, name='bathy')
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lat.shape[0]
    assert bathy.shape[1] == lon.shape[0]

def test_can_read_bathymetry_from_geotiff_file():
    path = path_to_assets + '/tif/CA2_4300N06000W.tif'
    bathy = read_geotiff(path=path)
    assert np.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.max(bathy) == pytest.approx(1.645, abs=0.001)
