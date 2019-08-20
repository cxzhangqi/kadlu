""" Unit tests for the the 'geospatial.load' module in the 'kadlu' package

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
import kadlu.geospatial.load as geo

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")

def test_load_single_chs_file():
    folder = os.path.join(path_to_assets, "tif")
    bathy, lats, lons = geo.load_bathy(storage_location=folder, south=43, west=-60, north=44, east=-59, source="CHS")
    assert np.ma.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.ma.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]

def test_load_multiple_chs_files():
    folder = os.path.join(path_to_assets, "tif")
    bathy1, lats1, lons1 = geo.load_bathy(storage_location=folder, south=43, west=-60, north=44, east=-59, source="CHS")
    bathy2, lats2, lons2 = geo.load_bathy(storage_location=folder, south=44, west=-60, north=45, east=-59, source="CHS")
    bathy12, lats12, lons12 = geo.load_bathy(storage_location=folder, south=43, west=-60, north=45, east=-59, source="CHS")
    assert bathy12.shape[0] == bathy1.shape[0] + bathy2.shape[0]
    assert np.min(lats12) == np.min(lats1)
    assert np.max(lats12) == np.max(lats2)
    assert min(np.min(lons1),np.min(lons2)) == np.min(lons12)
    assert max(np.max(lons1),np.max(lons2)) == np.max(lons12)

def test_load_partial_chs_file():
    folder = os.path.join(path_to_assets, "tif")
    bathy, lats, lons = geo.load_bathy(storage_location=folder, south=43, west=-60, north=43.5, east=-59.5, source="CHS")
    assert np.ma.min(bathy) == pytest.approx(-2695.4, abs=0.1)
    assert np.ma.max(bathy) == pytest.approx(-493.8, abs=0.1)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]
