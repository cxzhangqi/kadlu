""" Unit tests for the the 'geospatial.bathy_chs' module in the 'kadlu' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import pytest
import numpy as np
import os
import kadlu.geospatial.data_sources.chs as chs
from kadlu.utils import LatLon

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "assets")

def test_fetch_returns_two_files():
    folder = os.path.join(path_to_assets, "tif")
    paths = chs.fetch(storage_location=folder, south=42, west=-61, north=44.6, east=-59.1)
    assert len(paths) == 2
    assert os.path.basename(paths[0]) == "CA2_4300N06000W.tif"
    assert os.path.basename(paths[1]) == "CA2_4400N06000W.tif"

def test_filename():
    fname = chs.filename(south=50, west=-62)
    assert fname == "CA2_5000N06200W.tif"

def test_parse_sw_corner():
    s,w = chs.parse_sw_corner("CA2_5000N06200W.tif")
    assert s == 50
    assert w == -62

def test_latlon():
    lats, lons = chs.latlon(path="CA2_5000N06200W.tif")
    assert lats.shape[0] == 1001
    assert lons.shape[0] == 1001
    assert lats[0] == 50
    assert lats[1] == 50.001
    assert lons[0] == -62
    assert lons[1] == -61.999
    lats, lons = chs.latlon(path="CA2_5000N06200W.tif", num_lat=501)
    assert lats.shape[0] == 501
    assert lats[0] == 50
    assert lats[1] == 50.001

def test_load_single_chs_file():
    folder = os.path.join(path_to_assets, "tif")
    bathy, lats, lons = chs.load(storage_location=folder, south=43, west=-60, north=44, east=-59)
    assert np.ma.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.ma.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]

def test_load_multiple_chs_files():
    folder = os.path.join(path_to_assets, "tif")
    bathy1, lats1, lons1 = chs.load(storage_location=folder, south=43, west=-60, north=44, east=-59)
    bathy2, lats2, lons2 = chs.load(storage_location=folder, south=44, west=-60, north=45, east=-59)
    bathy12, lats12, lons12 = chs.load(storage_location=folder, south=43, west=-60, north=45, east=-59)
    assert bathy12.shape[0] == bathy1.shape[0] + bathy2.shape[0]
    assert np.min(lats12) == np.min(lats1)
    assert np.max(lats12) == np.max(lats2)
    assert min(np.min(lons1),np.min(lons2)) == np.min(lons12)
    assert max(np.max(lons1),np.max(lons2)) == np.max(lons12)

def test_load_partial_chs_file():
    folder = os.path.join(path_to_assets, "tif")
    bathy, lats, lons = chs.load(storage_location=folder, south=43, west=-60, north=43.5, east=-59.5)
    assert np.ma.min(bathy) == pytest.approx(-2695.4, abs=0.1)
    assert np.ma.max(bathy) == pytest.approx(-493.8, abs=0.1)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]
