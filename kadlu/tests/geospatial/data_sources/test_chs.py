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
import os
import kadlu.geospatial.data_sources.chs as chs
from kadlu.geospatial.bathy_reader import LatLon

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


    