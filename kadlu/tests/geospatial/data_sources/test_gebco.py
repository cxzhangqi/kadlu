""" Unit tests for the the 'geospatial.data_sources.gebco' module in the 'kadlu' package

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
import kadlu.geospatial.data_sources.gebco as geb

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "assets")

def test_load_bathymetry_for_bornholm():
    path = os.path.join(path_to_assets, "bornholm.nc")
    bathy, lats, lons = geb.load(storage_location=path)
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[1] == lons.shape[0]

def test_load_bathymetry_for_parts_of_bornholm():
    path = os.path.join(path_to_assets, "bornholm.nc")
    bathy1, lats1, lons1 = geb.load(storage_location=path, west=14.5, east=14.9)
    bathy2, lats2, lons2 = geb.load(storage_location=path, west=14.9, east=15.5)
    assert min(np.min(bathy1), np.min(bathy2)) == -100
    assert max(np.min(bathy1), np.min(bathy2)) > -100
    assert max(np.max(bathy1), np.max(bathy2)) == 159
    assert min(np.max(bathy1), np.max(bathy2)) < 159
    assert np.max(lons1) <= 14.9
    assert np.min(lons2) >= 14.9

