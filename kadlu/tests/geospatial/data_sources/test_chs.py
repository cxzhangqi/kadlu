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
from kadlu.geospatial.data_sources.chs import Chs
import numpy as np

source = Chs()

# gulf st lawrence - northumberland strait
south, west = 45.91, -63.81
north, east = 46.12, -62.92

northumberland_strait_deepest_point = -37
# http://fishing-app.gpsnauticalcharts.com/i-boating-fishing-web-app/fishing-marine-charts-navigation.html?title=Northumberland+Strait+boating+app#9/46.0018/-63.1677

def test_fetch_bathy():
    source.fetch_bathymetry(south=south, north=north, west=west, east=east)

def test_load_bathy():
    bathy, lat, lon = source.load_bathymetry(south=south, north=north, west=west, east=east)
    assert (len(bathy) == len(lat) == len(lon))
    assert (len(bathy) > 0)
    assert np.all(np.logical_and(lat >= south, lat <= north))
    assert np.all(np.logical_and(lon >= west, lon <= east))
    assert np.min(bathy) >= northumberland_strait_deepest_point 
    assert np.max(bathy) <= 2 
