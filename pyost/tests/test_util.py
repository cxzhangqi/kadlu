""" Unit tests for the the 'utility' module in the 'pyost' package

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
from pyost.util import LLtoXY, XYtoLL

def test_can_convert_single_point_from_ll_to_xy():
    lat_ref = 45
    lon_ref = 10
    lats = 45
    lons = 20
    x, y = LLtoXY(lats, lons, lat_ref, lon_ref)    
    assert x == pytest.approx(788E3, 1E3)
    assert y == pytest.approx(0, 1E3)
    lats = 46
    lons = 10
    x, y = LLtoXY(lats, lons, lat_ref, lon_ref)    
    assert x == pytest.approx(0, 1E3)
    assert y == pytest.approx(111E3, 1E3)

def test_can_convert_several_points_from_ll_to_xy():
    lat_ref = 45
    lon_ref = 10
    lats = [45, 45]
    lons = [20, 30]
    x, y = LLtoXY(lats, lons, lat_ref, lon_ref)    
    assert x[0] == pytest.approx(788E3, 1E3)
    assert y[0] == pytest.approx(0, 1E3)
    assert x[1] == 2*x[0]
    assert y[1] == pytest.approx(0, 1E3)

def test_can_convert_single_point_from_xy_to_ll():
    lat_ref = 45
    lon_ref = 10
    xs = 788E3
    ys = 0
    lat, lon = XYtoLL(xs, ys, lat_ref, lon_ref)    
    assert lon == pytest.approx(lon_ref+10, 0.1)
    assert lat == pytest.approx(lat_ref, 0.1)
    xs = 0
    ys = 111E3
    lat, lon = XYtoLL(xs, ys, lat_ref, lon_ref)    
    assert lon == pytest.approx(lon_ref, 0.1)
    assert lat == pytest.approx(lat_ref+1, 0.1)
