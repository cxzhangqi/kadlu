""" Unit tests for the the 'geospatial.ocean' module in the 'kadlu' package

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
import datetime
from kadlu.utils import LatLon
from kadlu.geospatial.ocean import Ocean

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")
now = datetime.datetime.now()

def test_null_ocean():
    """ Test that ocean is initialized with all variables set to 
        null (0) when default=False"""
    o = Ocean(default=False, cache=False)
    assert o.bathy() == 0
    assert o.temp() == 0
    assert o.salinity() == 0
    assert o.wavedir() == 0
    assert o.waveheight() == 0
    assert o.waveperiod() == 0
    assert o.windspeed() == 0
    assert o.origin == LatLon(0,0)
    assert o.SW == LatLon(-90,-180)
    assert o.NE == LatLon(90,180)

def test_uniform_bathy():
    """ Test that ocean can be initialized with uniform bathymetry"""
    o = Ocean(default=False, cache=False, load_bathymetry=-500.5)
    assert o.bathy() == -500.5
    assert o.temp() == 0

def test_interp_uniform_temp():
    """ Test that we can interpolate a uniform ocean temperature 
        on any set of coordinates"""
    o = Ocean(default=False, cache=False, load_temp=16.1)
    assert o.temp(x=1, y=2.2, z=-3.0) == 16.1
    assert o.temp(lat=41.2, lon=-66.0, z=-33.0) == 16.1
    assert np.all(o.temp(x=[5,20], y=[0,10], z=[-300,-400]) == [16.1, 16.1])

def test_uniform_bathy_deriv():
    """ Test that uniform bathy has derivative zero"""
    o = Ocean(default=False, cache=False, load_bathymetry=-500.5)
    assert o.bathy_deriv(x=1,y=17,axis='x') == 0

def test_chs_bathy():
    """ Test that ocean can be initialized with bathymetry data 
        from a CHS file with automatic fetching enabled"""
    o = Ocean(default=False, cache=False, fetch=True,
        load_bathymetry='chs', south=43.1, west=-59.8, 
        north=43.8, east=-59.2)
    (bathy,lats,lons) = o.bathy()
    assert len(bathy) > 0 #check that some data was retrieved
    assert  43.1 <= np.min(lats) and np.max(lats) <=  43.8 #check that lats are within limits
    assert -59.8 <= np.min(lons) and np.max(lons) <= -59.2 #check that lons are within limits
    # check that all nodes have meaningful bathymetry values
    assert np.all(bathy < 10000)
    assert np.all(bathy > -15000)
    # also check boundaries and origin
    assert o.origin == LatLon(43.45,-59.5)
    assert o.SW == LatLon(43.1,-59.8)
    assert o.NE == LatLon(43.8,-59.2)

def test_interp_chs_bathy():
    """ Test that we can interpolate bathymetry data 
        obtained from a CHS file"""
    o = Ocean(default=False, cache=False, fetch=True,
        load_bathymetry='chs', south=43.1, west=-59.8, 
        north=43.8, east=-59.2)
    b = o.bathy(x=1, y=2)
    assert isinstance(b, float)
    b = o.bathy(lat=[43.2,43.7], lon=[-59.3, -59.4])
    assert len(b) == 2

def test_hycom_temp():
    """ Test that ocean can be initialized with temperature data 
        from HYCOM with automatic fetching enabled"""
    o = Ocean(default=False, cache=False, fetch=True,
        load_temp='hycom', 
        south=43.1, west=-59.8, 
        north=43.8, east=-59.2,
        top=-100, bottom=3000,
        start=datetime.datetime(2014,1,1,0),
        end=datetime.datetime(2014,1,1,1))
    (temp,lats,lons,depths) = o.temp()
    assert len(temp) > 0 #check that some data was retrieved
    assert  43.1 <= np.min(lats) and np.max(lats) <=  43.8 #check that lats are within limits
    assert -59.8 <= np.min(lons) and np.max(lons) <= -59.2 #check that lons are within limits
    assert -3000 <= np.min(depths) and np.max(depths) <= 100 #check that depths are within limits

def test_array_bathy():
    """ Test that ocean can be initialized with bathymetry data 
        from arrays"""
    bathy = np.array([[-100., -200.],
                      [-100., -200.]])
    lats = np.array([44.5, 44.7])
    lons = np.array([-60.1, -59.5])
    o = Ocean(default=False, cache=False, fetch=True,
        load_bathymetry=(bathy,lats,lons))
    (b,la,lo) = o.bathy()
    assert np.all(b == bathy)
    assert np.all(la == lats)
    assert np.all(lo == lons)
    res = o.bathy(lat=44.5, lon=-60.1)
    assert res == -100
    res = o.bathy(lat=44.5, lon=-59.8)
    assert pytest.approx(res == -150., abs=1e-6)


""" Interactive testing
    south, west = 44, -59
    north, east = 46, -57
    start, end = datetime(2015, 1, 10), datetime(2015, 1, 10, 12)
    top, bottom = 0, 5000
"""
