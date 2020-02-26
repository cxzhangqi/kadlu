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
from datetime import datetime
#from kadlu.utils import LatLon
from kadlu.geospatial.ocean import Ocean
from kadlu.geospatial.data_sources.source_map import default_val

#path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")

bounds = dict(
        start=datetime(2015, 1, 9), end=datetime(2015, 1, 9, 3),
        south=44,                   west=-64.5, 
        north=46,                   east=-62.5, 
        top=0,                      bottom=5000
    )
test_lat, test_lon, test_depth = bounds['south'], bounds['west'], bounds['top']

def test_null_ocean():
    """ Test that ocean is initialized with all variables set to 
        null (0) when default=False"""
    #o = Ocean(default=False, cache=False)
    
    # changed to make ocean null by default
    o = Ocean(**bounds)
    
    assert o.bathy(test_lat, test_lon) == 0
    assert o.temp(test_lat, test_lon, test_depth) == 0
    assert o.salinity(test_lat, test_lon, test_depth) == 0
    assert o.wavedir(test_lat, test_lon) == 0
    assert o.waveheight(test_lat, test_lon) == 0
    assert o.waveperiod(test_lat, test_lon) == 0
    assert o.windspeed(test_lat, test_lon) == 0
    #assert o.origin == LatLon(0,0)
    #assert o.SW == LatLon(-90,-180)
    #assert o.NE == LatLon(90,180)

def test_uniform_bathy():
    """ Test that ocean can be initialized with uniform bathymetry"""
    #o = Ocean(default=False, cache=False, load_bathymetry=-500.5)
    o = Ocean(load_bathymetry=-500.5, **bounds)

    assert o.bathy(test_lat, test_lon) == -500.5
    assert o.temp(test_lat, test_lon, test_depth) == 0

def test_interp_uniform_temp():
    """ Test that we can interpolate a uniform ocean temperature 
        on any set of coordinates"""
    #o = Ocean(default=False, cache=False, load_temp=16.1)
    o = Ocean(load_temp=16.1, **bounds)
    assert o.temp(lat=41.2, lon=-66.0, depth=-33.0) == 16.1
    #assert o.temp_xy(x=1, y=2.2, z=-3.0) == 16.1
    #assert np.all(o.temp_xy(x=[5,20], y=[0,10], z=[-300,-400]) == [16.1, 16.1])

def test_uniform_bathy_deriv():
    """ Test that uniform bathy has derivative zero"""
    #o = Ocean(default=False, cache=False, load_bathymetry=-500.5)
    o = Ocean(load_bathymetry=-500.5, **bounds)
    assert o.bathy_deriv(lat=1,lon=17,axis='lon') == 0

def test_chs_bathy():
    """ Test that ocean can be initialized with bathymetry data 
        from a CHS file with automatic fetching enabled"""
    #o = Ocean(default=False, cache=False, fetch=True,
    #    load_bathymetry='chs', south=43.1, west=-59.8, 
    #    north=43.8, east=-59.2)
    bound_args = bounds.copy()
    bound_args['south'], bound_args['west'], bound_args['north'], bound_args['east'] = 43.1, -59.8, 43.8, -59.2
    o = Ocean(fetch=8, load_bathymetry='chs', **bound_args)
    #(bathy,lats,lons) = o.bathy()
    test_lat = [43.4, 43.5]
    test_lon = [-59.6, -59.5]
    bathy = o.bathy(test_lat, test_lon)
    assert len(bathy) > 0 #check that some data was retrieved
    assert  43.1 <= np.min(test_lat) and np.max(test_lat) <=  43.8 #check that lats are within limits
    assert -59.8 <= np.min(test_lon) and np.max(test_lon) <= -59.2 #check that lons are within limits
    # check that all nodes have meaningful bathymetry values
    assert np.all(bathy < 10000)
    assert np.all(bathy > -15000)
    # also check boundaries and origin
    #assert o.origin == LatLon(43.45,-59.5)
    #assert o.SW == LatLon(43.1,-59.8)
    #assert o.NE == LatLon(43.8,-59.2)

def test_interp_chs_bathy():
    """ Test that we can interpolate bathymetry data 
        obtained from a CHS file"""
    #o = Ocean(default=False, cache=False, fetch=True,
    #    load_bathymetry='chs', south=43.1, west=-59.8, 
    #    north=43.8, east=-59.2)
    #b = o.bathy(x=1, y=2)
    o = Ocean(load_bathymetry='chs', 
            south=43.1, west=-59.8, north=43.8, east=-59.2, 
            top=0, bottom=0, start=default_val['start'], end=default_val['end'])
    b = o.bathy(lat=1, lon=2)

    assert isinstance(b, float)
    b = o.bathy(lat=[43.2,43.7], lon=[-59.3, -59.4])
    assert len(b) == 2

def test_hycom_temp_time_interval():
    """ Test that ocean can be initialized with temperature data 
        from HYCOM with automatic fetching enabled and using 
        start/end args.
    """
    #o = Ocean(cache=False, fetch=True,
    #    load_temp='hycom', 
    #    south=43.1, west=-59.8, 
    #    north=43.8, east=-59.2,
    #    top=-100, bottom=3000,
    #    start=datetime(2015,1,1),
    #    end=datetime(2015,1,2)
    #    )
    o = Ocean(#fetch=True, #cache=False,
        load_temp='hycom', 
        south=43.1, west=-59.8, 
        north=43.8, east=-59.2,
        top=-100, bottom=3000,
        start=datetime(2015,1,1),
        end=datetime(2015,1,2)
        )
    lats = [43.4, 43.5]
    lons = [-59.6, -59.5]
    depths = [200, 300]
    #(temp,lats,lons,depths) = o.temp()
    temp = o.temp(lats, lons, depths)
    assert len(temp) > 0 #check that some data was retrieved
    assert  43.1 <= np.min(lats) and np.max(lats) <=  43.8 #check that lats are within limits
    assert -59.8 <= np.min(lons) and np.max(lons) <= -59.2 #check that lons are within limits
    #assert -3000 <= np.min(depths) and np.max(depths) <= 100 #check that depths are within limits
    assert 3000 >= np.max(depths) and np.min(depths) >= -100 #check that depths are within limits
    assert len(depths) > 1 #check that more than one depth was fetched

def test_hycom_temp_nearest_time():
    """ Test that ocean can be initialized with temperature data 
        from HYCOM with automatic fetching enabled and using the 
        time arg"""

    # low-priority feature, passing this test until we can justify implementation
    pass
    return

    o = Ocean(default=False, fetch=True,
        load_temp='hycom', 
        south=43.1, west=-59.8, 
        north=43.8, east=-59.2,
        top=-100, bottom=3000,
        time=datetime(2015,1,1))

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
    #o = Ocean(default=False, cache=False, fetch=True,
    #    load_bathymetry=(bathy,lats,lons))
    
    # note that fetching does nothing when supplying raw array data
    o = Ocean(load_bathymetry=(bathy, lats, lons), fetch=True, **bounds)

    #(b,la,lo) = o.bathy()
    la = lats #[44.55, 44.65]
    lo = lons #[-60, -59.75]
    b = o.bathy(lat=lats, lon=lons)
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
