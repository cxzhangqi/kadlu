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




## old tests below ...

def test_interpolate_chs_bathymetry():
    o = Ocean(bathy="CHS")
    storage = os.path.join(path_to_assets, 'tif')
    o.load_bathy(south=43, west=-60, north=44, east=-59, storage=storage)
    N = 10
    x = y = np.arange(N) + 1
    o.bathy(x,y)
    o.bathy_gradient(x,y,axis='x')
    o.bathy_gradient(x,y,axis='y')

def test_interpolate_uniform_bathymetry():
    o = Ocean(bathy=-2000)
    N = 10
    x = y = np.arange(N) + 1
    b = o.bathy(x,y)
    assert np.all(b == -2000)
    bx = o.bathy_gradient(x,y,axis='x')
    assert np.all(bx == 0)
    by = o.bathy_gradient(x,y,axis='y')
    assert np.all(by == 0)

def test_query_for_bathymetry_uniform_data():
    o = Ocean(bathy=-2000)
    b = o.bathy()
    assert b == -2000

def test_query_for_bathymetry_grid_data():
    lat = np.array([44, 45, 46, 47, 48])
    lon = np.array([60, 61, 62, 63])
    bathy = np.random.rand(len(lat),len(lon))
    o = Ocean(bathy=(bathy,lat,lon))
    b = o.bathy()
    assert isinstance(b, tuple)
    assert np.all(b[0] == bathy)
    assert np.all(b[1] == lat)
    assert np.all(b[2] == lon)

def test_query_for_temperature_uniform_data():
    o = Ocean(temp=4)
    t = o.temp()
    assert t == 4

def test_query_for_temp_grid_data():
    lat = np.array([44, 45, 46, 47, 48])
    lon = np.array([60, 61, 62, 63])
    z = np.array([1000, 2000])
    temp = np.random.rand(len(lat),len(lon),len(z))
    o = Ocean(temp=(temp,lat,lon,z))
    b = o.temp()
    assert isinstance(b, tuple)
    assert np.all(b[0] == temp)
    assert np.all(b[1] == lat)
    assert np.all(b[2] == lon)
    assert np.all(b[3] == z)


def test_interpolate_uniform_temperature():
    o = Ocean(temp=4)
    N = 10
    x = y = z = np.arange(N) + 1
    temp = o.temp(x, y, z)
    assert temp.shape[0] == N
    assert np.all(temp == 4)

def test_interpolate_uniform_temperature_on_grid():
    o = Ocean(temp=8)
    x = np.arange(10) + 1
    y = np.arange(11) + 1
    z = np.arange(12) + 1
    # planar coordinates
    temp = o.temp(x=x, y=y, z=z, grid=True)
    assert temp.shape[0] == 10
    assert temp.shape[1] == 11
    assert temp.shape[2] == 12
    assert np.all(temp == 8)
    # spherical coordinates
    lat = np.arange(10)
    lon = np.arange(11)
    temp = o.temp(x=lon, y=lat, z=z, grid=True, geometry='spherical')
    assert temp.shape[0] == 10
    assert temp.shape[1] == 11
    assert temp.shape[2] == 12
    assert np.all(temp == 8)

def test_interpolate_uniform_salinity():
    o = Ocean(salinity=35)
    N = 10
    x = y = z = np.arange(N) + 1
    x = x * 10000
    y = y * 10000
    z = z * 100
    salinity = o.salinity(x, y, z)
    assert salinity.shape[0] == N
    assert np.all(salinity == 35)

def test_interpolate_uniform_wave():
    o = Ocean(wave=1.5)
    N = 10
    x = y = np.arange(N) + 1
    wave = o.wave(x, y)
    assert wave.shape[0] == N
    assert np.all(wave == 1.5)


""" Interactive testing
    south, west = 44, -59
    north, east = 46, -57
    start, end = datetime(2015, 1, 10), datetime(2015, 1, 10, 12)
    top, bottom = 0, 5000


"""
