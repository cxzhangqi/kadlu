import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import fetch_util
from kadlu.geospatial.data_sources import wwiii
from kadlu.geospatial.data_sources.wwiii import Wwiii, Boundary, wwiii_regions, wwiii_global


test_fetch = True

start = datetime(2017, 2, 3, 0, 0, 0, 0)
end = datetime(2017, 2, 3, 0, 0, 0, 0)

# gulf st lawrence
south =  46
north =  52
west  = -70
east  = -56

def test_wwiii_ll2regionstr():
    
    # gulf st lawrence
    south =  46
    north =  52
    west  = -70
    east  = -56
    regions = wwiii.ll_2_regionstr(south, north, west, east, wwiii_regions, wwiii_global)
    assert(len(regions) == 1)
    assert(regions[0] == 'at_4m')

    # bering sea
    # test area intersection across antimeridian 
    south, north = 46, 67
    west, east = 158, -156
    east=-156
    regions = wwiii.ll_2_regionstr(south, north, west, east, wwiii_regions, wwiii_global)
    assert(len(regions) == 3)
    assert('ak_4m' in regions)
    assert('ao_30m' in regions)
    assert('wc_4m' in regions)

    # global 
    globe = wwiii.ll_2_regionstr(-90, 90, -180, 180, wwiii_regions, wwiii_global)
    assert(len(globe) == 5)


def test_wwiii_fetch_windwaveheight():
    if not test_fetch: return
    filenames = Wwiii().fetch_windwaveheight(south, north, west, east, start, end)


def test_wwiii_load_windwaveheight():
    wave, lat, lon, time = Wwiii().load_windwaveheight(south=43, west=-60, north=44, east=-59, start=start, end=end)
   

def test_wwiii_fetch_waveperiod():
    if not test_fetch: return
    filenames = Wwiii().fetch_waveperiod(south, north, west, east, start, end)


def test_wwiii_load_waveperiod():
    wave, lat, lon, time = Wwiii().load_waveperiod(south=43, west=-60, north=44, east=-59, start=start, end=end)
   

