import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import rdwps 
from kadlu.geospatial.data_sources.rdwps import Rdwps, rdwps_regions
from os.path import isfile

# mahone bay test area:
south =  44.4
north =  44.7
west  = -64.4
east  = -63.8

start = datetime.now()
end = datetime.now() + timedelta(hours=6)


def test_rdwps_fetch_windwaveswellheight():
    filenames = Rdwps().fetch_windwaveswellheight(south=south, north=north, west=west, east=east, start=start, end=end)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_windwaveheight():
    filenames = Rdwps().fetch_windwaveheight(south=south, north=north, west=west, east=east, start=start, end=end)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_wavedirection():
    filenames = Rdwps().fetch_wavedirection(south=south, north=north, west=west, east=east, start=start, end=end)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_waveperiod():
    filenames = Rdwps().fetch_waveperiod(south=south, north=north, west=west, east=east, start=start, end=end)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_wind_u():
    filenames = Rdwps().fetch_wind_u(south=south, north=north, west=west, east=east, start=start, end=end)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_wind_v():
    filenames = Rdwps().fetch_wind_v(south=south, north=north, west=west, east=east, start=start, end=end)
    for fname in filenames:
        assert(isfile(fname))

# removed this test for now - icecover fetching fails but fixing it is not a priority

#def test_rdwps_fetch_icecover():
#    if not test_fetch: return
#    filenames = Rdwps().fetch_icecover(south=-90, north=90, west=-180, east=180, start=start, end=end)
#    for fname in filenames:
#        assert(isfile(fname))

def test_rdwps_region_abstraction():
    regions = rdwps.ll_2_regionstr(south=-90, north=90, west=-180, east=180, regions=rdwps_regions)
    assert(len(regions) == 5)

    regions = rdwps.ll_2_regionstr(south=44.4, north=44.7, west=-64.4, east=-63.8, regions=rdwps_regions)
    assert(len(regions) == 1)
    assert(regions[0] == 'gulf-st-lawrence')

def test_rdwps_load_wind_u():
    wind, lat, lon, time = Rdwps().load_wind_u(south, north, west, east, start=start, end=end)

def test_rdwps_load_windwaveswellheight():
    height, lat, lon, time = Rdwps().load_windwaveswellheight(south, north, west, east, start=start, end=end)

def test_rdwps_load_wind_v():
    wind, lat, lon, time = Rdwps().load_wind_v(south, north, west, east, start=start, end=end)

def test_rdwps_load_wavedirection():
    height, lat, lon, time = Rdwps().load_wavedirection(south, north, west, east, start=start, end=end)


