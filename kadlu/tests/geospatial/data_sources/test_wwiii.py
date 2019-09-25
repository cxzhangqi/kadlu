import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import wwiii
from kadlu.geospatial.data_sources.wwiii import Wwiii
from kadlu.geospatial.data_sources import fetch_util


test_fetch = True

start = datetime(2017, 2, 3, 0, 0, 0, 0)
end = datetime(2017, 2, 3, 0, 0, 0, 0)

# gulf st lawrence
south =  46
north =  52
west  = -70
east  = -56

def test_wwiii_ll2regionstr():
    """
    south=-80
    north=84
    west=-180
    east=180
    regions = wwiii.ll_2_regionstr(south, north, west, east)
    assert(len(regions) == 1)
    assert(regions[0] == 'glo_30m')
    """
    # gulf st lawrence
    south =  46
    north =  52
    west  = -70
    east  = -56
    regions = wwiii.ll_2_regionstr(south, north, west, east)
    print(f"Regions for gulf st lawrence coords : {regions}")


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
   

