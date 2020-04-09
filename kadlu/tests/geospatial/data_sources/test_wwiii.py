import pytest
import logging
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import wwiii
from kadlu.geospatial.data_sources.wwiii import Wwiii, Boundary, wwiii_regions, wwiii_global


start   = datetime(2014, 2, 3, 0, 0, 0, 0)
end     = datetime(2014, 2, 3, 3, 0, 0, 0)

# gulf st lawrence
south, west = 47,   -64
north, east = 48,   -63


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

"""
def test_wwiii_fetch_windwaveheight():
    if not Wwiii().fetch_windwaveheight(south=south, north=north, west=west, east=east, start=start, end=end):
        print('wwiii query was fetched already, skipping...')
    return 

def test_wwiii_fetch_wavedirection():
    if not Wwiii().fetch_wavedirection(south=south, north=north, west=west, east=east, start=start, end=end):
        print('wwiii query was fetched already, skipping...')
    return 
def test_wwiii_fetch_waveperiod():
    if not Wwiii().fetch_waveperiod(south=south, north=north, west=west, east=east, start=start, end=end):
        print('wwiii query was fetched already, skipping...')
    return 

def test_wwiii_fetch_windwaveheight():
    if not Wwiii().fetch_windwaveheight(south=south, north=north, west=west, east=east, start=start, end=end):
        print('wwiii query was fetched already, skipping...')
    return 
"""

def test_wwiii_load_windwaveheight():
    wave, lat, lon, time = Wwiii().load_windwaveheight(south=south, west=west, north=north, east=east, start=start, end=end)


"""
def test_wwiii_fetch_wind():
    if not Wwiii().fetch_wind_uv(south=south, north=north, west=west, east=east, start=start, end=end):
        print('wwiii query was fetched already, skipping...')
    return 
"""

def test_wwiii_load_wind():
    try:
        wave, lat, lon, time = Wwiii().load_wind_uv(south=south, west=west, north=north, east=east, start=start, end=end)
    except AssertionError as err:
        logging.info(f'CAUGHT EXCEPTION: {str(err)}')
    except Exception as err:
        raise err

"""
def test_wwiii_fetch_waveperiod():
    if not Wwiii().fetch_waveperiod(south=south, north=north, west=west, east=east, start=start, end=end):
        print('wwiii query was fetched already, skipping...')
"""

def test_wwiii_load_waveperiod():
    wave, lat, lon, time = Wwiii().load_waveperiod(south=south, west=west, north=north, east=east, start=start, end=end)
  

