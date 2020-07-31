import pytest
import kadlu
import numpy as np
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources.era5 import Era5
from os.path import isfile
from kadlu.geospatial.data_sources.fetch_handler import fetch_handler

# gulf st lawrence
kwargs = dict(
        south=40, west=-65, north=41, east=-64,
        start=datetime(2018, 2, 1), end=datetime(2018, 2, 3)
    )

# note: see assertions within fetch function

def test_era5_fetch_windwaveswellheight():
    #if not Era5().fetch_windwaveswellheight(**kwargs):
    #    print('era5 query was fetched already, skipping...')
    fetch_handler('waveheight', 'era5', **kwargs)

def test_era5_fetch_wavedirection():
    #if not Era5().fetch_wavedirection(**kwargs):
    #    print('era5 query was fetched already, skipping...')
    fetch_handler('wavedir', 'era5', **kwargs)

def test_era5_fetch_waveperiod():
    #if not Era5().fetch_waveperiod(**kwargs):
    #    print('era5 query was fetched already, skipping...')
    fetch_handler('waveperiod', 'era5', **kwargs)

def test_era5_fetch_wind():
    #if not Era5().fetch_wind(**kwargs):
    #    print('era5 query was fetched already, skipping...')
    fetch_handler('wind_u', 'era5', **kwargs)
    fetch_handler('wind_v', 'era5', **kwargs)
    fetch_handler('wind_uv', 'era5', **kwargs)

def test_era5_load_windwaveswellheight():
    height, lat, lon, time = Era5().load_windwaveswellheight(**kwargs)
    assert(len(height)==len(lat)==len(lon))
    assert(len(lat) > 0)

def test_era5_load_wavedirection():
    wave, lat, lon, time = Era5().load_wavedirection(**kwargs)
    assert(len(wave)==len(lat)==len(lon))
    assert(len(lat) > 0)

def test_era5_load_waveperiod():
    wave, lat, lon, time = Era5().load_waveperiod(**kwargs)
    assert(len(wave)==len(lat)==len(lon))
    assert(len(lat) > 0)

def test_era5_load_wind_uv():
    kwargs = dict(
            south=40, west=-70, north=41, east=-69,
            start=datetime(2018, 2, 1), end=datetime(2018, 2, 2)
        )
    wind, lat, lon, time = Era5().load_wind_uv(**kwargs)
    assert(len(wind)==len(lat)==len(lon))
    assert(len(lat) > 0)


def test_era5_load_wind():
    ns_offset = 1
    ew_offset = 1

    uval,lat,lon,epoch = kadlu.load(source='era5', var='wind_u', 
                        start=datetime(2016, 3, 9) , end=datetime(2016,3,11),
                        south=44.5541333 - ns_offset, west=-64.17682 - ew_offset, 
                        north=44.5541333 + ns_offset, east=-64.17682 + ew_offset, 
                        top=0, bottom=0)

    vval,lat,lon,epoch = kadlu.load(source='era5', var='wind_v', 
                        start=datetime(2016, 3, 9) , end=datetime(2016,3,11),
                        south=44.5541333 - ns_offset, west=-64.17682 - ew_offset, 
                        north=44.5541333 + ns_offset, east=-64.17682 + ew_offset, 
                        top=0, bottom=0)

    uvval,lat,lon,epoch = kadlu.load(source='era5', var='wind_uv', 
                        start=datetime(2016, 3, 9) , end=datetime(2016,3,11),
                        south=44.5541333 - ns_offset, west=-64.17682 - ew_offset, 
                        north=44.5541333 + ns_offset, east=-64.17682 + ew_offset, 
                        top=0, bottom=0)


