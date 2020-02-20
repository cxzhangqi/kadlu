import pytest
import numpy as np
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import era5
from kadlu.geospatial.data_sources.era5 import Era5
from os.path import isfile

# gulf st lawrence
kwargs = dict(
        south=40, west=-70, north=59, east=-55,
        start=datetime(2018, 1, 1), end=datetime(2018, 1, 1)
    )

# note: see assertions within fetch function

def test_era5_fetch_windwaveswellheight():
    if not Era5().fetch_windwaveswellheight(**kwargs):
        print('era5 query was fetched already, skipping...')

def test_era5_fetch_wavedirection():
    if not Era5().fetch_wavedirection(**kwargs):
        print('era5 query was fetched already, skipping...')

def test_era5_fetch_waveperiod():
    if not Era5().fetch_waveperiod(**kwargs):
        print('era5 query was fetched already, skipping...')

def test_era5_fetch_wind():
    if not Era5().fetch_wind(**kwargs):
        print('era5 query was fetched already, skipping...')

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

def test_era5_load_wind():
    wind, lat, lon, time = Era5().load_wind(**kwargs)
    assert(len(wind)==len(lat)==len(lon))
    assert(len(lat) > 0)


