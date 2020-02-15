import pytest
import numpy as np
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import era5
from kadlu.geospatial.data_sources.era5 import Era5
from os.path import isfile

# gulf st lawrence
south, west = 40, -70
north, east = 59, -55

start = datetime(2018, 1, 1, 0, 0, 0, 0)
end   = datetime(2018, 1, 1, 0, 0, 0, 0)

# note: see assertions within fetch function

def test_era5_fetch_windwaveswellheight():
    if not Era5().fetch_windwaveswellheight(start=start, end=end):
        print('era5 query was fetched already, skipping...')

def test_era5_fetch_wavedirection():
    if not Era5().fetch_wavedirection(start=start, end=end):
        print('era5 query was fetched already, skipping...')

def test_era5_fetch_waveperiod():
    if not Era5().fetch_waveperiod(start=start, end=end):
        print('era5 query was fetched already, skipping...')

def test_era5_fetch_wind():
    if not Era5().fetch_wind(start=start, end=end):
        print('era5 query was fetched already, skipping...')

def test_era5_load_windwaveswellheight():
    height, lat, lon, time = Era5().load_windwaveswellheight(
            south=south, north=north, west=west, east=east, start=start, end=end
        )
    assert(len(height)==len(lat)==len(lon))
    assert(len(lat) > 0)

def test_era5_load_wavedirection():
    wave, lat, lon, time = Era5().load_wavedirection(
            south=south, north=north, west=west, east=east, start=start, end=end
        )
    assert(len(wave)==len(lat)==len(lon))
    assert(len(lat) > 0)

def test_era5_load_waveperiod():
    wave, lat, lon, time = Era5().load_waveperiod(
            south=south, north=north, west=west, east=east, start=start, end=end
        )
    assert(len(wave)==len(lat)==len(lon))
    assert(len(lat) > 0)

def test_era5_load_wind():
    wind, lat, lon, time = Era5().load_wind(
            south=south, north=south+1, west=west, east=west+1, start=start, end=end
        )
    assert(len(wind)==len(lat)==len(lon))
    assert(len(lat) > 0)


