import pytest
import numpy as np
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import era5
from kadlu.geospatial.data_sources.era5 import Era5
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
from os.path import isfile

source = Era5()

# gulf st lawrence
south, west = 40, -70
north, east = 59, -55

start = datetime(2018, 1, 1, 0, 0, 0, 0)
end   = datetime(2018, 1, 1, 0, 0, 0, 0)

# note: see assertions within fetch function

def test_era5_fetch_windwaveswellheight():
    filenames = source.fetch_windwaveswellheight(start=start, end=end)

def test_era5_fetch_wavedirection():
    filenames = source.fetch_wavedirection(start=start, end=end)

def test_era5_fetch_waveperiod():
    filenames = source.fetch_waveperiod(start=start, end=end)

def test_era5_load_windwaveswellheight():
    height, lat, lon, time = source.load_windwaveswellheight(south=south, north=north, west=west, east=east,
            start=start, end=end)
    assert(len(height)==len(lat)==len(lon))
    assert(len(lat) > 0)

def test_era5_load_wavedirection():
    wave, lat, lon, time = source.load_wavedirection(south=south, north=north, west=west, east=east,
            start=start, end=end)
    assert(len(wave)==len(lat)==len(lon))
    assert(len(lat) > 0)

def test_era5_load_waveperiod():
    wave, lat, lon, time = source.load_waveperiod(south=south, north=north, west=west, east=east,
            start=start, end=end)
    assert(len(wave)==len(lat)==len(lon))
    assert(len(lat) > 0)

