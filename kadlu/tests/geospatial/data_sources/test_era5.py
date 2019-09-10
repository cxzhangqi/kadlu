import pytest
import numpy as np
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import era5
from kadlu.geospatial.data_sources.era5 import Era5
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
from os.path import isfile

# time used for testing
time = datetime(2018, 1, 1, 0, 0, 0, 0)

# disable fetch testing when not necessary to avoid slamming the API with automated requests
test_fetch = True

def test_era5_fetch_windwaveswellheight():
    if not test_fetch: return
    filenames = Era5().fetch_windwaveswellheight(time)
    for fname in filenames:
        assert(isfile(fname))

def test_era5_fetch_wavedirection():
    if not test_fetch: return
    filenames = Era5().fetch_wavedirection(time)
    for fname in filenames:
        assert(isfile(fname))

def test_era5_fetch_waveperiod():
    if not test_fetch: return
    filenames = Era5().fetch_waveperiod(time)
    for fname in filenames:
        assert(isfile(fname))

def test_era5_load_windwaveswellheight():
    height, lat, lon = Era5().load_windwaveswellheight(time)
    assert(len(height)==len(lat)==len(lon))

def test_era5_load_wavedirection():
    wave, lat, lon = Era5().load_wavedirection(time)
    assert(len(wave)==len(lat)==len(lon))

def test_era5_load_waveperiod():
    wave, lat, lon = Era5().load_waveperiod(time)
    assert(len(wave)==len(lat)==len(lon))

def test_plot_era5_windwaveswellheight():
    filenames = Era5().fetch_windwaveswellheight(time)
    fname = filenames[0]
    for fname in filenames:
        assert(isfile(fname))
    wave, lat, lon = Era5().load_windwaveswellheight(time, plot="Sig. Height of Combined Wind, Wave, and Swell 2018-01-01")
    assert(len(wave) == len(lat) == len(lon))

