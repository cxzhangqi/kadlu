import pytest
import numpy as np
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources.era5 import Era5
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
from os.path import isfile

# time used for testing
time = datetime(2018, 1, 1, 0, 0, 0, 0)

def test_era5_fetch_windwaveswellheight():
    filenames = Era5().fetch_windwaveswellheight(time)
    for fname in filenames:
        assert(isfile(fname))

def test_era5_fetch_wavedirection():
    filenames = Era5().fetch_wavedirection(time)
    for fname in filenames:
        assert(isfile(fname))

def test_era5_fetch_waveperiod():
    filenames = Era5().fetch_waveperiod(time)
    for fname in filenames:
        assert(isfile(fname))

def test_era5_load_windwaveswellheight():
    pass

def test_era5_load_wavedirection():
    pass

def test_era5_load_waveperiod():
    pass

def test_plot_era5_windwaveswellheight():
    filenames = Era5().fetch_windwaveswellheight(time)
    fname = filenames[0]
    for fname in filenames:
        assert(isfile(fname))
    wave, lat, lon = Era5().load(fname, plot="Sig. Height of Combined Wind, Wave, and Swell 2018-01-01")
    #wave, lat, lon = Era5().load(filenames, plot=False)
    assert(len(wave) == len(lat) == len(lon))

