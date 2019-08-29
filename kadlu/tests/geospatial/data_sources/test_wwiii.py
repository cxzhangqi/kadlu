import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import wwiii
from kadlu.geospatial.data_sources import fetch_util

def test_fetch_wwiii():
    wwiii.fetch_var(wavevar=wwiii.waveSources['swh'], time=datetime(2017, 2, 3, 0, 0, 0, 0), region=wwiii.regions['global'])

def test_plot_wwiii():
    storage_location = fetch_util.instantiate_storage_config()
    #filename = "multi_1.glo_30m.hs.201702.grb2"
    fname = wwiii.fetchfile()
    filepath = f"{storage_location}{fname}"
    data = wwiii.load(filepath=filepath, plot=True)

def test_load_wwiii():
    wave, lats, lons = wwiii.load(south=43, west=-60, north=44, east=-59)
   
def test_fetch_wwiii_fetchall():
    for output in wwiii.fetch:
        datafiles = output(south, north, west, east, time=datetime(2017, 2, 3, 0, 0, 0, 0))

