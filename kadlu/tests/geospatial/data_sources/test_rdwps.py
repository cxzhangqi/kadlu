import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import rdwps 
from kadlu.geospatial.data_sources.rdwps import Rdwps
from kadlu.geospatial.data_sources import fetch_util
from os.path import isfile

time = datetime.now() - timedelta(hours=3)

def test_rdwps_fetch_windwaveswellheight():
    filenames = Rdwps().fetch_windwaveswellheight(south=-90, north=90, west=-180, east=180, time=time)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_windwaveheight():
    filenames = Rdwps().fetch_windwaveheight(south=-90, north=90, west=-180, east=180, time=time)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_wavedirection():
    filenames = Rdwps().fetch_wavedirection(south=-90, north=90, west=-180, east=180, time=time)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_waveperiod():
    filenames = Rdwps().fetch_waveperiod(south=-90, north=90, west=-180, east=180, time=time)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_wind_u():
    filenames = Rdwps().fetch_wind_u(south=-90, north=90, west=-180, east=180, time=time)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_wind_v():
    filenames = Rdwps().fetch_wind_v(south=-90, north=90, west=-180, east=180, time=time)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_fetch_icecover():
    filenames = Rdwps().fetch_icecover(south=-90, north=90, west=-180, east=180, time=time)
    for fname in filenames:
        assert(isfile(fname))

def test_rdwps_region_abstraction():
    regions = rdwps.ll_2_regionstr(south=-90, north=90, west=-180, east=180)
    assert(len(regions) == 5)
    # TODO: add more assertions testing a single region boundary for each: 
    #       gulf-st-lawrence, erie, ontario, huron-michigan, superior

#def test_rdwps_plot():
#    filenames = Rdwps().fetch_windwaveswellheight(south=-90, north=90, west=-180, east=180, time=time)
#    data = Rdwps().load(filenames, plot=f"Wind Wave Swell Height {time.strftime('%Y-%m')}" )

def test_rdwps_load_icecover():
    ice, lat, lon = Rdwps().load_icecover(south=-90, north=90, west=-180, east=180, time=datetime(2018, 2, 1), plot="yep it works")

