import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import rdwps as rdwps_module
from kadlu.geospatial.data_sources.rdwps import Rdwps
from kadlu.geospatial.data_sources import fetch_util

rdwps = Rdwps()
storage_location = fetch_util.storage_cfg()

def test_fetch_rdwps_windwaveswellheight():
    test = rdwps.fetch_windwaveswellheight(south=-90, north=90, west=180, east=-180, time=datetime.now()-timedelta(hours=3))

def test_fetch_rdwps_fetchall():
    for fcn in rdwps.fetch_functions:
        fcn(rdwps, south=-90, north=90, west=180, east=-180, time=datetime.now()-timedelta(hours=3))

def test_rdwps_region_abstraction():
    regions = rdwps_module.abstract_region(south=-90, north=90, west=180, east=-180)
    assert(len(regions) == 5)
    # TODO: add more assertions testing a single region boundary for each: 
    #       gulf-st-lawrence, erie, ontario, huron-michigan, superior

def test_plot_rdwps():
    time = datetime.now() - timedelta(hours=3)
    fnames = rdwps.fetch_windwaveswellheight(south=-90, north=90, west=180, east=-180, time=time)
    #data = rdwps.load(filepath=fnames[0], plot="Wind Wave Swell Height")

