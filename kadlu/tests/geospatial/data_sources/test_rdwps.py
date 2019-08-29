import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import rdwps
from kadlu.geospatial.data_sources import fetch_util

def test_fetch_rdwps_waveswellheight():
    rdwps.fetch_windwaveswellheight(south=-90, north=90, west=180, east=-180, time=datetime.now()-timedelta(hours=3))

def test_fetch_rdwps_fetchall():
    for output in rdwps.fetch:
        datafiles = output(south=-90, north=90, west=180, east=-180, time=datetime.now()-timedelta(hours=3))

def test_rdwps_region_abstraction():
    regions = rdwps.abstract_region(south=-90, north=90, west=180, east=-180)
    assert(len(regions) == 5)
    # TODO: add more assertions testing the region output for each: 
    # gulf-st-lawrence, erie, ontario, huron-michigan, superior

def test_plot_rdwps():
    storage_location = fetch_util.instantiate_storage_config()

    # filename for today's wave height data in the gulf of st lawrence
    time=datetime.now()-timedelta(hours=3)
    region = 'gulf-st-lawrence'
    wavevar = 'HTSGW'
    fname = rdwps.fetchname(wavevar, time, region)
    filepath = f"{storage_location}{fname}"

    data = rdwps.load(filepath=filepath, plot=True)

