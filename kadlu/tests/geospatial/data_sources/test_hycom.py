import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import hycom
from kadlu.geospatial.data_sources import fetch_util

def test_fetch_hycom():
    salinity, lat, lon = hycom.fetch(south=43, west=-60, north=44, east=-59)

def test_load_hycom():
    salinity, lat, lon = hycom.load(south=43, west=-60, north=44, east=-59)
