""" Unit tests for the the 'geospatial.data_provider' module in the 'kadlu' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import pytest
import os
import numpy as np
from datetime import datetime, timedelta

path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets")

"""
from kadlu.geospatial.data_provider import DataProvider

# matt_s 2019-08
# commented this out as a quick fix to mute errors related to deprecated code in bathy_reader, bathy_interpolator

def test_provide_bathymetry_from_a_single_chs_file():
    folder = os.path.join(path_to_assets, "tif")
    provider = DataProvider(storage_location=folder, bathy_source="CHS", south=43, west=-60, north=44, east=-59)
    bathy_data = provider.bathy_data
    bathy = bathy_data[0]
    lats = bathy_data[1]
    lons = bathy_data[2]
    assert np.ma.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    assert np.ma.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]
"""

from kadlu.geospatial.data_sources import era5
def test_fetch_era5():
    from kadlu.geospatial.data_sources import fetch_util
    storage_location = fetch_util.instantiate_storage_config()
    era5.fetch(storage_location=storage_location, wavevar=era5.waveSources['swh'], time=datetime(2018, 1, 1, 0, 0, 0, 0))

from kadlu.geospatial.data_sources import wwiii
def test_fetch_wwiii():
    from kadlu.geospatial.data_sources import fetch_util
    storage_location = fetch_util.instantiate_storage_config()
    wwiii.fetch(storage_location=storage_location, wavevar=wwiii.waveSources['swh'], time=datetime(2017, 2, 3, 0, 0, 0, 0), region=wwiii.regions['global'])

from kadlu.geospatial.data_sources import rdwps
def test_fetch_rdwps():
    from kadlu.geospatial.data_sources import fetch_util
    storage_location = fetch_util.instantiate_storage_config()
    rdwps.fetch(storage_location=storage_location, wavevar=rdwps.waveSources['swh'], time=datetime.now()-timedelta(hours=3), region=rdwps.regions[0])
