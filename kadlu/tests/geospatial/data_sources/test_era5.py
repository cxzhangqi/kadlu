import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import era5
from kadlu.geospatial.data_sources import fetch_util

def test_fetch_era5_waveheight():
    era5.fetch_waveheight(time=datetime(2018, 1, 1, 0, 0, 0, 0))

def test_fetch_era5_wavedirection():
    era5.fetch_wavedirection(time=datetime(2018, 1, 1, 0, 0, 0, 0))

def test_fetch_era5_waveperiod():
    era5.fetch_waveperiod(time=datetime(2018, 1, 1, 0, 0, 0, 0))

def test_plot_era5_waveheight():
    storage_location = fetch_util.instantiate_storage_config()
    filename = "ERA5_reanalysis_significant_height_of_combined_wind_waves_and_swell_2018-01-01_00h.grb2"
    filepath = f"{storage_location}{filename}"
    data = era5.load(filepath=filepath, plot=True)

"""
mahone bay test area:
north =  44.7
south =  44.4
west  = -64.4
east  = -63.8
"""

def test_fetch_era5_fetchall():
    # for fcn in era5.fetch:
    #     datafiles = fcn(time=datetime(2018, 1, 1, 0, 0, 0, 0))
    datafiles = map(lambda fcn : fcn(time=datetime(2018, 1, 1, 0, 0, 0, 0)), era5.fetch)

