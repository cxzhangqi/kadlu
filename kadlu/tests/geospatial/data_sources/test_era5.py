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

def test_plot_era5():
    storage_location = fetch_util.instantiate_storage_config()
    filename = "ERA5_reanalysis_significant_height_of_combined_wind_waves_and_swell_2018-01-01_00h.grb2"
    filepath = f"{storage_location}{filename}"
    data = era5.load(filepath=filepath, plot=True)

"""
mahone bay test area:
Lat: north: 44.7; south: 44.4
Long: west: -64.4; east: -63.8
"""

"""
def test_fetch_era5_refactored():
    for output in era5.fetch:
        print(f"Testing {output.__name__}")
        datafiles = output(south=44.4, north=44.7, west-64.4, east-63.8, time=datetime(2018, 1, 1, 0, 0, 0, 0))
"""

