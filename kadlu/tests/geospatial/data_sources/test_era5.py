import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import era5
from kadlu.geospatial.data_sources import fetch_util

def test_fetch_era5():
    storage_location = fetch_util.instantiate_storage_config()
    era5.fetch(storage_location=storage_location, wavevar=era5.waveSources['swh'], time=datetime(2018, 1, 1, 0, 0, 0, 0))

def test_plot_era5():
    storage_location = fetch_util.instantiate_storage_config()
    filename = "ERA5_reanalysis_significant_height_of_combined_wind_waves_and_swell_2018-01-01_00h.grb2"
    filepath = f"{storage_location}{filename}"
    data = era5.load(filepath=filepath, plot=True)
