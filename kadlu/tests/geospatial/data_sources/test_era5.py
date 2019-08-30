import pytest
from datetime import datetime, timedelta
#from kadlu.geospatial.data_sources import era5
from kadlu.geospatial.data_sources.era5 import Era5
from kadlu.geospatial.data_sources.fetch_util import storage_cfg

# time used for testing
time = datetime(2018, 1, 1, 0, 0, 0, 0)

def test_era5_fetch_waveheight():
    Era5().fetch_waveheight(time)

def test_era5_fetch_wavedirection():
    Era5().fetch_wavedirection(time)

def test_era5_fetch_waveperiod():
    Era5().fetch_waveperiod(time)

def test_era5_load_waveheight():
    pass

def test_era5_load_wavedirection():
    pass

def test_era5_load_waveperiod():
    pass

def test_plot_era5_waveheight():
    filename = "ERA5_reanalysis_significant_height_of_combined_wind_waves_and_swell_2018-01-01_00h.grb2"
    filepath = f"{storage_cfg()}{filename}"
    data = Era5().load(filepath=filepath, plot="Wave Height 2018-01-01")

