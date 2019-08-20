import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import era5

def test_fetch_era5():
    from kadlu.geospatial.data_sources import fetch_util
    storage_location = fetch_util.instantiate_storage_config()
    era5.fetch(storage_location=storage_location, wavevar=era5.waveSources['swh'], time=datetime(2018, 1, 1, 0, 0, 0, 0))

