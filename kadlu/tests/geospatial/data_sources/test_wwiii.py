import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import wwiii

def test_fetch_wwiii():
    from kadlu.geospatial.data_sources import fetch_util
    storage_location = fetch_util.instantiate_storage_config()
    wwiii.fetch(storage_location=storage_location, wavevar=wwiii.waveSources['swh'], time=datetime(2017, 2, 3, 0, 0, 0, 0), region=wwiii.regions['global'])

