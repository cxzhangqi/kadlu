import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import rdwps

def test_fetch_rdwps():
    from kadlu.geospatial.data_sources import fetch_util
    storage_location = fetch_util.instantiate_storage_config()
    rdwps.fetch(storage_location=storage_location, wavevar=rdwps.waveSources['swh'], time=datetime.now()-timedelta(hours=3), region=rdwps.regions[0])
