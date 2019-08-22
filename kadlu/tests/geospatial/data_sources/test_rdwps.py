import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import rdwps
from kadlu.geospatial.data_sources import fetch_util

def test_fetch_rdwps():
    rdwps.fetch(wavevar=rdwps.waveSources['swh'], time=datetime.now()-timedelta(hours=3), region=rdwps.regions[0])

def test_plot_rdwps():
    storage_location = fetch_util.instantiate_storage_config()

    # build the filename dynamically for today's data
    time=datetime.now()-timedelta(hours=3)
    region = 'gulf-st-lawrence'
    wavevar = 'HTSGW'
    level_ref = 'SFC_0'
    hour = f"{((time.hour % 24) // 6 * 6):02d}"
    predictionhour = '000'
    filename = f"CMC_rdwps_{region}_{wavevar}_{level_ref}_latlon0.05x0.05_{time.strftime('%Y%m%d')}{hour}_P{predictionhour}.grib2"
    filepath = f"{storage_location}{filename}"

    data = rdwps.load(filepath=filepath, plot=True)
