import numpy as np
from datetime import datetime, timedelta
import os
import urllib.request

#from kadlu.geospatial.data_sources import fetch_util
from kadlu.geospatial.data_sources.fetch_util import validate_wavesource

waveSources = {
    'swh' : 'HTSGW',
    'mwd' : 'WVDIR',
    'mwp' : 'PKPER'
}

regions = [
    # matt_s 2019-08
    # this will probably change in the future to something other than a list...
    'gulf-st-lawrence',
    'superior',
    'huron-michigan',
    'erie',
    'ontario'
]


def fetch(storage_location, wavevar=waveSources['swh'], time=datetime.now(), region=regions[0]):
    level_ref = 'SFC_0'
    if wavevar is 'WVDIR': level_ref = 'TGL_0'
    hour = f"{((time.hour % 24) // 6 * 6):02d}"  # better option?
    predictionhour = '000'  # any better than 0-hour?
    fetchname = f"CMC_rdwps_{region}_{wavevar}_{level_ref}_latlon0.05x0.05_{time.strftime('%Y%m%d')}{hour}_P{predictionhour}.grib2"
    fetchfile= f"{storage_location}{fetchname}"

    if os.path.isfile(fetchfile):
        print("File exists, skipping retrieval...")
        validate_wavesource(fetchfile, waveSources)
    else:
        print("Downloading file from the Regional Deterministic Wave Prediction System...")
        fetchurl = f"http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/{hour}/{fetchname}"
        urllib.request.urlretrieve(fetchurl, fetchfile)

def load():
    pass

