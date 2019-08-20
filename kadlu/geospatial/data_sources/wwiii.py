import numpy as np
from datetime import datetime, timedelta
import os
import urllib.request

#from kadlu.geospatial.data_sources import fetch_util
from kadlu.geospatial.data_sources.fetch_util import validate_wavesource

waveSources = {
    'swh' : 'hs',
    'mwd' : 'dp',
    'mwp' : 'tp'
}

regions = {
    # This may be replaced by a more modular mapping system later
    # depending on how the API works for [region]_[interval] formatting
    'global'    : 'glo_30m',
    'arctic'    : 'ao_30m',
    'pacific'   : 'ep_10m',
    'atlantic'  : {'10m' : 'at_10m', '4m' : 'at_4m'},
    'US west'   : {'10m' : 'wc_10m', '4m' : 'wc_4m'},
    'alaska'    : {'10m' : 'ak_10m', '4m' : 'ak_4m'}
}

def fetch(storage_location, wavevar=waveSources['swh'], time=datetime.now(), region=regions['global']):
    fetchname = f"multi_1.{region}.{wavevar}.{time.strftime('%Y%m')}.grb2"
    fetchfile = f"{storage_location}{fetchname}"

    if os.path.isfile(fetchfile):
        print("File exists, skipping retrieval...")
        validate_wavesource(fetchfile, waveSources)
    else:
        print("Downloading file from NOAA WaveWatch III...")
        fetchurl = f"https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/{time.strftime('%Y/%m')}/gribs/{fetchname}"
        urllib.request.urlretrieve(fetchurl, fetchfile)


def load():
    pass

