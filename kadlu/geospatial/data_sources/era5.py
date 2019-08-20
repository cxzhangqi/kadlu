import numpy as np
from datetime import datetime, timedelta
import os
import cdsapi

#from kadlu.geospatial.data_sources import fetch_util
from kadlu.geospatial.data_sources.fetch_util import validate_wavesource

waveSources = {
    'swh' : 'significant_height_of_combined_wind_waves_and_swell',
    'mwd' : 'mean_wave_direction',
    'mwp' : 'mean_wave_period'
}


def fetch(storage_location, wavevar=waveSources['swh'], time=datetime.now()):
    fetchfile = f"{storage_location}ERA5_reanalysis_{wavevar}_{time.strftime('%Y-%m-%d_%Hh')}.grb2"

    if os.path.isfile(fetchfile):
        print("File exists, skipping retrieval...")
        validate_wavesource(fetchfile, waveSources)
    else:
        print("Downloading file from Copernicus Climate Data Store...")
        c = cdsapi.Client()
        c.retrieve('reanalysis-era5-single-levels', {
            'product_type'  : 'reanalysis',
            'format'        : 'grib',
            'variable'      : wavevar,
            'year'          : time.strftime("%Y"),
            'month'         : time.strftime("%m"),
            'day'           : time.strftime("%d"),
            'time'          : time.strftime("%H:%M")
        }, fetchfile)


def load():
    pass

