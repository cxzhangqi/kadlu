import numpy as np
from datetime import datetime, timedelta
import os
import cdsapi
import pygrib

from kadlu.geospatial.data_sources import fetch_util

waveSources = {
    'swh' : 'significant_height_of_combined_wind_waves_and_swell',
    'mwd' : 'mean_wave_direction',
    'mwp' : 'mean_wave_period'
}


def fetch_era5(wavevar, time):
    storage_location = fetch_util.instantiate_storage_config() 
    fetchfile = f"{storage_location}ERA5_reanalysis_{wavevar}_{time.strftime('%Y-%m-%d_%Hh')}.grb2"

    if os.path.isfile(fetchfile):
        print("File exists, skipping retrieval...")
        # obsolete ???
        #validate_wavesource(fetchfile, waveSources)
    else:
        print(f"Downloading {wavevar} from Copernicus Climate Data Store...")
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

    return fetchfile

def fetch_waveheight(time=datetime.now()):      return fetch_era5(waveSources['swh'], time)
def fetch_wavedirection(time=datetime.now()):   return fetch_era5(waveSources['mwd'], time)
def fetch_waveperiod(time=datetime.now()):      return fetch_era5(waveSources['mwp'], time)


def load(filepath, plot=False): return fetch_util.loadgrib(filepath, plot)


"""
    # If no gribfile argument is provided, default to the fetched file.
    if grib is None:
        grib = self.fetch_filename

    # load grib structure from target.
    grbs=pygrib.open(grib)

    # Identify parameter and date for extraction.
    # Date field: validDate
    # Target date (not incl. time yet)
    if (target_date is None):
        target_date = self.fetch_datetimestamp.replace(minute=0, hour=0, second=0, microsecond=0)
    else:
### Should any filtering on time be added here to enforce valid time intervals?
        target_date = target_date

    # Fetch the indicated slice from the overall Grib file.
    grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]

    if(wavevar == ERA5Wavevar.significant_height_of_combined_wind_waves_and_swell):
        title_text = "ERA5 Sig. Wave + Swell Height from GRIB\n({}) ".format(target_date)
    elif(wavevar == ERA5Wavevar.mean_wave_period):
        title_text = "ERA5 Mean Wave Period from GRIB\n({}) ".format(target_date)
    elif(wavevar == ERA5Wavevar.mean_wave_direction):
        title_text = "ERA5 Mean Wave Direction from GRIB\n({}) ".format(target_date)
    else:
        title_text = "ERA5\nUnknown variable from GRIB\n({}) ".format(target_date)

    return (grb, title_text)
"""

