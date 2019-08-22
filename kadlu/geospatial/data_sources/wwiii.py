import numpy as np
from datetime import datetime, timedelta
import os
import urllib.request
import pygrib

from kadlu.geospatial.data_sources import fetch_util

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

def fetchname(wavevar, time, region): return f"multi_1.{region}.{wavevar}.{time.strftime('%Y%m')}.grb2"

def fetch(wavevar=waveSources['swh'], time=datetime.now(), region=regions['global']):
    """
        matt_s 2019-08
        Note that WWIII returns data from an entire month, regardless of
        the time given. In the future the load function should be updated
        to find the desired message within the grib file containing data for
        the desired date
    """
    storage_location = fetch_util.instantiate_storage_config() 
    fetchname = fetchname(wavevar, time, region)
    fetchfile = f"{storage_location}{fetchname}"

    if os.path.isfile(fetchfile):
        print("File exists, skipping retrieval...")
        # obsolete ???
        #validate_wavesource(fetchfile, waveSources)
    else:
        print("Downloading file from NOAA WaveWatch III...")
        fetchurl = f"https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/{time.strftime('%Y/%m')}/gribs/{fetchname}"
        urllib.request.urlretrieve(fetchurl, fetchfile)


    # matt_s 2019-08
    # the loadgrib function returns data for the whole month for the
    # WWIII source only. This should be updated here to parse grib
    # messages from the entire month, and return the 3-hour slice
    # with the desired data. 
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
    # DEBUG        date_valid = datetime(2018,11,2)
    if (target_date is None):
        target_date = self.fetch_datetimestamp.replace(minute=0, hour=0, second=0, microsecond=0)
    else:
    ### Should any filtering on time be added here to enforce valid time intervals?
        target_date = target_date

    # Fetch the indicated slice from the overall Grib file.
    grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]

    if(wavevar == WWIIIWavevar.hs):
        title_text = "WW III Sig. Wave + Swell Height from GRIB\n({}) ".format(target_date)
    elif(wavevar == WWIIIWavevar.tp):
        title_text = "WW III Mean Wave Period from GRIB\n({}) ".format(target_date)
    elif(wavevar == WWIIIWavevar.dp):
        title_text = "WW III Mean Wave Direction from GRIB\n({}) ".format(target_date)
    else:
        title_text = "WW III\nUnknown variable from GRIB\n({}) ".format(target_date)

    return (grb, title_text)
"""

