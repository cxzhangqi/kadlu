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

def fetch(storage_location, wavevar=waveSources['swh'], time=datetime.now(), region=regions['global']):
    fetchname = f"multi_1.{region}.{wavevar}.{time.strftime('%Y%m')}.grb2"
    fetchfile = f"{storage_location}{fetchname}"

    if os.path.isfile(fetchfile):
        print("File exists, skipping retrieval...")
        # obsolete ???
        #validate_wavesource(fetchfile, waveSources)
    else:
        print("Downloading file from NOAA WaveWatch III...")
        fetchurl = f"https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/{time.strftime('%Y/%m')}/gribs/{fetchname}"
        urllib.request.urlretrieve(fetchurl, fetchfile)


def load(filepath, plot=False):
    grib = pygrib.open(filepath)

    if plot: fetch_util.plotSampleGrib(grib[1], "testing")

    #data = [None for msg in range(grib.messages)]
    #for x in range(0, grib.messages):
    #    data[x] = grib[x+1].data()  # grib indexing starts at 1 for some reason

    #return np.array(data)
    return grib[1]

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

