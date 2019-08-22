import numpy as np
from datetime import datetime, timedelta
import os
import urllib.request
import pygrib

from kadlu.geospatial.data_sources import fetch_util

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


def fetch(wavevar=waveSources['swh'], time=datetime.now(), region=regions[0]):
    storage_location = fetch_util.instantiate_storage_config() 

    level_ref = 'SFC_0'
    if wavevar is 'WVDIR': level_ref = 'TGL_0'
    hour = f"{((time.hour % 24) // 6 * 6):02d}"  # better option?
    predictionhour = '000'  # any better than 0-hour?
    fetchname = f"CMC_rdwps_{region}_{wavevar}_{level_ref}_latlon0.05x0.05_{time.strftime('%Y%m%d')}{hour}_P{predictionhour}.grib2"
    fetchfile= f"{storage_location}{fetchname}"

    if os.path.isfile(fetchfile):
        print("File exists, skipping retrieval...")
        # obsolete ???
        #validate_wavesource(fetchfile, waveSources)
    else:
        print("Downloading file from the Regional Deterministic Wave Prediction System...")
        fetchurl = f"http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/{hour}/{fetchname}"
        urllib.request.urlretrieve(fetchurl, fetchfile)


def load(filepath, plot=False): return fetch_util.loadgrib(filepath, plot)

"""
    # If no gribfile argument is provided, default to the fetched file.
    if grib is None:
        grib = self.fetch_filename

    # If no date argument is provided, default to the module timestamp, truncated to nearest earlier 6-hour interval.
    if target_date is None:

        rep_hour = ((self.fetch_datetimestamp.hour % 24) // 6) * 6
        target_date = self.fetch_datetimestamp.replace(minute=0, hour=rep_hour, second=0, microsecond=0)
    else:
### Should similar filtering on time be added here to enforce 6 hour intervals?
        target_date = target_date

    # load grib structure from target.
    grbs=pygrib.open(grib)

    # Fetch slice of data corresponding to selected target date and wave variable.
        #swh	- significant height of combined wind waves and swell, metres (HTSGW)
        #mwp	- primary wave mean period, seconds (PKPER)
        #mwd	- primary wind wave direction, degrees true (i.e. 0 deg = North; proceeding clockwise) (WVDIR)
    grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]
    
    if(wavevar == RDWPSWavevar.HTSGW):
        title_text = "CMC-RDWPS Sig. Wave + Swell Height from GRIB\n({}) ".format(target_date)
    elif(wavevar == RDWPSWavevar.PKPER):
        title_text = "CMC-RDWPS Mean Wave Period from GRIB\n({}) ".format(target_date)
    elif(wavevar == RDWPSWavevar.WVDIR):
        title_text = "CMC-RDWPS Mean Wave Direction from GRIB\n({}) ".format(target_date)
    else:
        title_text = "CMC-RDWPS\nUnknown variable from GRIB\n({}) ".format(target_date)

    return (grb, title_text)
"""
    
