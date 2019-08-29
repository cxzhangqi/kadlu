"""
https://weather.gc.ca/grib/grib2_RDWPS_e.html
"""

import numpy as np
from datetime import datetime, timedelta
import os
import urllib.request
import pygrib
from kadlu.geospatial.data_sources import fetch_util
from collections import defaultdict

# dictionary to obtain reference level as per RDWPS nomenclature
# for more info see the following: https://weather.gc.ca/grib/grib2_RDWPS_e.html
regions = ['gulf-st-lawrence', 'superior', 'huron-michigan', 'erie', 'ontario']
default = defaultdict(lambda : 'SFC_0', key='default')
level_ref = {key : default for key in regions}
for reg in regions: 
    level_ref[reg]['UGRD'] = 'TGL_10'
    level_ref[reg]['VGRD'] = 'TGL_10'
level_ref['gulf-st-lawrence']['PKPER'] = 'TGL_0'
level_ref['gulf-st-lawrence']['PRMSL'] = 'MSL_0'


def fetch_rdwps(wavevar, time, regions, level_ref=level_ref):
    storage_location = fetch_util.instantiate_storage_config() 
    fetchfiles = []

    for region in regions:
        fname = fetchname(wavevar, time, region, level_ref)
        fetchfile= f"{storage_location}{fname}"
        if os.path.isfile(fetchfile):
            print(f"File {fname} exists, skipping retrieval...")
        else:
            print("Downloading from the Regional Deterministic Wave Prediction System...")
            fetchurl = f"http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/{((time.hour % 24) // 6 * 6):02d}/{fname}"
            urllib.request.urlretrieve(fetchurl, fetchfile)
        fetchfiles.append(fetchfile)

    return fetchfiles


def fetchname(wavevar, time, region, level_ref=level_ref):
    hour = f"{((time.hour % 24) // 6 * 6):02d}"  # better option?
    predictionhour = '000'  # any better than 0-hour?
    return f"CMC_rdwps_{region}_{wavevar}_{level_ref[region][wavevar]}_latlon0.05x0.05_{time.strftime('%Y%m%d')}{hour}_P{predictionhour}.grib2"


def abstract_region(south, north, east, west):
    # this function will eventually return a list of regions determined by the input boundaries
    regions = [
        'gulf-st-lawrence',
        'superior',
        'huron-michigan',
        'erie',
        'ontario'
    ]
    return


def fetch_windwaveswellheight(south=-90, north=90, east=-180, west=180, time=datetime.now()): return fetch_rdwps('HTSGW', time, abstract_region(south, north, east, west))
def fetch_windwaveheight    (south=-90, north=90, east=-180, west=180, time=datetime.now()): return fetch_rdwps('WVHGT', time, abstract_region(south, north, east, west))
def fetch_wavedirection     (south=-90, north=90, east=-180, west=180, time=datetime.now()): return fetch_rdwps('WVDIR', time, abstract_region(south, north, east, west))
def fetch_waveperiod        (south=-90, north=90, east=-180, west=180, time=datetime.now()): return fetch_rdwps('WVPER', time, abstract_region(south, north, east, west))
def fetch_wind_u            (south=-90, north=90, east=-180, west=180, time=datetime.now()): return fetch_rdwps('UGRD', time, abstract_region(south, north, east, west))
def fetch_wind_v            (south=-90, north=90, east=-180, west=180, time=datetime.now()): return fetch_rdwps('VGRD', time, abstract_region(south, north, east, west))
def fetch_icecover          (south=-90, north=90, east=-180, west=180, time=datetime.now()): return fetch_rdwps('ICEC', time, abstract_region(south, north, east, west))

fetch = [fetch_windwaveheight, fetch_windwaveswellheight, fetch_wavedirection, fetch_waveperiod, fetch_wind_u, fetch_wind_v, fetch_icecover]

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
    
