"""
    API for Regional Deterministic Wave Prediction System (RDWPS) data provided by Govt. of Canada

    Metadata regarding the dataset can be found here:
        https://weather.gc.ca/grib/grib2_RDWPS_e.html

    matt_s 2019-08
"""

import numpy as np
from datetime import datetime, timedelta
import os
import urllib.request
import pygrib
from kadlu.geospatial.data_sources import fetch_util 
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
from collections import defaultdict
import warnings

# dictionary to obtain reference level as per RDWPS nomenclature
# for more info see the following: https://weather.gc.ca/grib/grib2_RDWPS_e.html
regions = ['gulf-st-lawrence', 'superior', 'huron-michigan', 'erie', 'ontario']
default = defaultdict(lambda : 'SFC_0')
level_ref = {key : default for key in regions}
for reg in regions: 
    level_ref[reg]['UGRD'] = 'TGL_10'
    level_ref[reg]['VGRD'] = 'TGL_10'
level_ref['gulf-st-lawrence']['PKPER'] = 'TGL_0'
level_ref['gulf-st-lawrence']['PRMSL'] = 'MSL_0'


def fetch_rdwps(wavevar, time, regions, level_ref=level_ref):
    fetchfiles = []
    for region in regions:
        fname = fetchname(wavevar, time, region, level_ref)
        fetchfile = f"{storage_cfg()}{fname}"
        fetchfiles.append(fetchfile)
        if os.path.isfile(fetchfile):
            print(f"File {fname} exists, skipping retrieval...")
        else:
            print(f"Downloading {fname} from the Regional Deterministic Wave Prediction System...")
            fetchurl = f"http://dd.weather.gc.ca/model_wave/ocean/gulf-st-lawrence/grib2/{((time.hour % 24) // 6 * 6):02d}/{fname}"
            urllib.request.urlretrieve(fetchurl, fetchfile)

    return fetchfiles


def fetchname(wavevar, time, region, level_ref=level_ref):
    hour = f"{((time.hour % 24) // 6 * 6):02d}"  # better option?
    predictionhour = '000'  # any better than 0-hour?
    return f"CMC_rdwps_{region}_{wavevar}_{level_ref[region][wavevar]}_latlon0.05x0.05_{time.strftime('%Y%m%d')}{hour}_P{predictionhour}.grib2"


def abstract_region(south, north, west, east):
    # this function will eventually return a list of regions determined by the input boundaries
    warnings.warn("RDWPS region abstraction function is incomplete. Instead, you get the gulf of st lawrence")
    regions = [
        'gulf-st-lawrence',
        'superior',
        'huron-michigan',
        'erie',
        'ontario'
    ]
    return ['gulf-st-lawrence']

class Rdwps(): 
    def fetch_windwaveswellheight(self, south=-90, north=90, west=180, east=-180, time=datetime.now()): 
        return fetch_rdwps('HTSGW', time, abstract_region(south, north, west, east))

    def fetch_windwaveheight(self, south=-90, north=90, west=180, east=-180, time=datetime.now()):
        return fetch_rdwps('WVHGT', time, abstract_region(south, north, west, east))

    def fetch_wavedirection(self, south=-90, north=90, west=180, east=-180, time=datetime.now()):
        return fetch_rdwps('WVDIR', time, abstract_region(south, north, west, east))

    def fetch_waveperiod(self, south=-90, north=90, west=180, east=-180, time=datetime.now()):
        return fetch_rdwps('WVPER', time, abstract_region(south, north, west, east))

    def fetch_wind_u(self, south=-90, north=90, west=180, east=-180, time=datetime.now()):
        return fetch_rdwps('UGRD', time, abstract_region(south, north, west, east))

    def fetch_wind_v(self, south=-90, north=90, west=180, east=-180, time=datetime.now()):
        return fetch_rdwps('VGRD', time, abstract_region(south, north, west, east))

    def fetch_icecover(self, south=-90, north=90, west=180, east=-180, time=datetime.now()): 
        return fetch_rdwps('ICEC', time, abstract_region(south, north, west, east))

    header = "(south=-90, north=90, east=-180, west=180, time=datetime.now())"

    def fetchname(self, wavevar, time, region, level_ref=level_ref):
        return fetchname(wavevar, time, region, level_ref)

    def load(self, filepath, plot=False):
        """
        rep_hour = ((self.fetch_datetimestamp.hour % 24) // 6) * 6
        target_date = self.fetch_datetimestamp.replace(minute=0, hour=rep_hour, second=0, microsecond=0)

        grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]
        """
        return fetch_util.loadgrib(filepath, plot)

    def __str__(self):
        info = "RDWPS info goes here"
        args = "(south=-90, north=90, west=180, east=-180, time=datetime.now())"
        return fetch_util.str_def(self, info, args)
