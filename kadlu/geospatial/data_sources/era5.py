"""
    API for Era5 dataset from Copernicus Climate Datastore

    Metadata regarding the dataset can be found here:
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

    matt_s 2019-08
"""
import numpy as np
from datetime import datetime, timedelta
import os
import cdsapi
import pygrib
from kadlu.geospatial.data_sources import fetch_util 
from kadlu.geospatial.data_sources.fetch_util import storage_cfg 

def fetchname(wavevar, time):
    return f"ERA5_reanalysis_{wavevar}_{time.strftime('%Y-%m-%d_%Hh')}.grb2"


def fetch_era5(wavevar, time):
    fetchfile = f"{storage_cfg()}{fetchname(wavevar, time)}"
    if os.path.isfile(fetchfile):
        print(f"File {fetchname(wavevar, time)} exists, skipping retrieval...")
    else:
        print(f"Downloading {fetchname(wavevar, time)} from Copernicus Climate Data Store...")
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


class Era5():
    def fetch_windwaveswellheight(self, time=datetime.now()):      
        return fetch_era5('significant_height_of_combined_wind_waves_and_swell', time)
    
    def fetch_wavedirection(self, time=datetime.now()):   
        return fetch_era5('mean_wave_direction', time)

    def fetch_waveperiod(self, time=datetime.now()):      
        return fetch_era5('mean_wave_period', time)

    def load(self, filepath, plot=False): 
        """
        grbs=pygrib.open(grib)
        grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]
        """
        return fetch_util.loadgrib(filepath, plot)

    header = "(time=datetime.now())"

    def print_fcns(self): 
        for fcn in [self.fetch_waveheight, self.fetch_wavedirection, self.fetch_waveperiod]:
            print(fcn.__name__ + self.header)

