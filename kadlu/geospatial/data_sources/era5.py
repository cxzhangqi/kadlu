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
    def fetch_waveheight(self, time=datetime.now()):      
        return fetch_era5('significant_height_of_combined_wind_waves_and_swell', time)
    
    def fetch_wavedirection(self, time=datetime.now()):   
        return fetch_era5('mean_wave_direction', time)

    def fetch_waveperiod(self, time=datetime.now()):      
        return fetch_era5('mean_wave_period', time)

    def load(self, filepath, plot=False): return fetch_util.loadgrib(filepath, plot)

    header = "(time=datetime.now())"

    def printFcns(self): 
        for fcn in [self.fetch_waveheight, self.fetch_wavedirection, self.fetch_waveperiod]:
            print(fcn.__name__ + self.header)

"""
Era5().printFcns()
"""






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

