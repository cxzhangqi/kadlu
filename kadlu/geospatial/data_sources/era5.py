"""
    API for Era5 dataset from Copernicus Climate Datastore

    Metadata regarding the dataset can be found here:
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

    Oliver Kirsebom
    Casey Hilliard
    Matthew Smith
"""

import cdsapi
import os
from datetime import datetime, timedelta
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
    return [fetchfile]


class Era5():
    def fetch_windwaveswellheight(self, time=datetime.now()):      
        return fetch_era5('significant_height_of_combined_wind_waves_and_swell', time)
    def fetch_wavedirection(self, time=datetime.now()):   
        return fetch_era5('mean_wave_direction', time)
    def fetch_waveperiod(self, time=datetime.now()):      
        return fetch_era5('mean_wave_period', time)
    def load_windwaveswellheight(self, time=datetime.now(), plot=False):
        #grbs=pygrib.open(grib)
        #grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]
        return fetch_util.loadgrib(self.fetch_windwaveswellheight(time), plot)
    def load_wavedirection(self, time=datetime.now()):
        return fetch_util.loadgrib(self.fetch_wavedirection(time), plot)
    def load_waveperiod(self, time=datetime.now()):
        return fetch_util.loadgrib(self.fetch_waveperiod(time), plot)

    def __str__(self):
        info = "Era5 Global Dataset from Copernicus Climate Datastore"
        args = "(time=datetime(), plot=False)"
        return fetch_util.str_def(self, info, args)


"""
print(Era5())

time = datetime(2018, 1, 1)
fnames = Era5().fetch_windwaveswellheight(time)
wave, lat, lon = Era5().load_windwaveswellheight(time)

#wave, lat, lon = Era5().load_windwaveswellheight(time, "wind, wave, swell height")
wave, lat, lon = Era5().load_windwaveswellheight(time, plot=False)
"""

