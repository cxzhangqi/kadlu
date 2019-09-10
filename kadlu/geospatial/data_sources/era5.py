"""
    API for Era5 dataset from Copernicus Climate Datastore

    Metadata regarding the dataset can be found here:
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

    Oliver Kirsebom
    Casey Hilliard
    Matthew Smith
"""

import cdsapi
import numpy as np
import pygrib
import os
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import fetch_util 
from kadlu.geospatial.data_sources.fetch_util import storage_cfg 


def fetchname(wavevar, time):
    return f"ERA5_reanalysis_{wavevar}_{time.strftime('%Y-%m-%d_%Hh')}.grb2"


def fetch_era5(wavevar, time):
    fetchfile = f"{storage_cfg()}{fetchname(wavevar, time)}"
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


def load_era5(wavevar, time, plot):
    val = np.array([])
    lat = np.array([])
    lon = np.array([])
    fetchfiles = [f"{storage_cfg()}{fetchname(wavevar, time)}"]

    for fname in fetchfiles:
        if not os.path.isfile(fname): fetch_era5(wavevar, time)

        grib = pygrib.open(fname)
        for msg in grib:
            z, y, x = msg.data()
            val = np.append(val, z[~z.mask])
            lat = np.append(lat, y[~z.mask])
            lon = np.append(lon, x[~z.mask])

    if plot is not False: fetch_util.plot_sample_grib(fetchfiles, plot)
    return val, lat, (lon - 180)


class Era5():
    def fetch_windwaveswellheight(self, time=datetime.now()):      
        return fetch_era5('significant_height_of_combined_wind_waves_and_swell', time)
    def fetch_wavedirection(self, time=datetime.now()):   
        return fetch_era5('mean_wave_direction', time)
    def fetch_waveperiod(self, time=datetime.now()):      
        return fetch_era5('mean_wave_period', time)

    def load_windwaveswellheight(self, time=datetime.now(), plot=False):
        return load_era5('significant_height_of_combined_wind_waves_and_swell', time, plot)
    def load_wavedirection(self, time=datetime.now(), plot=False):
        return load_era5('mean_wave_direction', time, plot)
    def load_waveperiod(self, time=datetime.now(), plot=False):
        return load_era5('mean_wave_period', time, plot)

    def __str__(self):
        info = "Era5 Global Dataset from Copernicus Climate Datastore"
        args = "(time=datetime())"
        return fetch_util.str_def(self, info, args)


"""
print(Era5())

wavevar = 'significant_height_of_combined_wind_waves_and_swell'
time = datetime(2018, 1, 1)
fnames = Era5().fetch_windwaveswellheight(time)
wave, lat, lon = Era5().load_windwaveswellheight(time)

#wave, lat, lon = Era5().load_windwaveswellheight(time, "wind, wave, swell height")
wave, lat, lon = Era5().load_windwaveswellheight(time, plot=False)
"""

