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
import warnings


def fetchname(wavevar, time):
    return f"ERA5_reanalysis_{wavevar}_{time.strftime('%Y-%m-%d_%Hh')}.grb2"


def fetch_era5(wavevar, start, end):
    time = datetime(start.year, start.month, start.day, start.hour)
    fetchfiles = []

    while time <= end:
        fname= f"{storage_cfg()}{fetchname(wavevar, time)}"
        print(f"Downloading {fetchname(wavevar, time)} from Copernicus Climate Data Store...")
        c = cdsapi.Client()
        try:
            c.retrieve('reanalysis-era5-single-levels', {
                'product_type'  : 'reanalysis',
                'format'        : 'grib',
                'variable'      : wavevar,
                'year'          : time.strftime("%Y"),
                'month'         : time.strftime("%m"),
                'day'           : time.strftime("%d"),
                'time'          : time.strftime("%H:00")
            }, fname)
        except Exception:
            warnings.warn(f"No data for {time.ctime()}")
            time += timedelta(hours=1)
            continue

        fetchfiles.append(fname)
        time += timedelta(hours=1)

    return fetchfiles


def load_era5(wavevar, start, end, south, north, west, east, plot):
    fetchfiles = []  # used for plotting, this may be removed later
    val = np.array([])
    lat = np.array([])
    lon = np.array([])
    timestamps = np.array([])

    # adjust temporal resolution to hourly
    time = datetime(start.year, start.month, start.day, start.hour)
    assert(time <= end)
    print(f"Loading {wavevar} from {int((end-start).total_seconds() /60 /60)+1} hourly global data files")

    while time <= end:
        # get the filename
        fname = f"{storage_cfg()}{fetchname(wavevar, time)}" 
        fetchfiles.append(fname)  # used for plotting
        if not os.path.isfile(fname): fetch_era5(wavevar, time, time)

        # open the file and get the raw data
        grib = pygrib.open(fname)
        assert(grib.messages == 1)
        msg = grib[1]
        z, y, x = msg.data()
        x -= 180  # normalize longitudes

        # build index to collect points in area of interest
        latix = np.array([l >= south and l <= north for l in y[~z.mask]])
        lonix = np.array([l >= west and l <= east for l in x[~z.mask]])
        ix = latix & lonix

        # append points within AoI to return arrays
        val = np.append(val, z[~z.mask][ix])
        lat = np.append(lat, y[~z.mask][ix])
        lon = np.append(lon, x[~z.mask][ix])
        timestamps = np.append(timestamps, [time for x in range(sum(ix))])

        time += timedelta(hours=1)

    if plot is not False: fetch_util.plot_sample_grib(fetchfiles, plot)
    return val.data, lat, lon, timestamps


class Era5():
    def fetch_windwaveswellheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_era5('significant_height_of_combined_wind_waves_and_swell', start, end)
    def fetch_wavedirection(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_era5('mean_wave_direction', start, end)
    def fetch_waveperiod(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_era5('mean_wave_period', start, end)

    def load_windwaveswellheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now(), plot=False):
        return load_era5('significant_height_of_combined_wind_waves_and_swell', start, end, south, north, west, east, plot)
    def load_wavedirection(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now(), plot=False):
        return load_era5('mean_wave_direction', start, end, south, north, west, east, plot)
    def load_waveperiod(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now(), plot=False):
        return load_era5('mean_wave_period', start, end, south, north, west, east, plot)

    def __str__(self):
        info = '\n'.join([
                "Era5 Global Dataset from Copernicus Climate Datastore",
                "\thttps://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview"
            ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(), end=datetime())"
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

