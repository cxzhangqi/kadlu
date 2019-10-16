"""
    Kadlu API for the NOAA WaveWatch III Datastore

    User guides:
        https://github.com/NOAA-EMC/WW3/wiki/WAVEWATCH-III-User-Guide

    Data model description (boundary definitions, map visualizations, etc)
        https://polar.ncep.noaa.gov/waves/implementations.php

    Oliver Kirsebom
    Casey Hilliard
    Matthew Smith
"""

import numpy as np
from datetime import datetime, timedelta
import os
from urllib.request import urlretrieve
import pygrib
from kadlu.geospatial.data_sources import fetch_util
from kadlu.geospatial.data_sources.fetch_util import storage_cfg, Boundary, ll_2_regionstr


wwiii_global = Boundary(-90, 90, -180, 180, 'glo_30m')  # global
wwiii_regions = [
        """ region boundaries as defined in WWIII docs:
            https://polar.ncep.noaa.gov/waves/implementations.php
        """
        Boundary( 15,  47,  -99,  -60, 'at_4m'),    # atlantic
        Boundary( 15,  50, -165, -116, 'wc_4m'),    # US west
        Boundary( 48,  74,  140,  180, 'ak_4m'),    # alaska
        Boundary( 48,  74, -180, -120, 'ak_4m'),    # alaska
        Boundary( 65,  84, -180,  180, 'ao_30m'),   # arctic ocean
        Boundary(-20,  30,  130,  180, 'ep_10m'),   # pacific
        Boundary(-20,  30, -180, -145, 'ep_10m')    # pacific
    ]


def fetchname(wavevar, time, region):
    """ generate filename for given wave variable, time, and region """
    return f"multi_1.{region}.{wavevar}.{time.strftime('%Y%m')}.grb2"


def fetch_wwiii(wavevar, south, north, west, east, start, end):
    """ download wwiii data and return associated filepaths

    args:
        wavevar: string
            the variable short name of desired wave parameter according to WWIII docs
            the complete list of variable short names can be found here (under 'model output')
            https://polar.ncep.noaa.gov/waves/implementations.php
        south, north: float
            ymin, ymax coordinate boundaries (latitude). range: -90, 90
        west, east: float
            xmin, xmax coordinate boundaries (longitude). range: -180, 180
        start: datetime
            the start of the desired time range
        end: datetime
            the end of the desired time range

    return:
        filenames: list
            list of strings containing complete file paths of fetched data
    """
    regions = ll_2_regionstr(south, north, west, east, wwiii_regions, [str(wwiii_global)])
    time = datetime(start.year, start.month, 1)
    filenames = []

    while time <= end:
        for reg in regions:
            fname = fetchname(wavevar, time, reg)
            fetchfile = f"{storage_cfg()}{fname}"
            print(f"Downloading {fname} from NOAA WaveWatch III...")
            fetchurl = f"https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/{time.strftime('%Y/%m')}/gribs/{fname}"
            urlretrieve(fetchurl, fetchfile)
            filenames.append(fetchfile)

        # on this datasource, data is sorted per month
        # some months have more days than exactly 4 weeks
        time += timedelta(weeks=4)
        while (fetchname(wavevar, time, reg) == fname): time += timedelta(days=1)

    return filenames

def load_wwiii(wavevar, south, north, west, east, start, end, plot=False):
    """ return downloaded wwiii data for specified wavevar according to given time, lat, lon boundaries

    args:
        wavevar: string
            the variable short name of desired wave parameter according to WWIII docs
            the complete list of variable short names can be found here (under 'model output')
            https://polar.ncep.noaa.gov/waves/implementations.php
        south, north: float
            ymin, ymax coordinate boundaries (latitude). range: -90, 90
        west, east: float
            xmin, xmax coordinate boundaries (longitude). range: -180, 180
        start: datetime
            the start of the desired time range
        end: datetime
            the end of the desired time range
        plot: boolean
            if true a plot will be output (experimental feature)

    return:
        filenames: list
            list of strings containing complete file paths of fetched data
    """
    val = np.array([])
    lat = np.array([])
    lon = np.array([])
    timestamps = np.array([])
    regions = ll_2_regionstr(south, north, west, east, wwiii_regions, [str(wwiii_global)])
    time = datetime(start.year, start.month, 1)

    while time <= end:
        for reg in regions:
            fname = fetchname(wavevar, time, reg)
            fetchfile = f"{storage_cfg()}{fname}"
            if not os.path.isfile(fetchfile): fetch_wwiii(wavevar, south, north, west, east, start, end)

            grib = pygrib.open(fetchfile)
            for msg in grib:
                msgtime = msg.validDate
                if msgtime < start : continue
                if msgtime > end : break
                z, y, x = msg.data()
                x -= 360  # normalize longitude
                x[x < -180] += 360

                for slx in range(z.shape[0]):
                    latix = np.array([l >= south and l <= north for l in y[slx]])
                    lonix = np.array([l >= west and l <= east for l in x[slx]])
                    ix = latix & lonix

                    val = np.append(val, z[slx][ix])
                    lat = np.append(lat, y[slx][ix])
                    lon = np.append(lon, x[slx][ix])
                    timestamps = np.append(timestamps, [msgtime for x in range(sum(ix))])

        # on this datasource, data is sorted per month
        # some months have more days than exactly 4 weeks
        time += timedelta(weeks=4)
        while (fetchname(wavevar, time, reg) == fname): time += timedelta(days=1)

    return val, lat, lon, timestamps


class Wwiii():
    """ collection of module functions for fetching and loading. abstracted to include a seperate function for each variable """

    def fetch_windwaveheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_wwiii('hs', south, north, west, east, start, end)
    def fetch_wavedirection(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_wwiii('dp', south, north, west, east, start, end)
    def fetch_waveperiod(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_wwiii('tp', south, north, west, east, start, end)

    def load_windwaveheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_wwiii('hs', south, north, west, east, start, end)
    def load_wavedirection(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_wwiii('dp', south, north, west, east, start, end)
    def load_waveperiod(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_wwiii('tp', south, north, west, east, start, end)

    def __str__(self):
        info = '\n'.join([ "Wavewatch info goes here" ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(), end=datetime())"
        return fetch_util.str_def(self, info, args)

