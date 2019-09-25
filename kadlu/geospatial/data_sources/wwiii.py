"""
    Kadlu API for the NOAA WaveWatch III Datastore

    User guides:
        https://github.com/NOAA-EMC/WW3/wiki/WAVEWATCH-III-User-Guide

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
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
import utm 
from shapely.geometry import box


class Boundary():
    def __init__(self, south, north, west, east, fetchvar=''):
        self.west, self.south, self.east, self.north = self.ll_2_utm(south, north, west, east)
        self.fetchvar = fetchvar

    def __str__(self): return self.fetchvar

    def intersects(self, other):
        """
        #separating axis theorem: quickly check for intersecting regions
        #    https://en.wikipedia.org/wiki/Hyperplane_separation_theorem#Use_in_collision_detection
        return not (self.east < other.west or 
                    self.west > other.east or 
                    self.north < other.south or 
                    self.south > other.north)
        """

        return box(self.west, self.south, self.east, self.north).intersects(box(other.west, other.south, other.east, other.north))

    def ll_2_utm(self, south, north, west, east):
        llcrn = utm.from_latlon(south, west)
        urcrn = utm.from_latlon(north, east)
        # returns minx, miny, maxx, maxy -> this is the order used for shapely
        return llcrn[1], llcrn[0], urcrn[1], urcrn[0]


wwiii_global = Boundary(-80, 84, -180, 180, 'glo_30m')  # global
wwiii_regions = [
        Boundary( 15,  47,  -99,  -60, 'at_4m'),  # atlantic, also available : at_10m for 10 minute resolution
        Boundary( 15,  50, -165, -116, 'wc_4m'),  # US west, also available : wc_10m for 10 minute resolution
        Boundary( 48,  74,  140, -120, 'ak_4m'),  # alaska, also available : ak_10m for 10 minute resolution
        Boundary( 65,  84, -180,  180, 'ao_30m'), # arctic ocean
        Boundary(-20,  30,  130, -145, 'ep_10m'), # pacific
        wwiii_global
    ]

def ll_2_regionstr(south, north, west, east):
    """ convert input bounds to region strings using universal transverse mercator projection """
    regions = [str(reg) for reg in wwiii_regions if Boundary(south, north, west, east).intersects(reg)]
    # default to low 30m resolution for global queries
    if (str(wwiii_global) in regions and len(regions) > 4): return [str(wwiii_global)]
    ix = [reg != str(wwiii_global) for reg in regions]

    return np.array(regions)[ix]


def fetchname(wavevar, time, region): return f"multi_1.{region}.{wavevar}.{time.strftime('%Y%m')}.grb2"


def fetch_wwiii(wavevar, south, north, west, east, start, end):
    """
    reg = 'glo_30m'
    regions = [reg]
    wavevar = 'hs'
    start = datetime(2017, 2, 3)
    end = datetime(2017, 2, 4)
    """
    regions = ll_2_regionstr(south, north, west, east)
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
    """
    south=-90
    north=90
    west=-180
    east=180
    reg = 'glo_30m'
    """
    val = np.array([])
    lat = np.array([])
    lon = np.array([])
    timestamps = np.array([])
    regions = ll_2_regionstr(south, north, west, east)
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

