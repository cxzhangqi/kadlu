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
import requests
import shutil
import pygrib

from kadlu.geospatial.data_sources.fetch_util import \
storage_cfg, database_cfg, Boundary, ll_2_regionstr, dt_2_epoch, epoch_2_dt, str_def


conn, db = database_cfg()

wwiii_src = "https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/"
            
wwiii_global = Boundary(-90, 90, -180, 180, 'glo_30m')  # global
wwiii_regions = [
        #region boundaries as defined in WWIII docs:
        #    https://polar.ncep.noaa.gov/waves/implementations.php
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
    #if str(wwiii_global) not in regions: regions = np.append(regions, str(wwiii_global))
    time = datetime(start.year, start.month, 1)
    filenames = []

    while time <= end:
        for reg in regions:
            fname = fetchname(wavevar, time, reg)
            fetchfile = f"{storage_cfg()}{fname}"
            print(f"downloading {fname} from NOAA WaveWatch III...")
            fetchurl = f"{wwiii_src}{time.strftime('%Y/%m')}/{reg}/{fname}"
            with requests.get(fetchurl, stream=True) as payload:
                assert payload.status_code == 200, 'couldn\'t retrieve file'
                with open(fetchfile, 'wb') as f:
                    shutil.copyfileobj(payload.raw, f)
            filenames.append(fetchfile)

        # on this datasource, data is sorted per month
        # some months have more days than exactly 4 weeks
        time += timedelta(weeks=4)
        while (fetchname(wavevar, time, reg) == fname): time += timedelta(days=1)

    for fetchfile in filenames:
        print(f"preparing {fetchfile.split('/')[-1]} for the database...")
        grib = pygrib.open(fetchfile)
        assert(grib.messages > 0)
        val = np.array([])
        lat = np.array([])
        lon = np.array([])
        t   = np.array([])
        nulls = 0

        for msg in grib:
            print(f"processing messages from {msg.validDate}...")
            z, y, x = msg.data()
            val = np.append(val, z[~z.mask].data)
            lat = np.append(lat, y[~z.mask]) 
            lon = np.append(lon, x[~z.mask]-180)
            t = np.append(t, dt_2_epoch([msg.validDate for each in z[~z.mask].data]))
            nulls += sum(sum(z.mask))

        print("formatting and inserting...")
        src = np.array(['wwiii' for each in val])
        grid = list(map(tuple, np.vstack((val, lat, lon, t, src)).T))
        n1 = db.execute(f"SELECT COUNT(*) FROM {wavevar}").fetchall()[0][0]
        db.executemany(f"INSERT OR IGNORE INTO {wavevar} VALUES (?,?,?,?,?)", grid)
        n2 = db.execute(f"SELECT COUNT(*) FROM {wavevar}").fetchall()[0][0]
        db.execute("COMMIT")
        conn.commit()

        print(f"processed and inserted {n2-n1} rows. "
              f"{nulls} null values removed, "
              f"{len(grid) - (n2-n1)} duplicate rows ignored")

    return 

def load_wwiii(wavevar, south, north, west, east, start, end):
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

    return:
        val, lat, lon, time as np arrays
        (time is datetime)
    """

    db.execute(' AND '.join([f"SELECT * FROM {wavevar} WHERE lat >= ?",
                                                            "lat <= ?",
                                                            "lon >= ?",
                                                            "lon <= ?",
                                                           "time >= ?",
                                                           "time <= ?"]),
               tuple(map(str, [south, north, west, east, 
                               dt_2_epoch(start), dt_2_epoch(end)])))
    val, lat, lon, time, source = np.array(db.fetchall(), dtype=object).T

    return val, lat, lon, epoch_2_dt(time)


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
        return str_def(self, info, args)

