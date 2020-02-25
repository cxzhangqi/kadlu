"""
    Kadlu API for the NOAA WaveWatch III Datastore

    User guides:
        https://github.com/NOAA-EMC/WW3/wiki/WAVEWATCH-III-User-Guide

    Data model description (boundary definitions, map visualizations, etc)
        https://polar.ncep.noaa.gov/waves/implementations.php
"""

import numpy as np
from datetime import datetime, timedelta
import os
from os.path import isfile
import requests
import shutil
import pygrib
import warnings

from kadlu.geospatial.data_sources.data_util import                 \
        ll_2_regionstr,                                             \
        database_cfg,                                               \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized,                                                 \
        dt_2_epoch,                                                 \
        epoch_2_dt,                                                 \
        Boundary,                                                   \
        str_def


conn, db = database_cfg()
wwiii_src = "https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/"

# region boundaries as defined in WWIII docs:
#    https://polar.ncep.noaa.gov/waves/implementations.php
wwiii_varmap = dict(zip(
    ('hs','dp','tp', 'wind'),
    ('waveheight','wavedirection','waveperiod', 'windspeed')))

wwiii_global = Boundary(-90, 90, -180, 180, 'glo_30m')  # global
wwiii_regions = [
        Boundary( 15,  47,  -99,  -60, 'at_4m'),    # atlantic
        Boundary( 15,  50, -165, -116, 'wc_4m'),    # US west
        Boundary( 48,  74,  140,  180, 'ak_4m'),    # alaska (west)
        Boundary( 48,  74, -180, -120, 'ak_4m'),    # alaska (east)
        Boundary( 65,  84, -180,  180, 'ao_30m'),   # arctic ocean
        Boundary(-20,  30,  130,  180, 'ep_10m'),   # pacific (west)
        Boundary(-20,  30, -180, -145, 'ep_10m')]   # pacific (east) 



def fetch_wwiii(var, kwargs):
    """ download wwiii data and return associated filepaths

        args:
            var: string
                the variable name of desired parameter according to WWIII docs
                the complete list of variables can be found at the following 
                URL under 'model output'
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
            nothing. some status messages are printed to stdout
    """
    assert 6 == sum(map(lambda kw: kw in kwargs.keys(),
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'
    #assert kwargs['start'] < kwargs['end']
    t = datetime(kwargs['start'].year, kwargs['start'].month, 1)
    assert t.month == (kwargs['end']-timedelta(hours=1)).month, \
            'use fetch_handler for this'


    def fetchname(var, t, region):
        """ generate filename for given wave variable, time, and region """
        return f"multi_1.{region}.{var}.{t.strftime('%Y%m')}.grb2"


    # used for sanity checking the fetch_handler hash index
    # this can be removed later for increased performance
    """
    if 'lock' in kwargs.keys(): 
        kwargs['lock'].acquire()
        if serialized(kwargs, f'fetch_wwiii_{wwiii_varmap[var]}'): 
            print(f'WWIII DEBUG MSG: fetch_handler index check failed')
            return False
        kwargs['lock'].release()
    """

    warnings.warn("resolution selection not implemented yet. defaulting to 0.5°")
    regions = ['glo_30m']

    assert regions == ['glo_30m'], 'invalid region string'
    reg = regions[0]
    fname = fetchname(var, t, reg)
    fetchfile = f"{storage_cfg()}{fname}"
    if not isfile(fetchfile) and kwargs['start'].day != 1:
        while not isfile(fetchfile): pass  # hang thread until done downloading
        if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
        if 'lock' in kwargs.keys(): kwargs['lock'].release()
    elif not isfile(fetchfile) and kwargs['start'].day == 1: 
        if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
        print(f'WWIII {kwargs["start"].date().isoformat()} {var}: '
              f'downloading {fname} from NOAA WaveWatch III...')
        if reg == 'glo_30m' and t.year >= 2018:
            fetchurl = f"{wwiii_src}{t.strftime('%Y/%m')}/gribs/{fname}"
        else:
            fetchurl = f"{wwiii_src}{t.strftime('%Y/%m')}/{reg}/{fname}"
        with requests.get(fetchurl, stream=True) as payload:
            assert payload.status_code == 200, 'couldn\'t retrieve file'
            with open(fetchfile, 'wb') as f:
                shutil.copyfileobj(payload.raw, f)
        if 'lock' in kwargs.keys(): kwargs['lock'].release()

    grib = pygrib.open(fetchfile)
    assert grib.messages > 0, f'problem opening {fetchfile}'

    def insert(table, agg, null, kwargs):
        if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
        n1 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
        db.executemany(f"INSERT OR IGNORE INTO {table} VALUES (?,?,?,CAST(? AS INT),?)", agg.T)
        n2 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
        db.execute("COMMIT")
        conn.commit()
        if 'lock' in kwargs.keys(): kwargs['lock'].release()
        print(f"WWIII {kwargs['start'].date().isoformat()} {table}: "
              f"processed and inserted {n2-n1} rows. "
              f"{null} null values removed, "
              f"{len(agg[0]) - (n2-n1)} duplicates ignored")
    
    null = 0
    agg = np.array([[],[],[],[],[]])
    grbvar = grib[1]['name']
    for msg, num in zip(grib, range(1, grib.messages)):
        if msg['name'] != grbvar:
            table = f'{var}{msg["name"][0]}' if var == 'wind' else var
            insert(table, agg, null, kwargs)
            agg = np.array([[],[],[],[],[]])
            grbvar = msg['name']
            null = 0
        if msg.validDate < kwargs['start']: continue
        if msg.validDate > kwargs['end']:   continue
        z, y, x = msg.data()
        src = np.array(['wwiii' for each in z[~z.mask].data])
        grid = np.vstack((
                z[~z.mask].data, 
                y[~z.mask], 
                ((x[~z.mask] + 180) % 360 ) - 180, 
                dt_2_epoch([msg.validDate for each in z[~z.mask].data]), 
                src
            )).astype(object)
        agg = np.hstack((agg, grid))
        null += sum(sum(z.mask))
        
    table = f'{var}{msg["name"][0]}' if var == 'wind' else var
    insert(table, agg, null, kwargs)


    insert_hash(kwargs, f'fetch_wwiii_{wwiii_varmap[var]}')
    return True


def load_wwiii(var, kwargs):
    """ return downloaded wwiii data for specified wavevar according to given time, lat, lon boundaries

    args:
        var: string
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
    assert not 'time' in kwargs.keys(), 'nearest time search not implemented yet'
    assert 6 == sum(map(lambda kw: kw in kwargs.keys(),
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'

    db.execute(' AND '.join([
           f'SELECT * FROM {var} WHERE lat >= ?',
            'lat <= ?',
            'lon >= ?',
            'lon <= ?',
            'time >= ?',
            'time <= ? ']) + ' ORDER BY time, lat, lon ASC',
           tuple(map(str, [
               kwargs['south'], kwargs['north'],
               kwargs['west'],  kwargs['east'], 
               dt_2_epoch(kwargs['start'])[0], dt_2_epoch(kwargs['end'])[0]]))
       )

    slices = np.array(db.fetchall(), dtype=object).T
    assert len(slices) == 5, \
            "no data found, try adjusting query bounds or fetching some"
    val, lat, lon, time, source = slices

    return np.array((val, lat, lon, time), dtype=np.float)


class Wwiii():
    """ collection of module functions for fetching and loading. abstracted to include a seperate function for each variable """

    def fetch_windwaveheight(self,  **kwargs):  return fetch_wwiii('hs',    kwargs)
    def fetch_wavedirection(self,   **kwargs):  return fetch_wwiii('dp',    kwargs)
    def fetch_waveperiod(self,      **kwargs):  return fetch_wwiii('tp',    kwargs)
    def fetch_wind_u(self,          **kwargs):  return fetch_wwiii('wind',  kwargs)
    def fetch_wind_v(self,          **kwargs):  return fetch_wwiii('wind',  kwargs)
    def fetch_wind(self,            **kwargs):  return fetch_wwiii('wind',  kwargs)

    def load_windwaveheight(self,   **kwargs):  return load_wwiii('hs',     kwargs)
    def load_wavedirection(self,    **kwargs):  return load_wwiii('dp',     kwargs)
    def load_waveperiod(self,       **kwargs):  return load_wwiii('tp',     kwargs)
    def load_wind_u(self,           **kwargs):  return load_wwiii('windU',  kwargs)
    def load_wind_v(self,           **kwargs):  return load_wwiii('windV',  kwargs)
    def load_wind(self,             **kwargs):
        wind_u = load_wwiii('windU',  kwargs)
        wind_v = load_wwiii('windV',  kwargs)
        wind_uv = wind_u.copy()
        #wind_uv[0] = tuple(zip(wind_u[0], wind_v[0]))
        wind_uv[0] = np.sqrt(np.square(wind_u[0]), np.square(wind_v[0]))
        return wind_uv

    def __str__(self):
        info = '\n'.join([ "Wavewatch info goes here" ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(), end=datetime())"
        return str_def(self, info, args)

