"""
    Kadlu API for the NOAA WaveWatch III Datastore

    User guides:
        https://github.com/NOAA-EMC/WW3/wiki/WAVEWATCH-III-User-Guide

    Data model description (boundary definitions, map visualizations, etc)
        https://polar.ncep.noaa.gov/waves/implementations.php
"""

import os
import shutil
import logging
import requests
from os.path import isfile
from datetime import datetime, timedelta

import numpy as np
import pygrib

import kadlu.geospatial.data_sources.fetch_handler
from kadlu.geospatial.data_sources.data_util import                 \
        ll_2_regionstr,                                             \
        database_cfg,                                               \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized,                                                 \
        dt_2_epoch,                                                 \
        fmt_coords,                                                 \
        Boundary,                                                   \
        str_def


conn, db = database_cfg()
wwiii_src = "https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/"

# region boundaries as defined in WWIII docs:
#    https://polar.ncep.noaa.gov/waves/implementations.php
wwiii_varmap = dict(zip(
        ('hs','dp','tp', 'windU', 'windV', 'wind'),
        ('waveheight','wavedir','waveperiod', 'wind_u', 'wind_v', 'wind_uv')))

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
            True if new data was fetched, else False
    """
    assert 6 == sum(map(lambda kw: kw in kwargs.keys(),
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'
    t = datetime(kwargs['start'].year, kwargs['start'].month, 1)
    assert kwargs['end'] - kwargs['start'] <= timedelta(days=1), \
            'use fetch_handler for this'

    if serialized(kwargs, f'fetch_wwiii_{wwiii_varmap[var]}'): return False
    #print("WWIII NOTICE: resolution selection not implemented yet. defaulting to 0.5°")
    regions = ['glo_30m']

    assert regions == ['glo_30m'], 'invalid region string'
    reg = regions[0]
    fname = f"multi_1.{reg}.{var}.{t.strftime('%Y%m')}.grb2"
    fetchfile = f"{storage_cfg()}{fname}"

    # if file hasnt been downloaded, fetch it
    if not isfile(fetchfile):# and kwargs['start'].day == 1: 
        if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
        logging.info(f'WWIII {kwargs["start"].date().isoformat()} {var}: '
                     f'downloading {fname} from NOAA WaveWatch III...')
        if reg == 'glo_30m' and not 'wind' in var and t.year >= 2018:
            fetchurl = f"{wwiii_src}{t.strftime('%Y/%m')}/gribs/{fname}"
        else:
            fetchurl = f"{wwiii_src}{t.strftime('%Y/%m')}/{reg}/{fname}"
        with requests.get(fetchurl, stream=True) as payload:
            assert payload.status_code == 200, 'couldn\'t retrieve file'
            with open(fetchfile, 'wb') as f:
                shutil.copyfileobj(payload.raw, f)
        if 'lock' in kwargs.keys(): kwargs['lock'].release()

    # function to insert the parsed data to local database
    def insert(table, agg, null, kwargs):
        if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
        n1 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
        db.executemany(f"INSERT OR IGNORE INTO {table} VALUES (?,?,?,CAST(? AS INT),?)", agg.T)
        n2 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
        db.execute("COMMIT")
        conn.commit()
        insert_hash(kwargs, f'fetch_wwiii_{wwiii_varmap[var]}')
        if 'lock' in kwargs.keys(): kwargs['lock'].release()
        logging.info(f"WWIII {kwargs['start'].date().isoformat()} {table}: "
                f"processed and inserted {n2-n1} rows for region {fmt_coords(kwargs)}. "
                f"{null} null values removed, "
                f"{len(agg[0]) - (n2-n1)} duplicates ignored")

    # open the file, parse data, insert values
    grib = pygrib.open(fetchfile)
    assert grib.messages > 0, f'problem opening {fetchfile}'
    null = 0
    agg = np.array([[],[],[],[],[]])
    grbvar = grib[1]['name']
    table = f'{var}{grbvar[0]}' if var == 'wind' else var
    for msg, num in zip(grib, range(1, grib.messages)):
        if msg['name'] != grbvar:
            insert(table, agg, null, kwargs)
            table = f'{var}{msg["name"][0]}' if var == 'wind' else var
            agg = np.array([[],[],[],[],[]])
            grbvar = msg['name']
            null = 0
        if msg.validDate < kwargs['start']: continue
        if msg.validDate > kwargs['end']:   continue
        z, y, x = msg.data()
        src = np.array(['wwiii' for each in z[~z.mask].data])
        grid = np.vstack((z[~z.mask].data, 
                          y[~z.mask], 
                          ((x[~z.mask] + 180) % 360 ) - 180, 
                          dt_2_epoch([msg.validDate for each in z[~z.mask].data]), 
                          src)).astype(object)
        agg = np.hstack((agg, grid))
        null += sum(sum(z.mask))
    insert(table, agg, null, kwargs)

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
        val, lat, lon, epoch as np arrays of floats
    """
    if 'time' in kwargs.keys() and not 'start' in kwargs.keys():
        kwargs['start'] = kwargs['time']
        del kwargs['time']
    if not 'end' in kwargs.keys(): 
        kwargs['end'] = kwargs['start'] + timedelta(hours=3)

    assert 6 == sum(map(lambda kw: kw in kwargs.keys(),
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'

    kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
            wwiii_varmap[var], 'wwiii', parallel=1, **kwargs)

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
               dt_2_epoch(kwargs['start']), dt_2_epoch(kwargs['end'])]))
       )

    slices = np.array(db.fetchall(), dtype=object).T
    #assert len(slices) == 5, "no data found, try adjusting query bounds or fetching some"
    if len(slices) == 0:
        logging.warning(f'WWIII {var}: no data found in region {fmt_coords(kwargs)}, returning empty arrays')
        return np.array([[],[],[],[]])

    val, lat, lon, epoch, source = slices
    return np.array((val, lat, lon, epoch), dtype=np.float)


class Wwiii():
    """ collection of module functions for fetching and loading """

    def fetch_wavedirection(self,   **kwargs):  return fetch_wwiii('dp',    kwargs)
    def fetch_waveperiod(self,      **kwargs):  return fetch_wwiii('tp',    kwargs)
    def fetch_windwaveheight(self,  **kwargs):  return fetch_wwiii('hs',    kwargs)
    def fetch_wind_u(self,          **kwargs):  return fetch_wwiii('wind',  kwargs)
    def fetch_wind_v(self,          **kwargs):  return fetch_wwiii('wind',  kwargs)
    def fetch_wind_uv(self,         **kwargs):  return fetch_wwiii('wind',  kwargs)

    def load_wavedirection(self,    **kwargs):  return load_wwiii('dp',     kwargs)
    def load_waveperiod(self,       **kwargs):  return load_wwiii('tp',     kwargs)
    def load_windwaveheight(self,   **kwargs):  return load_wwiii('hs',     kwargs)
    def load_wind_u(self,           **kwargs):  return load_wwiii('windU',  kwargs)
    def load_wind_v(self,           **kwargs):  return load_wwiii('windV',  kwargs)
    def load_wind_uv(self,          **kwargs):

        kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
                wwiii_varmap['windU'], 'wwiii', parallel=1, **kwargs)
        kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
                wwiii_varmap['windV'], 'wwiii', parallel=1, **kwargs)

        sql = ' AND '.join(['SELECT windU.val, windU.lat, windU.lon, windU.time, windV.val FROM windU '\
            'INNER JOIN windV '\
            'ON windU.lat == windV.lat',
            'windU.lon == windV.lon',
            'windU.time == windV.time '\
            'WHERE windU.lat >= ?',
            'windU.lat <= ?',
            'windU.lon >= ?',
            'windU.lon <= ?',
            'windU.time >= ?',
            'windU.time <= ?']) + ' ORDER BY windU.time, windU.lat, windU.lon ASC'

        db.execute(sql, tuple(map(str, [
                kwargs['south'],                kwargs['north'], 
                kwargs['west'],                 kwargs['east'], 
                dt_2_epoch(kwargs['start']),    dt_2_epoch(kwargs['end'])
            ])))
        qry = np.array(db.fetchall()).T
        #assert len(qry) > 0, \
        #        f'no windspeed data found in region {fmt_coords(kwargs)}. consider expanding the region'
        if len(qry) == 0:
            logging.warning(f'ERA5 wind_uv: no data found in region {fmt_coords(kwargs)}, returning empty arrays')
            return np.array([[],[],[],[],[]])
        wind_u, lat, lon, epoch, wind_v = qry
        val = np.sqrt(np.square(wind_u.astype(float)), np.square(wind_v.astype(float)))
        return np.array((val, lat, lon, epoch)).astype(float)

    def __str__(self):
        info = '\n'.join(["WAVEWATCH III: a third generation wave height,",
            "water depth, and current hindcasting modeling framework.",
            "Developed for the community at NOAA/NCEP. About WWIII:",
            "\thttps://github.com/NOAA-EMC/WW3/wiki/About-WW3",
            "Model descriptions:",
            "\thttps://polar.ncep.noaa.gov/waves/implementations.php"])
        args = "(south, north, west, east, start, end)"
        return str_def(self, info, args)

