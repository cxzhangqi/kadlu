"""
    API for Era5 dataset from Copernicus Climate Datastore
     
    Metadata regarding the dataset can be found here:
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
"""

import os
import logging
import warnings
import configparser
from os.path import isfile, dirname
from datetime import datetime, timedelta

import cdsapi
import pygrib
import numpy as np

import kadlu.geospatial.data_sources.fetch_handler
from kadlu.geospatial.data_sources.data_util    import              \
        database_cfg,                                               \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized,                                                 \
        dt_2_epoch,                                                 \
        epoch_2_dt,                                                 \
        fmt_coords,                                                 \
        dev_null,                                                   \
        str_def,                                                    \
        cfg


logging.getLogger('cdsapi').setLevel(logging.WARNING)

conn, db = database_cfg()

era5_varmap = dict(zip(
        ('significant_height_of_combined_wind_waves_and_swell',
         'mean_wave_direction',
         'mean_wave_period',
         '10m_u_component_of_wind',
         '10m_v_component_of_wind',
         ),
        ('waveheight', 'wavedir', 'waveperiod', 'wind_u', 'wind_v')))


cfg = configparser.ConfigParser()       # read .ini into dictionary object
cfgfile = os.path.join(dirname(dirname(dirname(dirname(__file__)))), "config.ini")
cfg.read(cfgfile)

def era5_cfg(key=None, url=None):
    if 'cdsapi' not in cfg.sections():
        cfg.add_section('cdsapi')

    if key is not None:
        cfg.set('cdsapi', 'key', key)
        with open(cfgfile, 'w') as f:
            cfg.write(f)
    else:
        cfg.set('cdsapi', 'key', '20822:2d1c1841-7d27-4f72-bb8a-9680a073b4c3')
        with open(cfgfile, 'w') as f:
            cfg.write(f)

    if url is not None:
        cfg.set('cdsapi', 'url', url)
        with open(cfgfile, 'w') as f:
            cfg.write(f)
    else:
        cfg.set('cdsapi', 'url', 'https://cds.climate.copernicus.eu/api/v2')
        with open(cfgfile, 'w') as f:
            cfg.write(f)

    return 


def fetch_era5(var, kwargs):
    """ fetch global era5 data for specified variable and time range

        args:
            var: string
                the variable short name of desired wave parameter 
                according to ERA5 docs.  the complete list can be found 
                here (table 7 for wave params):
                https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation#ERA5datadocumentation-Temporalfrequency
            kwargs: dict
                keyword arguments passed from the Era5() class as a dictionary

        return:
            True if new data was fetched, else False 
    """
    # cleaner stack trace by raising outside of try/except
    err = False
    try: c = cdsapi.Client(url=cfg['cdsapi']['url'], key=cfg['cdsapi']['key'])
    except KeyError:
        try: c = cdsapi.Client()
        except Exception:
            err = True  

    if err: 
        raise KeyError('CDS API has not been configured for the ERA5 module. '
                       'obtain an API token from the following URL and run '
                       'kadlu.era5_cfg(url="URL_HERE", key="TOKEN_HERE"). '
                       'https://cds.climate.copernicus.eu/api-how-to')

    assert 6 == sum(map(lambda kw: kw in kwargs.keys(), 
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'
    t = datetime(kwargs['start'].year,   kwargs['start'].month,
                 kwargs['start'].day,    kwargs['start'].hour)
    assert kwargs['end'] - kwargs['start'] <= timedelta(days=1, hours=1), \
            'use fetch_handler for this instead'
        
    # check if data has been fetched already
    if serialized(kwargs, f'fetch_era5_{era5_varmap[var]}'): return False

    # fetch the data
    fname = f'ERA5_reanalysis_{var}_{t.strftime("%Y-%m-%d")}.grb2'
    fpath = f'{storage_cfg()}{fname}'
    if not isfile(fpath):
        with dev_null():
            c.retrieve('reanalysis-era5-single-levels', {
                       'product_type' : 'reanalysis',
                       'format'       : 'grib',
                       'variable'     : var,
                       'year'         : t.strftime("%Y"),
                       'month'        : t.strftime("%m"),
                       'day'          : t.strftime("%d"),
                       'time'         : [datetime(t.year, t.month, t.day, h)
                                         .strftime('%H:00') for h in range(24)]
                    }, fpath)
    
    # load the data file and insert it into the database
    assert isfile(fpath)
    grb = pygrib.open(fpath)
    agg = np.array([[],[],[],[],[]])
    table = var[4:] if var[0:4] == '10m_' else var

    for msg, num in zip(grb, range(1, grb.messages)):
        if msg.validDate < kwargs['start'] or msg.validDate > kwargs['end']: 
            continue

        # read grib data
        z, y, x = msg.data()
        if np.ma.is_masked(z):
            z2 = z[~z.mask].data
            y2 = y[~z.mask]
            x2 = x[~z.mask]
        else:  # wind data has no mask
            z2 = z.reshape(-1)
            y2 = y.reshape(-1)
            x2 = x.reshape(-1)

        # adjust latitude-zero to 180th meridian
        x3 = ((x2 + 180) % 360) - 180

        # index coordinates, select query range subset, aggregate results
        xix = np.logical_and(x3>=kwargs['west'],  x3<=kwargs['east'])
        yix = np.logical_and(y2>=kwargs['south'], y2<=kwargs['north'])
        idx = np.logical_and(xix, yix)
        agg = np.hstack((agg, [z2[idx],
                               y2[idx],
                               x3[idx],
                               dt_2_epoch([msg.validDate for i in z2[idx]]),
                               ['era5' for i in z2[idx]]]))

    # perform the insertion
    if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
    n1 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
    db.executemany(f"INSERT OR IGNORE INTO {table} "
                   f"VALUES (?,?,?,CAST(? AS INT),?)", agg.T)
    n2 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
    db.execute("COMMIT")
    conn.commit()
    insert_hash(kwargs, f'fetch_era5_{era5_varmap[var]}')
    if 'lock' in kwargs.keys(): kwargs['lock'].release()

    logging.info(f"ERA5 {msg.validDate.date().isoformat()} {var}: "
                 f"processed and inserted {n2-n1} rows in region {fmt_coords(kwargs)}. "
                 f"{len(agg[0])- (n2-n1)} duplicates ignored")

    return True


def load_era5(var, kwargs):
    """ load era5 data from local database

        args:
            var:
                variable to be fetched (string)
            kwargs:
                dictionary containing the keyword arguments used for the
                fetch request. must contain south, north, west, east
                keys as float values, and start, end keys as datetimes

        return:
            values:
                values of the fetched var
            lat:
                y grid coordinates
            lon:
                x grid coordinates
            epoch:
                timestamps in epoch hours since jan 1 2000
    """
    if 'time' in kwargs.keys() and not 'start' in kwargs.keys():
        kwargs['start'] = kwargs['time']
        del kwargs['time']
    if not 'end' in kwargs.keys(): 
        kwargs['end'] = kwargs['start'] + timedelta(hours=3)

    assert 6 == sum(map(lambda kw: kw in kwargs.keys(),
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'

    # check for missing data
    kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
            era5_varmap[var], 'era5', parallel=1, **kwargs)

    # load the data
    table = var[4:] if var[0:4] == '10m_' else var  # table cant start with int
    sql = ' AND '.join([f"SELECT * FROM {table} WHERE lat >= ?",
        'lat <= ?',
        'lon >= ?',
        'lon <= ?',
        'time >= ?',
        'time <= ?']) + ' ORDER BY time, lat, lon ASC'
    db.execute(sql, tuple(map(str, [
            kwargs['south'],                kwargs['north'], 
            kwargs['west'],                 kwargs['east'], 
            dt_2_epoch(kwargs['start']), dt_2_epoch(kwargs['end'])
        ])))
    rowdata = np.array(db.fetchall(), dtype=object).T
    #assert len(rowdata) > 0, "no data found for query"
    if len(rowdata) == 0:
        logging.warning(f'ERA5 {var}: no data found in region {fmt_coords(kwargs)}, returning empty arrays')
        return np.array([[],[],[],[]])

    val, lat, lon, epoch, source = rowdata 
    return np.array((val, lat, lon, epoch), dtype=np.float)


class Era5():
    """ collection of module functions for fetching and loading  """
    def fetch_windwaveswellheight(self, **kwargs):
        return fetch_era5('significant_height_of_combined_wind_waves_and_swell', kwargs)
    def fetch_wavedirection(self, **kwargs):
        return fetch_era5('mean_wave_direction', kwargs)
    def fetch_waveperiod(self, **kwargs):
        return fetch_era5('mean_wave_period', kwargs)
    def fetch_wind_u(self, **kwargs):
        return fetch_era5('10m_u_component_of_wind', kwargs)
    def fetch_wind_v(self, **kwargs):
        return fetch_era5('10m_v_component_of_wind', kwargs)
    def fetch_wind_uv(self, **kwargs):
        return fetch_era5('10m_u_component_of_wind', kwargs) and\
               fetch_era5('10m_v_component_of_wind', kwargs)

    def load_windwaveswellheight(self, **kwargs):
        return load_era5('significant_height_of_combined_wind_waves_and_swell', kwargs)
    def load_wavedirection(self, **kwargs):
        return load_era5('mean_wave_direction', kwargs)
    def load_waveperiod(self, **kwargs):
        return load_era5('mean_wave_period', kwargs)
    def load_wind_u(self, **kwargs):
        return load_era5('10m_u_component_of_wind', kwargs)
    def load_wind_v(self, **kwargs):
        return load_era5('10m_v_component_of_wind', kwargs)
    def load_wind_uv(self, **kwargs):
        """ an SQL join is used for loading wind data to deal with the
            edge case where u,v coordinates are not in matching pairs in
            the data. the JOIN ensures only values with a matching pair
            are loaded
        """
        kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
                era5_varmap['10m_u_component_of_wind'], 'era5', parallel=1, **kwargs)
        kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
                era5_varmap['10m_v_component_of_wind'], 'era5', parallel=1, **kwargs)

        sql = ' AND '.join(['SELECT u_component_of_wind.val, u_component_of_wind.lat, u_component_of_wind.lon, u_component_of_wind.time, v_component_of_wind.val FROM u_component_of_wind '\
            'INNER JOIN v_component_of_wind '\
            'ON u_component_of_wind.lat == v_component_of_wind.lat',
            'u_component_of_wind.lon == v_component_of_wind.lon',
            'u_component_of_wind.time == v_component_of_wind.time '\
            'WHERE u_component_of_wind.lat >= ?',
            'u_component_of_wind.lat <= ?',
            'u_component_of_wind.lon >= ?',
            'u_component_of_wind.lon <= ?',
            'u_component_of_wind.time >= ?',
            'u_component_of_wind.time <= ?']) + ' ORDER BY u_component_of_wind.time, u_component_of_wind.lat, u_component_of_wind.lon ASC'

        db.execute(sql, tuple(map(str, [
                kwargs['south'],                kwargs['north'], 
                kwargs['west'],                 kwargs['east'], 
                dt_2_epoch(kwargs['start']), dt_2_epoch(kwargs['end'])
            ])))

        wind_u, lat, lon, epoch, wind_v = np.array(db.fetchall()).T
        val = np.sqrt(np.square(wind_u), np.square(wind_v))
        return np.array((val, lat, lon, epoch)).astype(float)


    def __str__(self):
        info = '\n'.join([
                "Era5 Global Dataset from Copernicus Climate Datastore.",
                "Combines model data with observations from across",
                "the world into a globally complete and consistent dataset",
                "\thttps://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels"])
        args = "(south, north, west, east, datetime, end)"
        return str_def(self, info, args)

