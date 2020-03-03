"""
    API for Era5 dataset from Copernicus Climate Datastore
     
    Metadata regarding the dataset can be found here:
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
"""

import os
import warnings
from os.path import isfile, dirname
from configparser import ConfigParser
from datetime import datetime, timedelta

import cdsapi
import pygrib
import numpy as np

import kadlu.geospatial.data_sources.source_map
from kadlu.geospatial.data_sources.data_util    import              \
        database_cfg,                                               \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized,                                                 \
        dt_2_epoch,                                                 \
        epoch_2_dt,                                                 \
        dev_null,                                                   \
        str_def,                                                    \
        cfg


conn, db = database_cfg()
try: c = cdsapi.Client(url=cfg['cdsapi']['url'], key=cfg['cdsapi']['key'])
except KeyError:
    try: c = cdsapi.Client()
    except Exception:
        raise KeyError('CDS API has not been configured. obtain an API token '
                       'from the following URL and add it to kadlu/config.ini. '
                       'https://cds.climate.copernicus.eu/api-how-to')

era5_varmap = dict(zip(
        ('significant_height_of_combined_wind_waves_and_swell',
         'mean_wave_direction',
         'mean_wave_period',
         '10m_u_component_of_wind',
         '10m_v_component_of_wind'),
        ('waveheight', 'wavedir', 'waveperiod', 'windspeedU', 'windspeedV')))


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

    assert 6 == sum(map(lambda kw: kw in kwargs.keys(), 
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'
    t = datetime(kwargs['start'].year,   kwargs['start'].month,
                 kwargs['start'].day,    kwargs['start'].hour)
    assert kwargs['end'] - kwargs['start'] <= timedelta(days=1, hours=1), \
            'use fetch_handler for this instead'
        
    if serialized(kwargs, f'fetch_era5_{era5_varmap[var]}'): return False

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

    if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
    n1 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
    db.executemany(f"INSERT OR IGNORE INTO {table} "
                   f"VALUES (?,?,?,CAST(? AS INT),?)", agg.T)
    n2 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
    db.execute("COMMIT")
    conn.commit()
    insert_hash(kwargs, f'fetch_era5_{era5_varmap[var]}')
    if 'lock' in kwargs.keys(): kwargs['lock'].release()

    print(f"ERA5 {msg.validDate.date().isoformat()} {var}: "
          f"processed and inserted {n2-n1} rows. "
          f"{len(agg[0])- (n2-n1)} duplicates ignored")

    return True


def load_era5(var, kwargs):
    if 'time' in kwargs.keys():
        assert False, 'nearest time search not implemented'

    assert 6 == sum(map(lambda kw: kw in kwargs.keys(),
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'

    # check for missing data
    kadlu.geospatial.data_sources.source_map.fetch_handler(
            era5_varmap[var], 'era5', parallel=1, **kwargs)

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
    assert len(rowdata) > 0, "no data found for query"
    val, lat, lon, epoch, source = rowdata 

    return np.array((val, lat, lon, epoch), dtype=np.float)


class Era5():
    """ collection of module functions for fetching and loading. 
    
        fetch function kwargs:
            start:
                datetime for beginning of time range to be fetched
            end:
                datetime for end of time range to be fetched

        load function kwargs:
            start:
            end:
    """

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
        fetch_era5('10m_u_component_of_wind', kwargs)
        fetch_era5('10m_v_component_of_wind', kwargs)

        sql = ' AND '.join(['SELECT * FROM u_component_of_wind '\
            'INNER JOIN v_component_of_wind '\
            'ON u_component_of_wind.lat == v_component_of_wind.lat',
            'u_component_of_wind.lon == v_component_of_wind.lon',
            'u_component_of_wind.time == v_component_of_wind.time '\
            'WHERE u_component_of_wind.lat >= ?',
            'u_component_of_wind.lat <= ?',
            'u_component_of_wind.lon >= ?',
            'u_component_of_wind.lon <= ?',
            'u_component_of_wind.time >= ?',
            'u_component_of_wind.time <= ?']) + ' ORDER BY time, lat, lon ASC'
        db.execute(sql, tuple(map(str, [
                kwargs['south'],                kwargs['north'], 
                kwargs['west'],                 kwargs['east'], 
                dt_2_epoch(kwargs['start']), dt_2_epoch(kwargs['end'])
            ])))
        wind_u, lat, lon, epoch, _, wind_v, _, _, _, _ = np.array(db.fetchall()).T
        val = np.sqrt(np.square(wind_u.astype(float)), np.square(wind_v.astype(float)))
        return np.array((val, lat, lon, epoch)).astype(float)

    def __str__(self):
        info = '\n'.join([
                "Era5 Global Dataset from Copernicus Climate Datastore.",
                "Combines model data with observations from across",
                "the world into a globally complete and consistent dataset",
                "\thttps://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels"])
        args = "(south, north, west, east, datetime, end)"
        return str_def(self, info, args)

