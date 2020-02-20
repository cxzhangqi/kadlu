"""
    API for Era5 dataset from Copernicus Climate Datastore

    Metadata regarding the dataset can be found here:
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
"""

import cdsapi
import numpy as np
import pygrib
import os
from os.path import isfile, dirname
from datetime import datetime, timedelta
import warnings
from configparser import ConfigParser

from kadlu.geospatial.data_sources.data_util    import              \
        storage_cfg,                                                \
        database_cfg,                                               \
        str_def,                                                    \
        dt_2_epoch,                                                 \
        epoch_2_dt,                                                 \
        serialized,                                                 \
        insert_hash


conn, db = database_cfg()
cfg = ConfigParser() 
cfg.read(dirname(dirname(dirname(dirname(__file__))))+'/config.ini')
try:
    c = cdsapi.Client(url=cfg['cdsapi']['url'], key=cfg['cdsapi']['url'])
except KeyError:
    try:
        c = cdsapi.Client()
    except Exception:
        raise Exception('CDS API has not been configured. obtain an API key '
                    'from the following URL and add it to kadlu/config.ini. '
                    'https://cds.climate.copernicus.eu/api-how-to')
    raise KeyError('CDS API has not been configured. obtain an API key '
                'from the following URL and add it to kadlu/config.ini. '
                'https://cds.climate.copernicus.eu/api-how-to')



def fetch_era5(var, kwargs):
    """ fetch global era5 data for specified variable and time range

        args:
            var: string
                the variable short name of desired wave parameter according to ERA5 docs
                the complete list can be found here (table 7 for wave params)
                https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation#ERA5datadocumentation-Temporalfrequency
            kwargs: dict
                keyword arguments passed from the Era5() class as a dictionary

        return:
            nothing
    """
    assert 6 == sum(map(lambda kw: kw in kwargs.keys(), 
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'

    if kwargs['start'] == kwargs['end']: kwargs['end'] += timedelta(hours=3)
    if serialized(kwargs, f'fetch_era5_{var}'): return False

    t = datetime(
            kwargs['start'].year,   kwargs['start'].month, 
            kwargs['start'].day,    kwargs['start'].hour
        )
    filenames = []

    thismonth = t.month
    while t <= kwargs['end']:
        fname = f"{storage_cfg()}ERA5_reanalysis_{var}_{t.strftime('%Y-%m-%d')}.grb2"
        filenames.append(fname)

        if isfile(fname):
            while (t.month == thismonth): t += timedelta(days=1)
            continue
            
        print(f"downloading {fname} from Copernicus Climate Data Store...")
        eom = min(datetime(t.year, t.month+1, 1)-timedelta(days=1), kwargs['end'])
        ddelta = (eom - t).days
        days = ([datetime(t.year, t.month, day).strftime('%d') 
                for day in range(t.day, eom.day)] 
                if t.day != eom.day else t.strftime('%d'))
        if ddelta == 0:
            eod = min(datetime(t.year, t.month, t.day, 23), kwargs['end'])
            hours = [datetime(t.year, t.month, t.day, hour).strftime('%H:00') 
                     for hour in range(t.hour, eod.hour)]
        else: 
            hours = [datetime(t.year, t.month, t.day, hour).strftime('%H:00') 
                     for hour in range(24)]

        c.retrieve('reanalysis-era5-single-levels', {
                   'product_type'  : 'reanalysis',
                   'format'        : 'grib',
                   'variable'      : var,
                   'year'          : t.strftime("%Y"),
                   'month'         : t.strftime("%m"),
                   'day'           : days,#t.strftime("%d"),
                   'time'          : hours#t.strftime("%H:00")
                }, fname)


        if not os.path.isfile(fname):
            warnings.warn(f'error fetching era5 data for {t}')

        while (t.month == thismonth): t += timedelta(days=1)

    t = datetime(
            kwargs['start'].year,   kwargs['start'].month, 
            kwargs['start'].day,    kwargs['start'].hour
        )
    for fname in filenames:
        grib = pygrib.open(fname)
        table = var[4:] if var[0:4] == '10m_' else var
        n1 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
        msgnum = 1
        rows = 0
        for msg in grib:
            print(f'processing msg {msgnum}/{grib.messages} for {var} {msg.validDate}', end='\r')
            msgnum += 1
            if msg.validDate < kwargs['start'] or msg.validDate > kwargs['end']: continue
            z, y, x = msg.data()

            if np.ma.is_masked(z):
                z2 = z[~z.mask].data
                y2 = y[~z.mask]
                x2 = x[~z.mask]
            else:  # wind data has no mask
                z2 = z.reshape(-1)
                y2 = y.reshape(-1)
                x2 = x.reshape(-1)

            # index coordinates and select query range subset
            latix = np.logical_and(y2 >= kwargs['south'], y2 <= kwargs['north'])
            lonix = np.logical_and(x2 >= kwargs['west']+180,  x2 <= kwargs['east']+180)
            ix = np.logical_and(latix, lonix)

            # build coordinate grid and insert into db
            grid = np.empty((len(z2[ix]), 5), dtype=object)
            grid[:,0] = z2[ix]
            grid[:,1] = y2[ix]
            grid[:,2] = ((x2[ix] + 180) % 360) - 180
            grid[:,3] = dt_2_epoch([msg.validDate for item in z2[ix]])
            grid[:,4] = ['era5' for item in z2[ix]]
            db.executemany(f"INSERT OR IGNORE INTO {table} VALUES (?,?,?,CAST(? AS INT),?)", grid)
            rows += len(grid)

        n2 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
        db.execute("COMMIT")
        conn.commit()

        print(f"\nprocessed and inserted {n2-n1} rows for {msg.validDate}.\t"
              f"{rows - (n2-n1)} duplicate rows ignored")

        thismonth = t.month
        while (t.month == thismonth): t += timedelta(days=1)

    insert_hash(kwargs, f'fetch_era5_{var}')
    return True

"""
def msg_parser(var, grib, msg, kwargs)
    print(f'processing msg {msgnum}/{grib.messages} for {var} {msg.validDate}', end='\r')
    if msg.validDate < kwargs['start'] or msg.validDate > kwargs['end']: return
    z, y, x = msg.data()

    if np.ma.is_masked(z):
        z2 = z[~z.mask].data
        y2 = y[~z.mask]
        x2 = x[~z.mask]
    else:  # wind data has no mask
        z2 = z.reshape(-1)
        y2 = y.reshape(-1)
        x2 = x.reshape(-1)

    # build coordinate grid and insert into db
    grid = np.empty((len(z2), 5), dtype=object)
    grid[:,0] = z2
    grid[:,1] = y2
    grid[:,2] = ((x2 + 180) % 360) - 180
    grid[:,3] = dt_2_epoch([msg.validDate for item in z2])
    grid[:,4] = ['era5' for item in z2]
    db.executemany(f"INSERT OR IGNORE INTO {table} VALUES (?,?,?,CAST(? AS INT),?)", grid)
    conn.commit()

    ### end of testing ###
"""

def load_era5(var, kwargs):
    if 'time' in kwargs.keys():
        assert False, 'nearest time search not implemented'

    assert 6 == sum(map(lambda kw: kw in kwargs.keys(),
        ['south', 'north', 'west', 'east', 'start', 'end'])), 'malformed query'

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
            dt_2_epoch(kwargs['start'])[0], dt_2_epoch(kwargs['end'])[0]
        ])))
    rowdata = np.array(db.fetchall(), dtype=object).T
    #assert len(rowdata) > 0, "no data found for query"
    if len(rowdata) == 0:
        warnings.warn('no records found, returning empty arrays')
        return np.array([ [], [], [], [], [] ])

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
    def fetch_wind(self, **kwargs):
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
    def load_wind(self, **kwargs):
        #wind_u = load_era5('10m_u_component_of_wind', kwargs)
        #wind_v = load_era5('10m_v_component_of_wind', kwargs)
        #ix = np.unique(wind_u[1:] == wind_v[1:], axis=2)
        #wind_uv = wind_u.copy()
        #wind_uv[0] = tuple(zip(wind_u[0], wind_v[0]))
        #wind_uv[0] = np.sqrt(np.square(wind_u[0]), np.square(wind_v[0]))

        #table = var[4:] if var[0:4] == '10m_' else var  # table cant start with int
        sql = ' AND '.join(['SELECT * FROM u_component_of_wind '\
            'INNER JOIN v_component_of_wind '\
            'ON u_component_of_wind.lat == v_component_of_wind.lat',
            'u_component_of_wind.lon == v_component_of_wind.lon',
            'u_component_of_wind.time == v_component_of_wind.time WHERE u_component_of_wind.lat >= ?',
            'u_component_of_wind.lat <= ?',
            'u_component_of_wind.lon >= ?',
            'u_component_of_wind.lon <= ?',
            'u_component_of_wind.time >= ?',
            'u_component_of_wind.time <= ?']) + ' ORDER BY time, lat, lon ASC'
        db.execute(sql, tuple(map(str, [
                kwargs['south'],                kwargs['north'], 
                kwargs['west'],                 kwargs['east'], 
                dt_2_epoch(kwargs['start'])[0], dt_2_epoch(kwargs['end'])[0]
            ])))
        wind_u, lat, lon, epoch, _, wind_v, _, _, _, _ = np.array(db.fetchall()).T
        val = np.sqrt(np.square(wind_u.astype(float)), np.square(wind_v.astype(float)))
        return np.array((val, lat, lon, epoch)).astype(float)

    def __str__(self):
        info = '\n'.join([
                "Era5 Global Dataset from Copernicus Climate Datastore",
                "\thttps://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview"
            ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(), end=datetime())"
        return str_def(self, info, args)

