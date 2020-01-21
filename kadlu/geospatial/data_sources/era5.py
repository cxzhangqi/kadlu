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
import warnings

from kadlu.geospatial.data_sources.fetch_util import \
storage_cfg, database_cfg, str_def, plot_sample_grib, dt_2_epoch, epoch_2_dt


conn, db = database_cfg()


def fetchname(var, time):
    return f"ERA5_reanalysis_{var}_{time.strftime('%Y-%m-%d_%Hh')}.grb2"


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
    assert 2 == sum(map(lambda kw: kw in kwargs.keys(), 
        ['start', 'end'])), 'malformed query'

    t = datetime(
            kwargs['start'].year,   kwargs['start'].month, 
            kwargs['start'].day,    kwargs['start'].hour
        )

    while t <= kwargs['end']:
        fname = f"{storage_cfg()}{fetchname(var, t)}"
        print(f"downloading {fetchname(var, t)} from Copernicus Climate Data Store...")
        c = cdsapi.Client()
        try:
            c.retrieve('reanalysis-era5-single-levels', {
                'product_type'  : 'reanalysis',
                'format'        : 'grib',
                'variable'      : var,
                'year'          : t.strftime("%Y"),
                'month'         : t.strftime("%m"),
                'day'           : t.strftime("%d"),
                'time'          : t.strftime("%H:00")
            }, fname)
        except Exception:
            warnings.warn(f"No data for {t.ctime()}")
            t += timedelta(hours=1)
            continue

        grib = pygrib.open(fname)
        assert(grib.messages == 1)
        msg = grib[1]
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
        grid[:,3] = dt_2_epoch([t for item in z2])
        grid[:,4] = ['era5' for item in z2]

        table = var[4:] if var[0:4] == '10m_' else var
        n1 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
        db.executemany(f"INSERT OR IGNORE INTO {table} VALUES (?,?,?,CAST(? AS INT),?)", grid)
        n2 = db.execute(f"SELECT COUNT(*) FROM {table}").fetchall()[0][0]
        db.execute("COMMIT")
        conn.commit()

        print(f"processed and inserted {n2-n1} rows. "
              f"{len(grid) - (n2-n1)} duplicate rows ignored\n")

        t += timedelta(hours=1)

    return 


def load_era5(var, kwargs):
    if 'time' in kwargs.keys():
        assert False, 'nearest time search not implemented yet'

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

    return np.array((val, lat, lon, epoch_2_dt(epoch)), dtype=object)


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
        fetch_era5('significant_height_of_combined_wind_waves_and_swell', kwargs)
        return
    def fetch_wavedirection(self, **kwargs):
        fetch_era5('mean_wave_direction', kwargs)
        return
    def fetch_waveperiod(self, **kwargs):
        fetch_era5('mean_wave_period', kwargs)
        return
    def fetch_wind_u(self, **kwargs):
        fetch_era5('10m_u_component_of_wind', kwargs)
        return
    def fetch_wind_v(self, **kwargs):
        fetch_era5('10m_v_component_of_wind', kwargs)
        return
    def fetch_wind(self, **kwargs):
        fetch_era5('10m_u_component_of_wind', kwargs)
        fetch_era5('10m_v_component_of_wind', kwargs)
        return

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
        wind_u = load_era5('10m_u_component_of_wind', kwargs)
        wind_v = load_era5('10m_v_component_of_wind', kwargs)
        wind_uv = wind_u.copy()
        wind_uv[0] = tuple(zip(wind_u[0], wind_v[0]))
        return wind_uv

    def __str__(self):
        info = '\n'.join([
                "Era5 Global Dataset from Copernicus Climate Datastore",
                "\thttps://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview"
            ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(), end=datetime())"
        return str_def(self, info, args)

