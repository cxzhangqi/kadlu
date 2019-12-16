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


def fetchname(wavevar, time):
    return f"ERA5_reanalysis_{wavevar}_{time.strftime('%Y-%m-%d_%Hh')}.grb2"


def fetch_era5(wavevar, start, end):
    """ return list of filenames containing era5 data for time range

    args:
        wavevar: string
            the variable short name of desired wave parameter according to ERA5 docs
            the complete list can be found here (table 7 for wave params)
            https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation#ERA5datadocumentation-Temporalfrequency
        start: datetime
            the beginning of the desired time range
        end: datetime
            the end of the desired time range

    return:
        fetchfiles: list
            a list of strings describing complete filepaths of downloaded data
    """
    time = datetime(start.year, start.month, start.day, start.hour)

    while time <= end:
        fname= f"{storage_cfg()}{fetchname(wavevar, time)}"
        print(f"downloading {fetchname(wavevar, time)} from Copernicus Climate Data Store...")
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

        grib = pygrib.open(fname)
        assert(grib.messages == 1)
        msg = grib[1]
        z, y, x = msg.data()
        #x -= 180

        # build coordinate grid and insert into db
        src = ['era5' for item in z[~z.mask]]
        t = dt_2_epoch([time for item in z[~z.mask]])
        grid = list(map(tuple, np.vstack((
                z[~z.mask].data, 
                y[~z.mask], 
                ((x[~z.mask] + 180) % 360 ) - 180, 
                t,
                src
            )).T
        ))
        n1 = db.execute(f"SELECT COUNT(*) FROM {wavevar}").fetchall()[0][0]
        db.executemany(f"INSERT OR IGNORE INTO {wavevar} VALUES (?,?,?,?,?)", grid)
        n2 = db.execute(f"SELECT COUNT(*) FROM {wavevar}").fetchall()[0][0]
        db.execute("COMMIT")
        conn.commit()

        print(f"{fname.split('/')[-1]} processed and inserted {n2-n1} rows. "
              f"{(z.shape[0]*z.shape[1])-len(z[~z.mask])} null values removed, "
              f"{len(grid) - (n2-n1)} duplicate rows ignored")

        time += timedelta(hours=1)

    return 

def load_era5(wavevar, start, end, south, north, west, east):
    qry = ' AND '.join([f"SELECT * FROM {wavevar} WHERE lat >= ?",
                                                       "lat <= ?",
                                                       "lon >= ?",
                                                       "lon <= ?",
                                                       "time >= ?",
                                                       "time <= ?"])
    db.execute(qry, tuple(map(str, 
            [south, north, west, east, dt_2_epoch(start)[0], dt_2_epoch(end)[0]]
        )))
    slices = np.array(db.fetchall(), dtype=object).T
    assert len(slices) > 0, "no data found for query"
    val, lat, lon, time, source = slices
    return val, lat, lon, epoch_2_dt(time)


class Era5():
    """ collection of module functions for fetching and loading. abstracted to include a seperate function for each variable """

    def fetch_windwaveswellheight(self, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_era5('significant_height_of_combined_wind_waves_and_swell', start, end)

    def fetch_wavedirection(self, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_era5('mean_wave_direction', start, end)

    def fetch_waveperiod(self, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_era5('mean_wave_period', start, end)

    def load_windwaveswellheight(self, south, north, west, east, 
            start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_era5('significant_height_of_combined_wind_waves_and_swell', start, end, south, north, west, east)
    
    def load_wavedirection(self, south, north, west, east,
            start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_era5('mean_wave_direction', start, end, south, north, west, east)
    
    def load_waveperiod(self, south, north, west, east,
            start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_era5('mean_wave_period', start, end, south, north, west, east)

    def __str__(self):
        info = '\n'.join([
                "Era5 Global Dataset from Copernicus Climate Datastore",
                "\thttps://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview"
            ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(), end=datetime())"
        return str_def(self, info, args)

"""

wavevar = 'significant_height_of_combined_wind_waves_and_swell'
start = datetime(2018, 1, 1, 0, 0, 0, 0)
end   = datetime(2018, 1, 1, 0, 0, 0, 0)

"""
