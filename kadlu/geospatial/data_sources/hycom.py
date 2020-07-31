"""
    Kadlu API for HYCOM data source

    data source:
        https://www.hycom.org/data/glbv1pt08

    web interface for manual hycom data retrieval:
        https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html
"""

import time
import logging
import requests
import warnings
from functools import reduce
from datetime import datetime, timedelta
from os.path import isfile

import numpy as np

import kadlu.geospatial.data_sources.fetch_handler
from kadlu.geospatial.data_sources.data_util        import          \
        database_cfg,                                               \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized,                                                 \
        dt_2_epoch,                                                 \
        epoch_2_dt,                                                 \
        fmt_coords,                                                 \
        str_def,                                                    \
        index


hycom_src = "https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data"


hycom_varmap = dict(zip(
        ('salinity', 'water_temp', 'water_u', 'water_v'),
        ('salinity',       'temp', 'water_u', 'water_v')))


# database config
conn, db = database_cfg()
hycom_tables = ['hycom_salinity', 'hycom_water_temp', 'hycom_water_u', 'hycom_water_v']
for var in hycom_tables:
    db.execute(f'CREATE TABLE IF NOT EXISTS {var}'
                '( val     REAL NOT NULL,' 
                '  lat     REAL NOT NULL,' 
                '  lon     REAL NOT NULL,' 
                '  time    INT  NOT NULL,' 
                '  depth   INT  NOT NULL,' 
                '  source  TEXT NOT NULL )')

    db.execute(f'CREATE UNIQUE INDEX IF NOT EXISTS '
               f'idx_{var} on {var}(time, lon, lat, depth, val, source)')


def slices_str(var, slices, steps=(1, 1, 1, 1)):
    """ build the query to slice the data from the dataset """
    slicer = lambda tup, step : f"[{tup[0]}:{step}:{tup[1]}]"
    sliced = ''.join(map(slicer, slices, steps))
    return f"{var}{sliced}"


def fetch_grid():
    """ download lat/lon/time arrays for grid indexing """

    logging.info("fetching hycom lat/lon grid arrays...")
    url = f"{hycom_src}/2015.ascii?lat%5B0:1:3250%5D,lon%5B0:1:4499%5D"
    grid_netcdf = requests.get(url)
    assert(grid_netcdf.status_code == 200)

    meta, data = grid_netcdf.text.split\
    ("---------------------------------------------\n")
    lat_csv, lon_csv = data.split("\n\n")[:-1]
    lat = np.array(lat_csv.split("\n")[1].split(", "), dtype=np.float)
    lon = np.array(lon_csv.split("\n")[1].split(", "), dtype=np.float)

    np.save(f"{storage_cfg()}hycom_lats.npy", lat, allow_pickle=False)
    np.save(f"{storage_cfg()}hycom_lons.npy", lon, allow_pickle=False)

    ### END XY GRIDS ###
    """ fetch timestamps from hycom (epoch hours since 2000-01-01 00:00) """

    epoch = {}

    for year in map(str, range(1994, 2016)):
        url = f"{hycom_src}/{year}.ascii?time"
        time_netcdf = requests.get(url)
        assert(time_netcdf.status_code == 200)
        meta, data = time_netcdf.text.split\
        ("---------------------------------------------\n")
        csv = data.split("\n\n")[:-1][0]
        epoch[year] = np.array(csv.split("\n")[1].split(', ')[1:], dtype=float)
        time.sleep(0.5)

    np.save(f"{storage_cfg()}hycom_epoch.npy", epoch)

    ### END TIME GRID ### 

    return


def load_grid():
    """ put spatial grid into memory """
    if not isfile(f"{storage_cfg()}hycom_lats.npy"): fetch_grid()
    return (np.load(f"{storage_cfg()}hycom_lats.npy"),
            np.load(f"{storage_cfg()}hycom_lons.npy"))


#def fetch_times():
#    return


def load_times():
    """ put timestamps into memory """
    if not isfile(f"{storage_cfg()}hycom_epoch.npy"): fetch_grid()
    return np.load(f"{storage_cfg()}hycom_epoch.npy", allow_pickle=True).item()


def load_depth():
    """ return depth values array for indexing """
    return np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0,
        30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0,
        150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0, 800.0,
        900.0, 1000.0, 1250.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0])


def fetch_hycom(self, var, year, slices, kwargs):
    """ download data from hycom, prepare it, and load into database

        args:
            year: string
                string value between 1994 and 2016
            slices: list of tuples
                correct ordering for tuples is [epoch, depth, lon, lat]
                each tuple contains the start and end grid index of the
                variable to be sliced. an example of the slices list:
                slices = [
                    (0, 2),         # time: start, end 
                    (0, 3),         # depth: top, bottom
                    (800, 840),     # x grid index: xmin, xmax (lon)
                    (900, 1000)     # y grid index: ymin, ymax (lat)
                ]
            var: string
                variable to be fetched. complete list of variables here
                https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html
            lat: array
                the first array returned by load_grid()
            lon: array
                the second array returned by load_grid()
            epoch: dictionary
                dictionary containing temporal grid arrays
                used to convert epoch index to datetime
                a year string key between 1994 and 2015 holds a numpy array
                of datetimes
            depth: array
                array returned by load_depth()

        return: nothing
    """

    # generate request
    t1 = datetime.now()
    url = f"{hycom_src}/{year}.ascii?{slices_str(var, slices)}"
    with requests.get(url, stream=True) as payload_netcdf:
        assert payload_netcdf.status_code == 200, "couldn't access hycom server"
        meta, data = payload_netcdf.text.split\
        ("---------------------------------------------\n")

    t2 = datetime.now()

    # parse response into numpy array
    arrs = data.split("\n\n")[:-1]
    shape_str, payload = arrs[0].split("\n", 1)
    shape = tuple([int(x) for x in shape_str.split("[", 1)[1][:-1].split("][")])
    cube = np.ndarray(shape, dtype=np.float)

    for arr in payload.split("\n"):
        ix_str, row_csv = arr.split(", ", 1)
        a, b, c = [int(x) for x in ix_str[1:-1].split("][")]
        cube[a][b][c] = np.array(row_csv.split(", "), dtype=np.int)

    # build coordinate grid, populate with values, adjust scaling, remove nulls
    flatten = reduce(np.multiply, map(lambda s : s[1] - s[0] +1, slices))
    add_offset =  20 if 'salinity' in var or 'water_temp' in var else 0
    null_value = -10 if 'salinity' in var or 'water_temp' in var else -30
    grid = np.array([(None, y, x, t, d, 'hycom') 
            for t in self.epoch[year][slices[0][0] : slices[0][1] +1]
            for d in self.depth      [slices[1][0] : slices[1][1] +1]
            for y in self.ygrid      [slices[2][0] : slices[2][1] +1]
            for x in self.xgrid      [slices[3][0] : slices[3][1] +1]])
    grid[:,0] = np.reshape(cube, flatten) * 0.001 + add_offset
    grid = grid[grid[:,0] != null_value]

    # batch database insertion ignoring duplicates
    if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
    n1 = db.execute(f"SELECT COUNT(*) FROM hycom_{var}").fetchall()[0][0]
    db.executemany(f"INSERT OR IGNORE INTO hycom_{var} VALUES "
                    "(?, ?, ?, CAST(? AS INT), CAST(? AS INT), ?)", grid)
    n2 = db.execute(f"SELECT COUNT(*) FROM hycom_{var}").fetchall()[0][0]
    db.execute("COMMIT")
    conn.commit()
    insert_hash(kwargs, f'fetch_hycom_{hycom_varmap[var]}')
    if 'lock' in kwargs.keys(): kwargs['lock'].release()

    t3 = datetime.now()

    logging.info(f"HYCOM {epoch_2_dt([self.epoch[year][slices[0][0]]])[0].date().isoformat()} "
          f"{var}: downloaded {int(len(payload_netcdf.content)/8/1000)} Kb "
          f"in {(t2-t1).seconds}.{str((t2-t1).microseconds)[0:3]}s. "
          f"parsed and inserted {n2 - n1} rows in "
          f"{(t3-t2).seconds}.{str((t3-t2).microseconds)[0:3]}s. "
          f"{flatten - len(grid)} null values removed, "
          f"{len(grid) - (n2 - n1)} duplicates ignored")

    return


def load_hycom(self, var, kwargs):
    """ load hycom data from local database

        args:
            var:
                variable to be fetched. complete list of variables here
                https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html
            south, north: float
                ymin, ymax coordinate values. range: -90, 90
            west, east: float
                xmin, xmax coordinate values. range: -180, 180
            start, end: datetime
                temporal boundaries in datetime format

        return: 
            values: array
                values of the fetched variable
            lat: array
                y grid coordinates
            lon: array
                x grid coordinates
            epoch: array
                timestamps in epoch hours since jan 1 2000 
            depth: array
                measured in meters
    """
    # check if grids are initialized
    if not self.grids:
        self.ygrid, self.xgrid = load_grid()
        self.epoch = load_times()
        self.depth = load_depth()
        self.grids = [self.ygrid, self.xgrid, self.epoch, self.depth]

    # recursive function call for queries spanning antimeridian
    if (kwargs['west'] > kwargs['east']): 
        kwargs1 = kwargs.copy()
        kwargs2 = kwargs.copy()
        kwargs1['west'] = self.xgrid[0]
        kwargs2['east'] = self.xgrid[-1]
        return np.hstack((load_hycom(self, var, kwargs1), 
                          load_hycom(self, var, kwargs2)))

    # check for missing data
    kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
            hycom_varmap[var], 'hycom', **kwargs)

    # validate and execute query
    assert 8 == sum(map(lambda kw: kw in kwargs.keys(),
            ['south', 'north', 'west', 'east',
             'start', 'end', 'top', 'bottom'])), 'malformed query'

    assert kwargs['start'] <= kwargs['end']
    sql = ' AND '.join([f"SELECT * FROM hycom_{var} WHERE lat >= ?",
            'lat <= ?',
            'lon >= ?',
            'lon <= ?',
            'time >= ?',
            'time <= ?',
            'depth >= ?',
            'depth <= ?',
            "source == 'hycom' "]

        ) + 'ORDER BY time, depth, lat, lon ASC'
    db.execute(sql, tuple(map(str, [
            kwargs['south'],                kwargs['north'], 
            kwargs['west'],                 kwargs['east'],
            dt_2_epoch(kwargs['start']), dt_2_epoch(kwargs['end']),
            kwargs['top'],                  kwargs['bottom']
        ])))
    rowdata = np.array(db.fetchall(), dtype=object).T

    #assert len(rowdata) > 0, f'no data for query: {kwargs}'
    if len(rowdata) == 0:
        logging.warning(f'HYCOM {var}: no data found in region {fmt_coords(kwargs)}, returning empty arrays')
        return np.array([[],[],[],[],[]])

    return rowdata[0:5].astype(float)


def fetch_idx(self, var, kwargs): 
    """ convert user query to grid index slices, handle edge cases """

    def _idx(self, var, year, kwargs): 
        """ build indices for query and call fetch_hycom """
        haystack = np.array([self.epoch[year], self.depth, self.ygrid, self.xgrid])
        needles1 = np.array([
                dt_2_epoch(kwargs['start']),
                kwargs['top'], 
                kwargs['south'],
                kwargs['west']
            ])
        needles2 = np.array([
                dt_2_epoch(kwargs['end']),
                kwargs['bottom'],
                kwargs['north'],
                kwargs['east']
            ])
        slices = list(zip(
                map(index, needles1, haystack), 
                map(index, needles2, haystack)
            ))

        n = reduce(np.multiply, map(lambda s : s[1] - s[0] +1, slices))
        assert n > 0, f"{n} records available within query boundaries: {kwargs}"

        logging.info(f"HYCOM {kwargs['start'].date().isoformat()} "
              f"downloading {n} {var} values in region {fmt_coords(kwargs)}...")
        fetch_hycom(self=self, slices=slices, var=var, year=year, kwargs=kwargs)
        return

    assert kwargs['start'] <= kwargs['end']
    assert kwargs['start'] >  datetime(1994, 1, 1), 'data not available in this range'
    assert kwargs['end']   <  datetime(2016, 1, 1), 'data not available in this range'
    assert kwargs['south'] <= kwargs['north']
    assert kwargs['top']   <= kwargs['bottom']
    assert kwargs['start'] >= datetime(1994, 1, 1)
    assert kwargs['end']   <  datetime(2016, 1, 1)
    assert kwargs['end'] - kwargs['start'] <= timedelta(days=1), \
            "use fetch handler for this"

    # query local database for existing checksums
    if serialized(kwargs, f'fetch_hycom_{hycom_varmap[var]}'): return False
    if not serialized(seed='fetch_hycom_grid'):
        fetch_grid()
        insert_hash(seed='fetch_hycom_grid')

    if not self.grids:
        self.ygrid, self.xgrid = load_grid()
        self.epoch = load_times()
        self.depth = load_depth()
        self.grids = [self.ygrid, self.xgrid, self.epoch, self.depth]

    # if query spans antimeridian, make two seperate fetch requests
    year = str(kwargs['start'].year)
    if kwargs['west'] > kwargs['east']:
        logging.debug('splitting request')
        kwargs1, kwargs2 = kwargs.copy(), kwargs.copy()
        kwargs1['east'] = self.xgrid[-1]
        kwargs2['west'] = self.xgrid[0]
        if not serialized(kwargs1, f'fetch_hycom_{hycom_varmap[var]}'):
            _idx(self, var, year, kwargs1)
        if not serialized(kwargs2, f'fetch_hycom_{hycom_varmap[var]}'):
            _idx(self, var, year, kwargs2)
    else:
        _idx(self, var, year, kwargs)

    return True


class Hycom():
    """ collection of module functions for fetching and loading. 

        attributes:
            lat, lon: arrays
                spatial grid arrays. used to convert grid index to 
                coordinate value 
            epoch: dictionary
                dictionary containing temporal grid arrays
                used to convert epoch index to datetime
                a year string key between 1994 and 2015 holds a numpy array
                of datetimes
            depth: array
                array of depths. used to convert depth index to value
    """

    def __init__(self):
        #self.ygrid, self.xgrid = load_grid()
        #self.epoch = load_times()
        #self.depth = load_depth()
        self.grids = None
        pass

    def fetch_salinity(self, **kwargs): return fetch_idx(self,  'salinity',   kwargs)
    def fetch_temp    (self, **kwargs): return fetch_idx(self,  'water_temp', kwargs)
    def fetch_water_u (self, **kwargs): return fetch_idx(self,  'water_u',    kwargs)
    def fetch_water_v (self, **kwargs): return fetch_idx(self,  'water_v',    kwargs)
    def fetch_water_uv(self, **kwargs): return fetch_idx(self,  'water_u',    kwargs) and\
                                               fetch_idx(self,  'water_v',    kwargs)

    def load_salinity (self, **kwargs): return load_hycom(self, 'salinity',   kwargs)
    def load_temp     (self, **kwargs): return load_hycom(self, 'water_temp', kwargs)
    def load_water_u  (self, **kwargs): return load_hycom(self, 'water_u',    kwargs)
    def load_water_v  (self, **kwargs): return load_hycom(self, 'water_v',    kwargs)
    def load_water_uv (self, **kwargs):
        kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
                hycom_varmap['water_u'], 'hycom', parallel=1, **kwargs)
        kadlu.geospatial.data_sources.fetch_handler.fetch_handler(
                hycom_varmap['water_v'], 'hycom', parallel=1, **kwargs)

        sql = ' AND '.join(['SELECT hycom_water_u.val, hycom_water_u.lat, hycom_water_u.lon, hycom_water_u.time, hycom_water_u.depth, hycom_water_v.val FROM hycom_water_u '\
                'INNER JOIN hycom_water_v '\
                'ON hycom_water_u.lat == hycom_water_v.lat',
                'hycom_water_u.lon == hycom_water_v.lon',
                'hycom_water_u.time == hycom_water_v.time '\
                'WHERE hycom_water_u.lat >= ?',
                'hycom_water_u.lat <= ?',
                'hycom_water_u.lon >= ?',
                'hycom_water_u.lon <= ?',
                'hycom_water_u.time >= ?',
                'hycom_water_u.time <= ?',
                'hycom_water_u.depth >= ?',
                'hycom_water_u.depth <= ?',
            ]) + ' ORDER BY hycom_water_u.time, hycom_water_u.lat, hycom_water_u.lon ASC'

        db.execute(sql, tuple(map(str, [
                kwargs['south'],                kwargs['north'], 
                kwargs['west'],                 kwargs['east'], 
                dt_2_epoch(kwargs['start']),    dt_2_epoch(kwargs['end']),
                kwargs['top'],                  kwargs['bottom'],
            ])))
        qry = np.array(db.fetchall()).T

        if len(qry) == 0:
            logging.warning(f'HYCOM water_uv: no data found in region {fmt_coords(kwargs)}, returning empty arrays')
            return np.array([[],[],[],[],[]])

        logging.debug(f'{qry.shape}  {qry[:,0]}')
        water_u, lat, lon, epoch, depth, water_v = qry
        val = np.sqrt(np.square(water_u.astype(float)), np.square(water_v.astype(float)))
        return np.array((val, lat, lon, epoch, depth)).astype(float)


    def __str__(self):
        info = '\n'.join([
                "Native hycom .[ab] data converted to NetCDF at the Naval",
                "Research Laboratory, interpolated to 0.08° grid between",
                "40°S-40°N (0.04° poleward) containing 40 z-levels.",
                "Availability: 1994 to 2015",
                "\thttps://www.hycom.org/data/glbv0pt08" ])

        args = ("(south, north, west, east, "
                "start, end, top, bottom)")

        return str_def(self, info, args)

