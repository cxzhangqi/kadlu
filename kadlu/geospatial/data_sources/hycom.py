"""
    Kadlu API for HYCOM data source

    data source:
        https://www.hycom.org/data/glbv1pt08
    web interface for manual hycom data retrieval:
        https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html

    Oliver Kirsebom
    Casey Hilliard
    Matthew Smith 
"""

import numpy as np
from functools import reduce
import requests
import time
from datetime import datetime, timedelta
from os.path import isfile
import warnings

from kadlu.geospatial.data_sources.fetch_util import \
storage_cfg, database_cfg, dt_2_epoch, epoch_2_dt, str_def


hycom_src = "https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data"
conn, db = database_cfg()  # database connection and cursor objects 


def slices_str(var, slices, steps=(1, 1, 1, 1)):
    """ build the query to slice the data from the dataset """
    slicer = lambda tup, step : f"[{tup[0]}:{step}:{tup[1]}]"
    sliced = ''.join(map(slicer, slices, steps))
    return f"{var}{sliced}"


def index(val, sorted_arr):
    """ converts value in coordinate array to grid index """
    if val > sorted_arr[-1]: return len(sorted_arr) - 1
    return np.nonzero(sorted_arr >= val)[0][0]
  

def fetch_grid():
    """ download lat/lon arrays for grid indexing """
    print("fetching hycom lat/lon grid arrays...")
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
    return


def load_grid():
    """ put spatial grid into memory """
    if not isfile(f"{storage_cfg()}hycom_lats.npy"): fetch_grid()
    return (np.load(f"{storage_cfg()}hycom_lats.npy"),
            np.load(f"{storage_cfg()}hycom_lons.npy"))


def fetch_times():
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
    return


def load_times():
    """ put timestamps into memory """
    if not isfile(f"{storage_cfg()}hycom_epoch.npy"): fetch_times()
    return np.load(f"{storage_cfg()}hycom_epoch.npy", allow_pickle=True).item()


def load_depth():
    """ return depth values array for indexing """
    return np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0,
        30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0,
        150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0, 800.0,
        900.0, 1000.0, 1250.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0])


def fetch_hycom(*args, year, slices, var, lat, lon, epoch, depth, **kwargs):
    """ download data from hycom, prepare it, and load into db

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
                used as a Hycom() class attribute for optimization
            lon: array
                the second array returned by load_grid()
                used as a Hycom() class attribute for optimization
            epoch: dictionary
                dictionary containing temporal grid arrays
                used to convert epoch index to datetime
                a year string key between 1994 and 2015 holds a numpy array
                of datetimes
            depth: array
                array returned by load_depth()
                used as a Hycom() class attribute for optimization

        stores data in geospatial database and returns nothing.
        displays status message to standard output
    """
    # generate request
    n = reduce(np.multiply, map(lambda s : s[1] - s[0] +1, slices))
    assert n > 0, f"{n} records available within query boundaries {slices}"
    print(f"downloading {n} {var} values from hycom...")
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

    # build coordinate grid, populate with values, remove null entries
    flatten = reduce(np.multiply, map(lambda s : s[1] - s[0] +1, slices))
    grid = np.array([(None, y, x, t, d, 'hycom') 
            for t in epoch[year][slices[0][0] : slices[0][1] +1]
            for d in depth      [slices[1][0] : slices[1][1] +1]
            for y in lat        [slices[2][0] : slices[2][1] +1]
            for x in lon        [slices[3][0] : slices[3][1] +1]])
    grid[:,0] = np.reshape(cube, flatten)
    grid = grid[grid[:,0] > -30000]

    # batch database insertion ignoring duplicates
    n1 = db.execute(f"SELECT COUNT(*) FROM {var}").fetchall()[0][0]
    db.executemany(f"INSERT OR IGNORE INTO {var} VALUES (?,?,?,?,?,?)", grid)
    n2 = db.execute(f"SELECT COUNT(*) FROM {var}").fetchall()[0][0]
    db.execute("COMMIT")
    conn.commit()

    t3 = datetime.now()

    print(f"downloaded {len(payload_netcdf.content)/8/1000:.1f}Kb in "
          f"{(t2-t1).seconds}.{str((t2-t1).microseconds)[0:3]}s. "
          f"parsed and inserted {n2 - n1} rows in "
          f"{(t3-t2).seconds}.{str((t3-t2).microseconds)[0:3]}s\n"
          f"{n - len(grid)} null values removed, "
          f"{len(grid) - (n2 - n1)} duplicate rows ignored\n")

    return


def load_hycom(*args, var, south, north, west, east, start, end, top, bottom, **kwargs):
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
            time: array
                timestamps in datetime format
            depth: array
                measured in meters
    """
    south, north, west, east = map(float, [south, north, west, east])
    assert(south < north)
    assert(top <= bottom)
    assert(start < end)
    if 'limit' not in kwargs.keys() : kwargs['limit'] = '50000;--infinity'

    # recursive function call for queries spanning antimeridian
    if (west > east): return np.hstack(
            load_hycom(var, south, north, self.lon[0], east, 
                       start, end, top, bottom, limit=kwargs['limit']), 
            load_hycom(var, south, north, west, self.lon[-1], 
                       start, end, top, bottom, limit=kwargs['limit']))

    db.execute(' AND '.join([
            f"SELECT * FROM {var} WHERE lat >= ?",
                                       "lat <= ?",
                                       "lon >= ?",
                                       "lon <= ?",
                                       "time >= ?",
                                       "time <= ?",
                                       "depth >= ?",
                                       "depth <= ?",
                                      f"source == 'hycom' LIMIT {kwargs['limit']}"]),
            tuple(map(str, 
                [south, north, west, east, 
                dt_2_epoch(start)[0], dt_2_epoch(end)[0], 
                top, bottom]   )     )       )

    # transpose grid and convert epochs to datetime
    data = np.array(db.fetchall(), dtype=object).T
    assert len(data[0]) > 0, "no records found"

    return data[0:5]


def fetch_idx(self, var, qry): 
    """ convert user query to slices and handle edge cases """

    def _idx(self, var, year, qry): 
        """ build indices for query and call fetch_hycom """
        needles1 = np.array([dt_2_epoch(qry['start'])[0], qry['top'],
                             qry['south'], qry['west']])
        needles2 = np.array([dt_2_epoch(qry['end'])[0], qry['bottom'],
                             qry['north'], qry['east']])
        haystack = np.array([self.epoch[year], self.depth, self.lat, self.lon])
        slice1 = map(index, needles1, haystack)
        slice2 = map(index, needles2, haystack)
        slices = list(zip(slice1, slice2))
        
        return fetch_hycom(slices=slices, var=var, year=year, lat=self.lat,
                lon=self.lon, epoch=self.epoch, depth=self.depth)

    south, north, west, east = \
    map(float, [qry['south'], qry['north'], qry['west'], qry['east']])

    assert(south <= north)
    assert(qry['start'] >= datetime(1994, 1, 1))
    assert(qry['end']   <  datetime(2016, 1, 1))
    assert(qry['start'] <= qry['end'])
    assert(qry['top']   <= qry['bottom'])

    # TODO: 
    # if start.year != end.year:
    #     call _idx once per year
    assert qry['start'].year == qry['end'].year, \
            "hycom queries spanning multiple years are not supported yet"
    year = str(qry['start'].year)

    if west > east:
        qry1, qry2 = [qry.copy(), qry.copy()]
        qry1['east'] = self.lon[-1]
        qry2['west'] = self.lon[0]
        print('partitioning query boundaries at antimeridian')
        for qr in [qry1, qry2]: _idx(self, var, year, qr)
        return

    return _idx(self, var, year, qry)


class Hycom():
    """ collection of module functions for fetching and loading. 
    abstracted to include a fetch and load function for every variable
    fetched from the source

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
        self.lat, self.lon = load_grid()
        self.epoch = load_times()
        self.depth = load_depth()

    def fetch_salinity(self, **qry): return fetch_idx(self, 'salinity', qry)

    def fetch_temp    (self, **qry): return fetch_idx(self, 'water_temp', qry)

    def fetch_water_u (self, **qry): return fetch_idx(self, 'water_u', qry)

    def fetch_water_v (self, **qry): return fetch_idx(self, 'water_v', qry)

    def load_salinity(self,
            south=-90, north=90, west=-180, east=180, 
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=5000): 

        return load_hycom(var='salinity', 
                south=south, north=north, west=west, east=east, 
                start=start, end=end, top=top, bottom=bottom)

    def load_temp(self,
            south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=5000): 

        return load_hycom(var='water_temp',
                south=south, north=north, west=west, east=east, 
                start=start, end=end, top=top, bottom=bottom)

    def load_water_u(self,
            south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=5000): 

        return load_hycom(var='water_u',
                south=south, north=north, west=west, east=east, 
                start=start, end=end, top=top, bottom=bottom)

    def load_water_v(self,
            south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=5000): 

        return load_hycom(var='water_v',
                south=south, north=north, west=west, east=east, 
                start=start, end=end, top=top, bottom=bottom)

    def __str__(self):
        info = '\n'.join([
                "Native hycom .[ab] data converted to NetCDF at the Naval",
                "Research Laboratory, interpolated to a uniform 0.08° between",
                "40°S-40°N (0.04° poleward of these latitudes), and",
                "interpolated to 40 standard z-levels.",
                "Historical data available from 1994 to 2015 (inclusive).",
                "\thttps://www.hycom.org/data/glbv0pt08" ])

        args = ("(south, north, west, east, "
                "start, end, top, bottom)")

        return str_def(self, info, args)

