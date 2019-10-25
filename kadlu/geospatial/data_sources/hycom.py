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


def fetchname(fetchvar, slices, steps=(1, 1, 1, 1)):
    """ build the query to slice the data from the dataset """
    slicer = lambda tup, step : f"[{tup[0]}:{step}:{tup[1]}]"
    sliced = ''.join(map(slicer, slices, steps))
    return f"{fetchvar}{sliced}"


def index(val, sorted_arr):
    """ converts value in coordinate array to grid index """
    if val > sorted_arr[-1]: return len(sorted_arr) - 1
    return np.nonzero(sorted_arr >= val)[0][0]


def dt_2_tslice(start, end, epoch):
    """ converts datetime range to hycom time slice """
    assert(start >= datetime(1994, 1, 1))
    assert(end < datetime(2016, 1, 1))
    assert(start.year == end.year)
    return (index(dt_2_epoch(start), epoch[str(start.year)]), 
            index(dt_2_epoch(end),   epoch[str(end.year)]))


def fetch_grid():
    """ download lat/lon arrays for grid indexing """
    print("Fetching Hycom lat/lon grid arrays...")
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
    """ fetch timestamps from hycom """
    timestamps = {}

    for year in map(lambda y: f"{y}", range(1994, 2016)):
        url = f"{hycom_src}/{year}.ascii?time"
        time_netcdf = requests.get(url)
        assert(time_netcdf.status_code == 200)
        meta, data = time_netcdf.text.split\
        ("---------------------------------------------\n")
        time_csv = data.split("\n\n")[:-1][0]
        timestamps[year] = np.array(time_csv.split("\n")[1].split(', ')[1:], dtype=float)
        time.sleep(0.5)

    np.save(f"{storage_cfg()}hycom_times_dict.npy", timestamps)
    return


def load_times():
    """ put timestamps into memory """
    if not isfile(f"{storage_cfg()}hycom_times_dict.npy"): fetch_times()
    return np.load(f"{storage_cfg()}hycom_times_dict.npy", allow_pickle=True).item()


def load_depth():
    """ return depth values array for indexing """
    return np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0,
        30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0,
        150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0, 800.0,
        900.0, 1000.0, 1250.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0])


def fetch_hycom(year, slices, fetchvar, lat, lon, epoch, depth):
    """ download data from hycom , prepare it, and load into db

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
            fetchvar: string
                variable to be fetched. complete list of variables here
                https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html
            lat: array
                the first array returned by load_grid()
                used as a Hycom() class attribute for optimization
            lon: array
                the second array returned by load_grid()
                used as a Hycom() class attribute for optimization
            epoch: dictionary
                dictionary of timestamps. a year string passed as
                dictionary key returns a numpy array of datetimes.
                this dict is returned by load_times()
            depth: array
                array returned by load_depth()
                used as a Hycom() class attribute for optimization

        stores data in geospatial database and returns nothing.
        displays status message to standard output
    """

    t1 = datetime.now()

    # generate request
    n = reduce(np.multiply, map(lambda s : s[1] - s[0] +1, slices))
    assert n > 0, f"{n} records available within query boundaries {slices}"
    print(f"downloading {n} {fetchvar} values from hycom...")
    src = f"{hycom_src}/{year}.ascii?"
    payload_netcdf= requests.get(f"{src}{fetchname(fetchvar, slices)}")
    assert payload_netcdf.status_code == 200, "couldn't access hycom server"

    t2 = datetime.now()

    # parse response into numpy array
    meta, data = payload_netcdf.text.split\
    ("---------------------------------------------\n")
    arrs = data.split("\n\n")[:-1]
    shape_str, payload = arrs[0].split("\n", 1)
    assert(shape_str[0:len(fetchvar)] == fetchvar)
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
            for x in lon        [slices[2][0] : slices[2][1] +1]
            for y in lat        [slices[3][0] : slices[3][1] +1]])
    grid[:,0] = np.reshape(cube, flatten)
    grid = grid[grid[:,0] > -30000]

    # batch database insertion ignoring duplicates
    n1 = db.execute(f"SELECT COUNT(*) FROM {fetchvar}").fetchall()[0][0]
    db.executemany(f"INSERT OR IGNORE INTO {fetchvar} VALUES (?,?,?,?,?,?)", grid)
    n2 = db.execute(f"SELECT COUNT(*) FROM {fetchvar}").fetchall()[0][0]
    db.execute("COMMIT")
    conn.commit()

    t3 = datetime.now()

    print(f"downloaded in {(t2-t1).seconds}.{str((t2-t1).microseconds)[0:3]}s. "
          f"parsed and inserted {n2 - n1} rows in "
          f"{(t3-t2).seconds}.{str((t3-t2).microseconds)[0:3]}s\n"
          f"{n - len(grid)} null values removed, "
          f"{len(grid) - (n2 - n1)} duplicate rows ignored\n")

    return


def load_hycom(fetchvar, south, north, west, east, start, end, top, bottom):
    """ load hycom data from local database

        args:
            fetchvar:
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
    # validate and sanitize user input
    south, north, west, east = list(map(float, [south, north, west, east]))
    assert(south < north)
    assert west < east, "queries spanning antimeridian have not been implemented yet"
    assert(top < bottom)
    assert(start < end)

    # query db
    db.execute(' AND '.join([
            f"SELECT * FROM {fetchvar} WHERE lat <= {north}",
            f"lat >= {south}",
            f"lon >= {west}",
            f"lon <= {east}",
            f"time >= {dt_2_epoch(start)[0] }",
            f"time <= {dt_2_epoch(end)[0]   }",
            #f"depth >= {top}",
            #f"depth <= {bottom}",
            f"source == 'hycom'"]))

    # transpose grid and convert epochs to datetime
    data = np.array(db.fetchall(), dtype=object).T
    assert len(data) > 0, "no records found"
    data[3] = epoch_2_dt(data[3])

    return data[0], data[1], data[2], data[3], data[4]


class Hycom():
    """ collection of module functions for fetching and loading. 
    abstracted to include a fetch and load function for every variable
    fetched from the source

        attributes:
            lat, lon: arrays
                spatial grid arrays. used to convert grid index to 
                coordinate value 
            times_dict: dictionary
                dictionary containing temporal grid arrays
                used to convert time index to datetime
                a year string key between 1994 and 2015 holds a numpy array
                of datetimes
            depth: array
                array of depths. used to convert depth index to value
    """

    def __init__(self):
        self.lat, self.lon = load_grid()
        self.times_dict = load_times()
        self.depth = load_depth()

    def fetch_salinity(self, south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=0): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.times_dict),
                        (index(top,   self.depth), index(bottom, self.depth)),
                        (index(west,  self.lon),   index(east,   self.lon)),
                        (index(south, self.lat),   index(north,  self.lat))
                    ],
                    fetchvar='salinity',
                    lat=self.lat,
                    lon=self.lon,
                    epoch=self.times_dict,
                    depth=self.depth
                )

    def fetch_temp(self, south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=0): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.times_dict),
                        (index(top,   self.depth), index(bottom, self.depth)),
                        (index(west,  self.lon),   index(east,   self.lon)),
                        (index(south, self.lat),   index(north,  self.lat))
                    ],
                    fetchvar='water_temp',
                    lat=self.lat,
                    lon=self.lon,
                    epoch=self.times_dict,
                    depth=self.depth)

    def fetch_water_u(self, south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=0): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.times_dict),
                        (index(top,   self.depth), index(bottom, self.depth)),
                        (index(west,  self.lon),   index(east,   self.lon)),
                        (index(south, self.lat),   index(north,  self.lat))
                    ],
                    fetchvar='water_u',
                    lat=self.lat,
                    lon=self.lon,
                    epoch=self.times_dict,
                    depth=self.depth)

    def fetch_water_v(self, south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=0): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.times_dict),
                        (index(top,   self.depth), index(bottom, self.depth)),
                        (index(west,  self.lon),   index(east,  self.lon)),
                        (index(south, self.lat),   index(north, self.lat))
                    ],
                    fetchvar='water_v',
                    lat=self.lat,
                    lon=self.lon,
                    epoch=self.times_dict,
                    depth=self.depth)

    def load_salinity(self, south=-90, north=90, west=-180, east=180, 
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=5000): 
        return load_hycom(fetchvar='salinity', 
                south=south, north=north, west=west, east=east, 
                start=start, end=end, top=top, bottom=bottom)

    def load_temp(self, south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=5000): 
        return load_hycom(fetchvar='water_temp',
                south=south, north=north, west=west, east=east, 
                start=start, end=end, top=top, bottom=bottom)

    def load_water_u(self, south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=5000): 
        return load_hycom(fetchvar='water_u',
                south=south, north=north, west=west, east=east, 
                start=start, end=end, top=top, bottom=bottom)

    def load_water_v(self, south=-90, north=90, west=-180, east=180,
            start=datetime(2000, 1, 1), end=datetime(2000, 1, 2),
            top=0, bottom=5000): 
        return load_hycom(fetchvar='water_v',
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

