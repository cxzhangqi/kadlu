"""
    Kadlu API for HYCOM data source

    Data source:
        https://www.hycom.org/data/glbv0pt08
    Web interface for hycom data retrieval:
        https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html
    Example of GET query for salinity:
        https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.ascii?salinity[0:1:2][0:1:3][800:1:830][900:1:940]

    Oliver Kirsebom
    Casey Hilliard
    Matthew Smith 
"""

import numpy as np
import requests
from os.path import isfile
from kadlu.geospatial.data_sources import fetch_util
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
from datetime import datetime, timedelta
import time
import warnings


def fetchname(fetchvar, slices):
    """ build the query to slice the data from the dataset """
    slicer = lambda tup, step=1 : f"[{tup[0]}:{step}:{tup[1]}]"
    return f"{fetchvar}{''.join(map(slicer, slices))}"


def index(val, sorted_arr):
    """ converts lat/lon/time values to grid index (next nearest index >= value) """
    if val > sorted_arr[-1]: return len(sorted_arr) - 1
    return np.nonzero(sorted_arr >= val)[0][0]


def dt_2_tslice(start, end, times_dict):
    """ converts datetime range to hycom time slice """
    assert(start >= datetime(1994, 1, 1))
    assert(end < datetime(2016, 1, 1))
    assert(start.year == end.year)
    return (index(start, times_dict[str(start.year)]), index(end, times_dict[str(end.year)]))


def fetch_grid():
    """ download the lat/lon arrays for grid indexing """
    print("Fetching Hycom lat/lon grid arrays...")
    url = "https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.ascii?lat%5B0:1:3250%5D,lon%5B0:1:4499%5D"
    grid_ascii = requests.get(url)
    assert(grid_ascii.status_code == 200)

    meta, data = grid_ascii.text.split("---------------------------------------------\n")
    lat_csv, lon_csv = data.split("\n\n")[:-1]
    lat = np.array(lat_csv.split("\n")[1].split(", "), dtype=np.float)
    lon = np.array(lon_csv.split("\n")[1].split(", "), dtype=np.float)

    np.save(f"{storage_cfg()}hycom_lats.npy", lat, allow_pickle=False)
    np.save(f"{storage_cfg()}hycom_lons.npy", lon, allow_pickle=False)
    return


def load_grid():
    """ put downloaded grid into memory """
    if not isfile(f"{storage_cfg()}hycom_lats.npy"): fetch_grid()
    return np.load(f"{storage_cfg()}hycom_lats.npy"), np.load(f"{storage_cfg()}hycom_lons.npy")


def fetch_times():
    """ fetches timestamps from server """
    print("Fetching Hycom timestamp arrays...")
    timestamps = {}
    for year in map(lambda y: f"{y}", range(1994, 2016)):
        url = f"https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/{year}.ascii?time"
        time_ascii = requests.get(url)
        assert(time_ascii.status_code == 200)
        meta, data = time_ascii.text.split("---------------------------------------------\n")
        time_csv = data.split("\n\n")[:-1][0]
        basetime = datetime(2000, 1, 1)
        tdelta = np.array(time_csv.split("\n")[1].split(", "), dtype=float)
        timestamps[year] = np.array([basetime + timedelta(hours=x) for x in tdelta])
        time.sleep(0.1)

    np.save(f"{storage_cfg()}hycom_timestamps_dict.npy", timestamps)
    return


def load_times():
    """ put timestamps into memory """
    if not isfile(f"{storage_cfg()}hycom_timestamps_dict.npy"): fetch_times()
    return np.load(f"{storage_cfg()}hycom_timestamps_dict.npy", allow_pickle=True).item()


def load_depth():
    """ return depth values for indexing """
    return np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0])


def fetch_hycom(year, slices, fetchvar):
    """ download data from the hycom server

    args:
        year: string
            string value between 2011? and 2015
        slices: list of tuples
            correct ordering for tuples is [time, depth, lon, lat]
            each tuple contains the start and end index of the variable to be sliced
            optionally, step size can be included as a third value of the tuple
            An example of what this list may look like:
            [
                (0, 2),         # time: start, end 
                (0, 3, 2),      # depth: top, bottom, step (this is the index of array returned by load_depth() )
                (800, 840),     # x grid index: lon min, lon max
                (900, 1000)     # y grid index: lat min, lat max
            ]
        fetchvar: string
            the variable to be fetched. a complete list of variables is found here
            https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html

    return: 
        filenames: list
            list of strings describing complete path of fetched data
    """
    # generate request
    source = f"https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/{year}.ascii?"
    payload_ascii = requests.get(f"{source}{fetchname(fetchvar, slices)}")
    assert(payload_ascii.status_code == 200)
    fname = f"hycom_{year}_{fetchname(fetchvar, slices)}.npy"
    print(f"Downloading {fname} from Hycom")

    # parse response into numpy array
    meta, data = payload_ascii.text.split("---------------------------------------------\n")
    arrs = data.split("\n\n")[:-1]
    shape_str, payload = arrs[0].split("\n", 1)
    assert(shape_str[0:len(fetchvar)] == fetchvar)
    shape = tuple([int(x) for x in shape_str.split("[", 1)[1][:-1].split("][")])  # black magic list coercion
    output = np.ndarray(shape, dtype=np.float)
    for arr in payload.split("\n"):
        ix_str, row_csv = arr.split(", ", 1)
        a, b, c = [int(x) for x in ix_str[1:-1].split("][")]
        output[a][b][c] = np.array(row_csv.split(", "), dtype=np.int)

    np.save(f"{storage_cfg()}{fname}", output, allow_pickle=False)
    return [f"{storage_cfg()}{fname}"]


def load_hycom(year, slices, fetchvar, lat, lon, time, depth):
    """ load local data into memory as np arrays 

    args:
        year: string
            string value between 2011? and 2015
        slices: list of tuples
            correct ordering for tuples is [time, depth, lon, lat]
            each tuple contains the start and end index of the variable to be sliced
            optionally, step size can be included as a third value of the tuple
            An example of what this list may look like:
            [
                (0, 2),         # time: start, end 
                (0, 3, 2),      # depth: top, bottom, step (this is the index of array returned by load_depth() )
                (800, 840),     # x grid index: lon min, lon max
                (900, 1000)     # y grid index: lat min, lat max
            ]
        fetchvar: string
            the variable to be fetched. a complete list of variables is found here
            https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.html
        lat: np.array
            the first array returned by load_grid()
            this is used as a Hycom() class attribute for optimization
        lon: np.array
            the second array returned by load_grid()
            this is used as a Hycom() class attribute for optimization
        time: dictionary
            the dictionary of time arrays returned by load_times()
            this is used as a Hycom() class attribute for optimization
        depth: np.array
            the array returned by load_depth()
            this is used as a Hycom() class attribute for optimization

    return: 
       tuple 
            the tuple contains:
            - a 4D array containing the requested variable data slice
            - latitudes describing the first dimension of the 4D array
            - longitudes describing the second dimension
            - datetimes describing the third dimension
            - depths describing the fourth dimension
    """
    fname = f"{storage_cfg()}hycom_{year}_{fetchname(fetchvar, slices)}.npy"
    if not isfile(fname): fetch_hycom(year, slices, fetchvar)
    data = np.load(fname)
    data[data <= -29999] = None  # Convert missing values to NaN

    return (data,       # 4 dimensional array 
            lat         [slices[3][0] : slices[3][1] +1], 
            lon         [slices[2][0] : slices[2][1] +1], 
            time[year]  [slices[0][0] : slices[0][1] +1],
            depth       [slices[1][0] : slices[1][1] +1])


class Hycom():
    """ collection of module functions for fetching and loading. abstracted to include a seperate function for each variable 

    attributes:
        lat, lon: np.array
            latitude and longitude grid arrays, used to convert grid index to lat/lon
        time: dictionary
            dictionary containing time grid arrays for each year available, used to convert time index to datetime
        depth: np.array
            array of depths, used to convert depth index to depth value
    """

    def __init__(self):
        self.lat, self.lon = load_grid()
        self.time = load_times()
        self.depth = load_depth()

    def fetch_salinity(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return fetch_hycom(
                    year=str(start.year),                                   # each calendar year requires separate query
                    slices=[
                        dt_2_tslice(start, end, self.time),                 # time range (must be same year)
                        (0, 39),                                            # depth range: get entire water column
                        (index(west,  self.lon), index(east,  self.lon)),   # tuple: (xmin, xmax)
                        (index(south, self.lat), index(north, self.lat))    # tuple: (ymin, ymax)
                    ],
                    fetchvar='salinity'
                )
    def fetch_temp(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.time),
                        (0, 39),
                        (index(west,  self.lon), index(east,  self.lon)),
                        (index(south, self.lat), index(north, self.lat))
                    ],
                    fetchvar='water_temp'
                )
    def fetch_water_u(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.time),
                        (0, 39),
                        (index(west,  self.lon), index(east,  self.lon)),
                        (index(south, self.lat), index(north, self.lat))
                    ],
                    fetchvar='water_u'
                )
    def fetch_water_v(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.time),
                        (0,39),
                        (index(west,  self.lon),  index(east,  self.lon)),
                        (index(south, self.lat),  index(north, self.lat))
                    ],
                    fetchvar='water_v'
                )

    def load_salinity(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return load_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.time),
                        (0, 39),
                        (index(west,  self.lon),  index(east,  self.lon)),
                        (index(south, self.lat),  index(north, self.lat))
                    ],
                    fetchvar='salinity',
                    lat=self.lat,
                    lon=self.lon,
                    time=self.time,
                    depth=self.depth
                )
    def load_temp(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return load_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.time),
                        (0, 39),
                        (index(west,  self.lon),  index(east,  self.lon)),
                        (index(south, self.lat),  index(north, self.lat))
                    ],
                    fetchvar='water_temp',
                    lat=self.lat,
                    lon=self.lon,
                    time=self.time,
                    depth=self.depth
                )
    def load_water_u(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return load_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.time),
                        (0, 39),
                        (index(west,  self.lon),  index(east,  self.lon)),
                        (index(south, self.lat),  index(north, self.lat))
                    ],
                    fetchvar='water_u',
                    lat=self.lat,
                    lon=self.lon,
                    time=self.time,
                    depth=self.depth
                )
    def load_water_v(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return load_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end, self.time),
                        (0, 39),
                        (index(west,  self.lon),  index(east,  self.lon)),
                        (index(south, self.lat),  index(north, self.lat))
                    ],
                    fetchvar='water_v',
                    lat=self.lat,
                    lon=self.lon,
                    time=self.time,
                    depth=self.depth
                )

    def __str__(self):
        info = '\n'.join([
                "Native hycom .[ab] data converted to NetCDF at the Naval Research Laboratory,",
                "interpolated to a uniform 0.08° between 40°S-40°N (0.04° poleward of these latitudes),",
                "and interpolated to 40 standard z-levels.",
                "Historical data available from 1994 to 2015 (inclusive).",
                "\thttps://www.hycom.org/data/glbv0pt08"
            ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31))"
        return fetch_util.str_def(self, info, args)

