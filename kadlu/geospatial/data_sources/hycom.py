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
from datetime import datetime
import warnings


def fetchname(fetchvar, slices):
    """ build the query to slice the data from the dataset """
    slicer = lambda tup, step=1 : f"[{tup[0]}:{step}:{tup[1]}]"
    return f"{fetchvar}{''.join(map(slicer, slices))}"


def ll_2_xy(val, arr):
    """ converts lat/lon values to grid index (next nearest corner southwest of point) """
    if val > arr[-1]: return len(arr) - 1
    return np.max([0, np.nonzero(arr >= val)[0][0] - 1])


def dt_2_tslice(start, end):
    """ converts datetime to hycom time slice """
    assert(start >= datetime(1994, 1, 1))
    assert(end < datetime(2016, 1, 1))
    assert(start.year == end.year)
    warnings.warn("Hycom time conversion has a 3-hour margin of error. In the future, this should be mapped")
    
    def dt_2_t(time):
        seconds_delta = (time - datetime(time.year, 1, 1)).total_seconds()
        seconds_per_slice = 31536000.0 / 2860  # seconds per year / num of slices
        return int(seconds_delta / seconds_per_slice)
        
    return (dt_2_t(start), dt_2_t(end))


def fetch_grid():
    """ download the lat/lon arrays for grid indexing """
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
    def _load_grid(storage_cfg): 
        return np.load(f"{storage_cfg}hycom_lats.npy"), np.load(f"{storage_cfg}hycom_lons.npy")
    try:
        return _load_grid(storage_cfg())
    except FileNotFoundError:
        fetch_grid()
        return _load_grid(storage_cfg())


def fetch_hycom(year, slices, fetchvar):
    """
    download data from the hycom server

    year = '2015'       # string value between 2011? and 2015
    slices = [
        (0, 2),         # time: start, end 
        (0, 3),         # depth: top?, bottom?
        (800, 840),     # x grid index: lon min, lon max
        (900, 1000)     # y grid index: lat min, lat max
    ]
    fetchvar = 'salinity'

    returns: list of fetched filenames
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


def load_hycom(year, slices, fetchvar, lat, lon):
    """ load local data into memory as np arrays """
    fname = f"{storage_cfg()}hycom_{year}_{fetchname(fetchvar, slices)}.npy"
    if not isfile(fname): fetch_hycom(year, slices, fetchvar)
    data = np.load(fname)
    data[data <= -29999] = None  # Convert missing values to NaN
    return data, lat, lon


class Hycom():
    def fetch_salinity(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return fetch_hycom(
                    year=str(start.year),                                       # each calendar year requires separate query
                    slices=[
                        dt_2_tslice(start, end),                                # time range (must be same year)
                        (0, 39),                                                # depth range: get entire water column
                        (ll_2_xy(west,  self.lon), ll_2_xy(east,  self.lon)),   # tuple: (xmin, xmax)
                        (ll_2_xy(south, self.lat), ll_2_xy(north, self.lat))    # tuple: (ymin, ymax)
                    ],
                    fetchvar='salinity'
                )
    def fetch_temp(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end),
                        (0, 39),
                        (ll_2_xy(west,  self.lon), ll_2_xy(east,  self.lon)),
                        (ll_2_xy(south, self.lat), ll_2_xy(north, self.lat))
                    ],
                    fetchvar='water_temp'
                )
    def fetch_water_u(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end),
                        (0, 39),
                        (ll_2_xy(west,  self.lon), ll_2_xy(east,  self.lon)),
                        (ll_2_xy(south, self.lat), ll_2_xy(north, self.lat))
                    ],
                    fetchvar='water_u'
                )
    def fetch_water_v(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return fetch_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end),
                        (0,39),
                        (ll_2_xy(west,  self.lon),  ll_2_xy(east,  self.lon)),
                        (ll_2_xy(south, self.lat),  ll_2_xy(north, self.lat))
                    ],
                    fetchvar='water_v'
                )

    def load_salinity(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return load_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end),
                        (0, 39),
                        (ll_2_xy(west,  self.lon),  ll_2_xy(east,  self.lon)),
                        (ll_2_xy(south, self.lat),  ll_2_xy(north, self.lat))
                    ],
                    fetchvar='salinity',
                    lat=self.lat,
                    lon=self.lon
                )
    def load_temp(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return load_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end),
                        (0, 39),
                        (ll_2_xy(west,  self.lon),  ll_2_xy(east,  self.lon)),
                        (ll_2_xy(south, self.lat),  ll_2_xy(north, self.lat))
                    ],
                    fetchvar='water_temp',
                    lat=self.lat,
                    lon=self.lon
                )
    def load_water_u(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return load_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end),
                        (0, 39),
                        (ll_2_xy(west,  self.lon),  ll_2_xy(east,  self.lon)),
                        (ll_2_xy(south, self.lat),  ll_2_xy(north, self.lat))
                    ],
                    fetchvar='water_u',
                    lat=self.lat,
                    lon=self.lon
                )
    def load_water_v(self, south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31)): 
        return load_hycom(
                    year=str(start.year),
                    slices=[
                        dt_2_tslice(start, end),
                        (0, 39),
                        (ll_2_xy(west,  self.lon),  ll_2_xy(east,  self.lon)),
                        (ll_2_xy(south, self.lat),  ll_2_xy(north, self.lat))
                    ],
                    fetchvar='water_v',
                    lat=self.lat,
                    lon=self.lon
                )

    def __str__(self):
        info = '\n'.join([
                "Native hycom .[ab] data converted to NetCDF at the Naval Research Laboratory,",
                "interpolated to a uniform 0.08° between 40°S-40°N (0.04° poleward of these latitudes),",
                "and interpolated to 40 standard z-levels.",
                "\thttps://www.hycom.org/data/glbv0pt08"
            ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(1994, 1, 1), end=datetime(2015, 12, 31))"
        return fetch_util.str_def(self, info, args)
    def __init__(self):
        self.lat, self.lon = load_grid()


"""
print(Hycom())

mahone bay test area:
south =  44.4
north =  44.7
west  = -64.4
east  = -63.8
time  = datetime(2015, 1, 1, 0, 0)

source = Hycom()
salinity, lat, lon  = Hycom().load_salinity(south=south, north=north, west=west, east=east, time=datetime) 
temp, lat, lon      = source.load_temp(south=south, north=north, west=west, east=east, time=datetime) 
water_u, lat, lon   = Hycom().load_water_u(south=south, north=north, west=west, east=east, time=datetime) 
water_v, lat, lon   = Hycom().load_water_v(south=south, north=north, west=west, east=east, time=datetime) 
"""

