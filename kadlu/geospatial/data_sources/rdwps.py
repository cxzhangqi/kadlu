"""
    API for Regional Deterministic Wave Prediction System (RDWPS) data provided by Govt. of Canada

    Metadata regarding the dataset can be found here:
        https://weather.gc.ca/grib/grib2_RDWPS_e.html

    Oliver Kirsebom
    Casey Hilliard
    Matthew Smith
"""

import numpy as np
from datetime import datetime, timedelta
import os
import urllib.request
import pygrib
from kadlu.geospatial.data_sources import fetch_util 
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
import warnings


# dictionary to obtain reference level as per RDWPS nomenclature
# for more info see the following: https://weather.gc.ca/grib/grib2_RDWPS_e.html
from collections import defaultdict
region_strs = ['gulf-st-lawrence', 'superior', 'huron-michigan', 'erie', 'ontario']
default = defaultdict(lambda : 'SFC_0')
level_ref = {key : default for key in region_strs}
for reg in region_strs: 
    level_ref[reg]['UGRD'] = 'TGL_10'
    level_ref[reg]['VGRD'] = 'TGL_10'
level_ref['gulf-st-lawrence']['PKPER'] = 'TGL_0'
level_ref['gulf-st-lawrence']['PRMSL'] = 'MSL_0'


class Boundary():
    def __init__(self, south, north, west, east):
        self.south = south
        self.north = north
        self.west = west
        self.east = east
    def __eq__(self, other):
        """
        separating axis theorem to quickly check for intersecting regions
            https://en.wikipedia.org/wiki/Hyperplane_separation_theorem#Use_in_collision_detection
        """
        return not (self.east < other.west or 
                    self.west > other.east or 
                    self.north < other.south or 
                    self.south > other.north)


class Region():
    def __init__(self, ymin, xmin, nx, ny, ystep, xstep, name):
        self.south = ymin
        self.north = ymin + (ny * ystep)
        self.west = xmin
        self.east = xmin + (nx * xstep)
        self.name = name
    def __str__(self): return self.name


regions = [
        # region parameters as defined in RDWPS docs
        #   https://weather.gc.ca/grib/grib2_RDWPS_e.html
        Region(46.2590, -92.3116, 658, 318, 0.0090, 0.0124, 'superior'),
        Region(41.4260, -88.1452, 698, 573, 0.0090, 0.0124, 'huron-michigan'),
        Region(41.2190, -83.6068, 398, 210, 0.0090, 0.0124, 'erie'),
        Region(43.0640, -79.9736, 348, 158, 0.0090, 0.0124, 'ontario'),
        Region(44.0750, -70.9250, 331, 160, 0.0500, 0.0500, 'gulf-st-lawrence')
    ]


def ll_2_regionstr(south, north, west, east):
    """ convert input bounds to region string using the separating axis theorem """
    bounds = Boundary(south, north, west, east)
    return [str(reg) for reg in regions if bounds == reg]


def fetchname(wavevar, time, region):
    hour = f"{((time.hour % 24) // 6 * 6):02d}"  # better option?
    predictionhour = '000'  # any better than 0-hour?
    if 'gulf' not in region: 
        regionstr = 'lake-' + region
        grid = 'latlon0.0090x0.0124'
    else:
        regionstr = region
        grid = 'latlon0.05x0.05'
    return f"CMC_rdwps_{regionstr}_{wavevar}_{level_ref[region][wavevar]}_{grid}_{time.strftime('%Y%m%d')}{hour}_P{predictionhour}.grib2"


def fetch_rdwps(wavevar, time, regions):
    """
    fetchfiles = []
    for region in regions:
        fname = fetchname(wavevar, time, region)
        fetchfile = f"{storage_cfg()}{fname}"
        fetchfiles.append(fetchfile)
        if os.path.isfile(fetchfile):
            print(f"File {fname} exists, skipping retrieval...")
        else:
   """
    filenames = []
    for region in regions:
        fname = fetchname(wavevar, time, region)
        fetchfile = f"{storage_cfg()}{fname}"
        directory = 'great_lakes'
        if "gulf-st-lawrence" in fname:
            fetchurl = f"http://dd.weather.gc.ca/model_wave/ocean/{region}/grib2/{((time.hour % 24) // 6 * 6):02d}/{fname}"
        else:
            fetchurl = f"http://dd.weather.gc.ca/model_wave/great_lakes/{region}/grib2/{((time.hour % 24) // 6 * 6):02d}/{fname}"
        print(f"Downloading {fname} from the Regional Deterministic Wave Prediction System...")
        urllib.request.urlretrieve(fetchurl, fetchfile)
        filenames.append(fetchfile)

    return filenames


def load_rdwps(wavevar, time, regions, plot):
    filenames = []
    val = np.array([])
    lat = np.array([])
    lon = np.array([])
    for region in regions:
        fname = fetchname(wavevar, time, region)
        fetchfile = f"{storage_cfg()}{fname}"
        filenames.append(fetchfile)
        if not os.path.isfile(fetchfile):
            fetch_rdwps(wavevar, time, region)

        grib = pygrib.open(fetchfile)
        if plot is not False: fetch_util.plot_sample_grib(grib[1], plot)

        for msg in grib:
            val = np.append(val, msg.data()[0])
            lat = np.append(lat, msg.data()[1])
            lon = np.append(lon, msg.data()[2])
    """
    for msg in grib:
        lat, lon = np.array(msg.latlons())
        vals, lat, lon = msg.data()
        mask = vals.mask
        data = vals.data
        print(lat.shape)
        break
    """

    #grib[1].data()


    return val, lat, lon


class Rdwps(): 
    def fetch_windwaveswellheight(self, south=-90, north=90, west=-180, east=180, time=datetime.now()): 
        return fetch_rdwps('HTSGW', time, ll_2_regionstr(south, north, west, east))
    def fetch_windwaveheight(self, south=-90, north=90, west=-180, east=180, time=datetime.now()):
        return fetch_rdwps('WVHGT', time, ll_2_regionstr(south, north, west, east))
    def fetch_wavedirection(self, south=-90, north=90, west=-180, east=180, time=datetime.now()):
        return fetch_rdwps('WVDIR', time, ll_2_regionstr(south, north, west, east))
    def fetch_waveperiod(self, south=-90, north=90, west=-180, east=180, time=datetime.now()):
        return fetch_rdwps('WVPER', time, ll_2_regionstr(south, north, west, east))
    def fetch_wind_u(self, south=-90, north=90, west=-180, east=180, time=datetime.now()):
        return fetch_rdwps('UGRD', time, ll_2_regionstr(south, north, west, east))
    def fetch_wind_v(self, south=-90, north=90, west=-180, east=180, time=datetime.now()):
        return fetch_rdwps('VGRD', time, ll_2_regionstr(south, north, west, east))
    def fetch_icecover(self, south=-90, north=90, west=-180, east=180, time=datetime.now()): 
        return fetch_rdwps('ICEC', time, ll_2_regionstr(south, north, west, east))

    def load_windwaveswellheight(self, south=-90, north=90, west=-180, east=180, time=datetime.now(), plot=False):
        return load_rdwps('HTSGW', time, ll_2_regionstr(south, north, west, east), plot=plot)
    def load_windwaveheight(self, south=-90, north=90, west=-180, east=180, time=datetime.now(), plot=False):
        return load_rdwps('WVHGT', time, ll_2_regionstr(south, north, west, east), plot=plot)
    def load_wavedirection(self, south=-90, north=90, west=-180, east=180, time=datetime.now(), plot=False):
        return load_rdwps('WVDIR', time, ll_2_regionstr(south, north, west, east), plot=plot)
    def load_waveperiod(self, south=-90, north=90, west=-180, east=180, time=datetime.now(), plot=False):
        return load_rdwps('WVPER', time, ll_2_regionstr(south, north, west, east), plot=plot)
    def load_wind_u(self, south=-90, north=90, west=-180, east=180, time=datetime.now(), plot=False):
        return load_rdwps('UGRD', time, ll_2_regionstr(south, north, west, east), plot=plot)
    def load_wind_v(self, south=-90, north=90, west=-180, east=180, time=datetime.now(), plot=False):
        return load_rdwps('VGRD', time, ll_2_regionstr(south, north, west, east), plot=plot)
    def load_icecover(self, south=-90, north=90, west=-180, east=180, time=datetime.now(), plot=False):
        return load_rdwps('ICEC', time, ll_2_regionstr(south, north, west, east), plot=plot)

    def __str__(self):
        info = "RDWPS info goes here"
        args = "(south=-90, north=90, west=-180, east=180, time=datetime.now())"
        return fetch_util.str_def(self, info, args)


"""
# mahone bay test area:
    south =  44.4
    north =  44.7
    west  = -64.4
    east  = -63.8

time = datetime.now() - timedelta(hours=3)
"""


