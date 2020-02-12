"""
    API for Regional Deterministic Wave Prediction System (RDWPS)
    Data provided by Govt. of Canada

    Metadata regarding the dataset can be found here:
        https://weather.gc.ca/grib/grib2_RDWPS_e.html
"""

import numpy as np
from datetime import datetime, timedelta
import os
#import urllib.request
import requests
import pygrib
from kadlu.geospatial.data_sources import data_util
from kadlu.geospatial.data_sources.data_util import storage_cfg, Boundary, ll_2_regionstr
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


class Region():
    """ utility class for computing intersecting regions. compatible with ll_2_regionstr() in fetch_util module """
    def __init__(self, ymin, xmin, nx, ny, ystep, xstep, name):
        self.south = ymin
        self.north = ymin + (ny * ystep)
        self.west = xmin
        self.east = xmin + (nx * xstep)
        self.name = name

    def __str__(self): return self.name


rdwps_regions = [
        #region boundary parameters as defined in RDWPS docs:
        #    https://weather.gc.ca/grib/grib2_RDWPS_e.html
        Region(46.2590, -92.3116, 658, 318, 0.0090, 0.0124, 'superior'),
        Region(41.4260, -88.1452, 698, 573, 0.0090, 0.0124, 'huron-michigan'),
        Region(41.2190, -83.6068, 398, 210, 0.0090, 0.0124, 'erie'),
        Region(43.0640, -79.9736, 348, 158, 0.0090, 0.0124, 'ontario'),
        Region(44.0750, -70.9250, 331, 160, 0.0500, 0.0500, 'gulf-st-lawrence')
    ]


def fetchname(wavevar, time, region):
    """ generate a filename for given wave variable, time, and region """
    timestr = datetime.now().strftime('%Y%m%d')
    hour = f"{((datetime.now()-timedelta(hours=3)).hour // 6) * 6:02d}"
    predictionhour = f"{(time.hour // 3) * 3:03d}"

    if 'gulf' not in region: 
        regionstr = 'lake-' + region
        grid = 'latlon0.0090x0.0124'
    else:
        regionstr = region
        grid = 'latlon0.05x0.05'
    return (f"CMC_rdwps_{regionstr}_{wavevar}_{level_ref[region][wavevar]}"
            f"_{grid}_{timestr}{hour}_P{predictionhour}.grib2")


def fetch_rdwps(wavevar, start, end, regions):
    """ download rdwps data and return associated filepaths 

    args:
        wavevar: string
            the variable short name of desired wave parameter according to RDWPS docs
            the complete list of variable short names can be found here
            https://weather.gc.ca/grib/grib2_RDWPS_e.html
        start: datetime
            the start of the desired time range
        end: datetime
            the end of the desired time range
        regions: list
            list of strings containing regions to fetch

    return:
        filenames: list
            list of strings containing complete file paths of fetched data
    """
    filenames = []
    time = datetime(start.year, start.month, start.day, (start.hour // 3 * 3))

    # RDWPS is a prediction service - requested times must be in the 
    # following 48 hours from the current time
    assert(time >= datetime.now() - timedelta(hours=6))
    assert(end <= datetime.now() + timedelta(hours=48))
    assert(time <= end)
    
    while time <= end:
        for reg in regions:
            fname = fetchname(wavevar, time, reg)
            fetchfile = f"{storage_cfg()}{fname}"
            directory = 'ocean' if 'gulf-st-lawrence' in fname else 'great_lakes'
            hour = f"{(((datetime.now()-timedelta(hours=3)).hour) // 6 * 6):02d}"
            fetchurl = f"http://dd.weather.gc.ca/model_wave/{directory}/{reg}/grib2/{hour}/{fname}"
            print(f"Downloading {fname} from the Regional Deterministic Wave Prediction System...")
            #urllib.request.urlretrieve(fetchurl, fetchfile)
            grib = requests.get(fetchurl)
            assert(grib.status_code == 200)
            with open(fetchfile, 'wb') as f: f.write(grib.content)
            filenames.append(fetchfile)

        time += timedelta(hours=3)

    return filenames


def load_rdwps(wavevar, start, end, south, north, west, east, plot):
    """ return downloaded rdwps data for specified wavevar according to given time, lat, lon boundaries

    args:
        wavevar: string
            the variable short name of desired wave parameter according to RDWPS docs
            the complete list of variable short names can be found here
            https://weather.gc.ca/grib/grib2_RDWPS_e.html
        start: datetime
            the start of the desired time range
        end: datetime
            the end of the desired time range
        south, north: float
            ymin, ymax coordinate boundaries (latitude). range: -90, 90
        west, east: float
            xmin, xmax coordinate boundaries (longitude). range: -180, 180
        plot: boolean
            if true a plot will be output (experimental feature)

    return:
        filenames: list
            list of strings containing complete file paths of fetched data
    """
    filenames = []
    val = np.array([])
    lat = np.array([])
    lon = np.array([])
    timestamps = np.array([])
    time = datetime(start.year, start.month, start.day, (start.hour // 3 * 3))
    regions = ll_2_regionstr(south, north, west, east, rdwps_regions)

    # RDWPS is a prediction service - requested times must be in the 
    # following 48 hours from the current time
    assert(time >= datetime.now() - timedelta(hours=6))
    assert(end <= datetime.now() + timedelta(hours=48))
    assert(time <= end)

    while time <= end:
        for reg in regions:
            # get the filename
            fname = fetchname(wavevar, time, reg)
            fetchfile = f"{storage_cfg()}{fname}"
            filenames.append(fetchfile)
            if not os.path.isfile(fetchfile): fetch_rdwps(wavevar, start, end, regions)

            # open the file and get the raw data
            grib = pygrib.open(fetchfile)
            assert(grib.messages == 1)
            msg = grib[1]
            z_grid, y_grid, x_grid = msg.data()

            for slx in range(z_grid.shape[0]):
                z = z_grid[slx]
                y = y_grid[slx]
                x = x_grid[slx]
                x -= 360  # normalize longitudes

                # build index to collect points in area of interest
                latix = np.array([l >= south and l <= north for l in y])
                lonix = np.array([l >= west and l <= east for l in x])
                ix = latix & lonix

                # append points within AoI to return arrays
                val = np.append(val, z[ix])
                lat = np.append(lat, y[ix])
                lon = np.append(lon, x[ix])
                timestamps = np.append(timestamps, [time for x in range(sum(ix))])

        time += timedelta(hours=3)

    if plot is not False: fetch_util.plot_sample_grib(filenames, plot)
    
    return val, lat, lon, timestamps


class Rdwps(): 
    """ collection of module functions for fetching and loading. abstracted to include a seperate function for each variable """

    def fetch_windwaveswellheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now()): 
        return fetch_rdwps('HTSGW', start, end, ll_2_regionstr(south, north, west, east, rdwps_regions))

    def fetch_windwaveheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now()):
        return fetch_rdwps('WVHGT', start, end, ll_2_regionstr(south, north, west, east, rdwps_regions))
    
    def fetch_wavedirection(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now()):
        return fetch_rdwps('WVDIR', start, end, ll_2_regionstr(south, north, west, east, rdwps_regions))
    
    def fetch_waveperiod(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now()):
        return fetch_rdwps('WVPER', start, end, ll_2_regionstr(south, north, west, east, rdwps_regions))
    
    def fetch_wind_u(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now()):
        return fetch_rdwps('UGRD', start, end, ll_2_regionstr(south, north, west, east, rdwps_regions))
    
    def fetch_wind_v(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now()):
        return fetch_rdwps('VGRD', start, end, ll_2_regionstr(south, north, west, east, rdwps_regions))
    
#    def fetch_icecover(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now()): 
#        return fetch_rdwps('ICEC', start, end, ll_2_regionstr(south, north, west, east, rdwps_regions))

    def load_windwaveswellheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now(), plot=False):
        return load_rdwps('HTSGW', start, end, south, north, west, east, plot=plot)
    
    def load_windwaveheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now(), plot=False):
        return load_rdwps('WVHGT', start, end, south, north, west, east, plot=plot)
    
    def load_wavedirection(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now(), plot=False):
        return load_rdwps('WVDIR', start, end, south, north, west, east, plot=plot)
    
    def load_waveperiod(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now(), plot=False):
        return load_rdwps('WVPER', start, end, south, north, west, east, plot=plot)
    
    def load_wind_u(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now(), plot=False):
        return load_rdwps('UGRD', start, end, south, north, west, east, plot=plot)
    
    def load_wind_v(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now(), plot=False):
        return load_rdwps('VGRD', start, end, south, north, west, east, plot=plot)
    
#    def load_icecover(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=3), end=datetime.now(), plot=False):
#        return load_rdwps('ICEC', start, end, south, north, west, east, plot=plot)

    def __str__(self):
        info = '\n'.join([
                "API for Regional Deterministic Wave Prediction System (RDWPS)",
                "Data provided by Govt. of Canada.",
                "This dataset provides a prediction service - it is necessary to use start",
                "and end times within the next 48 hours to return RDWPS predictions.",
                "Metadata regarding the dataset can be found here:",
                "\thttps://weather.gc.ca/grib/grib2_RDWPS_e.html"
            ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(), end=datetime())"
        return fetch_util.str_def(self, info, args)

