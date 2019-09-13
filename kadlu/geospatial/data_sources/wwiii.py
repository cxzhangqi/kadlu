"""
    User guides:
        https://github.com/NOAA-EMC/WW3/wiki/WAVEWATCH-III-User-Guide
"""

import numpy as np
from datetime import datetime, timedelta
import os
from urllib.request import urlretrieve
import pygrib
from kadlu.geospatial.data_sources import fetch_util
from kadlu.geospatial.data_sources.fetch_util import storage_cfg

"""
waveSources = {
    'swh' : 'hs',
    'mwd' : 'dp',
    'mwp' : 'tp'
}

regions = {
    # This may be replaced by a more modular mapping system later
    # depending on how the API works for [region]_[interval] formatting
    'global'    : 'glo_30m',
    'arctic'    : 'ao_30m',
    'pacific'   : 'ep_10m',
    'atlantic'  : {'10m' : 'at_10m', '4m' : 'at_4m'},
    'US west'   : {'10m' : 'wc_10m', '4m' : 'wc_4m'},
    'alaska'    : {'10m' : 'ak_10m', '4m' : 'ak_4m'}
}

def fetch(wavevar=waveSources['swh'], time=datetime.now(), region=regions['global']):
    #    matt_s 2019-08
    #    Note that WWIII returns data from an entire month, regardless of
    #    the time given. In the future the load function should be updated
    #    to find the desired message within the grib file containing data for
    #    the desired date

    storage_location = fetch_util.instantiate_storage_config() 
    fname = fetchname(wavevar, time, region)
    fetchfile = f"{storage_location}{fname}"

    if os.path.isfile(fetchfile):
        print("File exists, skipping retrieval...")
        # obsolete ???
        #validate_wavesource(fetchfile, waveSources)
    else:
        print("Downloading file from NOAA WaveWatch III...")
        fetchurl = f"https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/{time.strftime('%Y/%m')}/gribs/{fname}"
        urllib.request.urlretrieve(fetchurl, fetchfile)
    return fetchfile
"""

wwiii_global = Region(, 'glo_30m'),  # global
wwiii_regions = [
        Region(, 'at_4m'),  # atlantic, also available : at_10m for 10 minute resolution
        Region(, 'wc_4m'),  # US west, also available : wc_10m for 10 minute resolution
        Region(, 'ak_4m'),  # alaska, also available : ak_10m for 10 minute resolution
        Region(, 'ao_30m'),  # arctic 
        Region(, 'ep_10m'),  # pacific
        wwiii_global
    ]


class Boundary():
    def __init__(self, south, north, west, east):
        self.south = south
        self.north = north
        self.west = west
        self.east = east

    def intersects(self, other):
        """
        separating axis theorem: quickly check for intersecting regions
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


def ll_2_regionstr(south, north, west, east):
    """ convert input bounds to region strings using the separating axis theorem """
    regions = [str(reg) for reg in wwiii_regions if Boundary(south, north, west, east).intersects(reg)]

    # default to low 30m resolution for global queries
    if (len(regions) == len(wwiii_regions)): return [str(wwiii_global)]
    return regions


def fetchname(wavevar, time, region): return f"multi_1.{region}.{wavevar}.{time.strftime('%Y%m')}.grb2"


def fetch_wwiii(wavevar, south, north, west, east, start, end):
    """
    reg = 'glo_30m'
    wavevar = 'hs'
    start = datetime(2017, 2, 3)
    end = datetime(2017, 2, 4)
    """
    regions = ll_2_regionstr(south, north, west, east)
    time = datetime(start.year, start.month, 1)
    filenames = []

    while time <= end:
        for reg in regions:
            fname = fetchname(wavevar, time, reg)
            fetchfile = f"{storage_cfg()}{fname}"
            print(f"Downloading {fname} from NOAA WaveWatch III...")
            fetchurl = f"https://data.nodc.noaa.gov/thredds/fileServer/ncep/nww3/{time.strftime('%Y/%m')}/gribs/{fname}"
            urlretrieve(fetchurl, fetchfile)
            filenames.append(fetchfile)

        # on this datasource, data is sorted per month
        # some months have more days than exactly 4 weeks
        time += timedelta(weeks=4)
        while (fetchname(wavevar, time, reg) == fname): time += timedelta(days=1)

    return filenames

def load_wwiii(wavevar, south, north, west, east, start, end, plot=False):
    """
    south=-90
    north=90
    west=-180
    east=180
    reg = 'glo_30m'
    """
    val = np.array([])
    lat = np.array([])
    lon = np.array([])
    timestamps = np.array([])
    regions = ll_2_regionstr(south, north, west, east)
    time = datetime(start.year, start.month, 1)

    while time <= end:
        for reg in regions:
            fname = fetchname(wavevar, time, reg)
            fetchfile = f"{storage_cfg()}{fname}"
            if not os.path.isfile(fetchfile): fetch_wwiii(wavevar, south, north, west, east, start, end)

            grib = pygrib.open(fetchfile)
            for msg in grib:
                msgtime = msg.validDate
                if msgtime < start : continue
                if msgtime > end : break
                z, y, x = msg.data()
                x -= 180  # normalize longitudes

                for slx in range(z.shape[0]):
                    latix = np.array([l >= south and l <= north for l in y[slx]])
                    lonix = np.array([l >= west and l <= east for l in x[slx]])
                    ix = latix & lonix

                    val = np.append(val, z[slx][ix])
                    lat = np.append(lat, y[slx][ix])
                    lon = np.append(lon, x[slx][ix])
                    timestamps = np.append(timestamps, [msgtime for x in range(sum(ix))])

        # on this datasource, data is sorted per month
        # some months have more days than exactly 4 weeks
        time += timedelta(weeks=4)
        while (fetchname(wavevar, time, reg) == fname): time += timedelta(days=1)

    return val, lat, lon, timestamps


class Wwiii():
    def fetch_windwaveheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_wwiii('hs', south, north, west, east, start, end)
    def fetch_wavedirection(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_wwiii('dp', south, north, west, east, start, end)
    def fetch_waveperiod(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return fetch_wwiii('tp', south, north, west, east, start, end)

    def load_windwaveswellheight(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_wwiii('hs', south, north, west, east, start, end)
    def load_wavedirection(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_wwiii('dp', south, north, west, east, start, end)
    def load_waveperiod(self, south=-90, north=90, west=-180, east=180, start=datetime.now()-timedelta(hours=24), end=datetime.now()):
        return load_wwiii('tp', south, north, west, east, start, end)

    def __str__(self):
        info = '\n'.join([ "Wavewatch info goes here" ])
        args = "(south=-90, north=90, west=-180, east=180, start=datetime(), end=datetime())"
        return fetch_util.str_def(self, info, args)

"""
    # If no gribfile argument is provided, default to the fetched file.
    if grib is None:
        grib = self.fetch_filename

    # load grib structure from target.
    grbs=pygrib.open(grib)

    # Identify parameter and date for extraction.
    # Date field: validDate
    # Target date (not incl. time yet)
    # DEBUG        date_valid = datetime(2018,11,2)
    if (target_date is None):
        target_date = self.fetch_datetimestamp.replace(minute=0, hour=0, second=0, microsecond=0)
    else:
    ### Should any filtering on time be added here to enforce valid time intervals?
        target_date = target_date

    # Fetch the indicated slice from the overall Grib file.
    grb = grbs.select(validDate=target_date,shortNameECMF=wavevar.value)[0]

    if(wavevar == WWIIIWavevar.hs):
        title_text = "WW III Sig. Wave + Swell Height from GRIB\n({}) ".format(target_date)
    elif(wavevar == WWIIIWavevar.tp):
        title_text = "WW III Mean Wave Period from GRIB\n({}) ".format(target_date)
    elif(wavevar == WWIIIWavevar.dp):
        title_text = "WW III Mean Wave Direction from GRIB\n({}) ".format(target_date)
    else:
        title_text = "WW III\nUnknown variable from GRIB\n({}) ".format(target_date)

    return (grb, title_text)
"""

