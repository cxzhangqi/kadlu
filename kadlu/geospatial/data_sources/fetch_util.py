import os
from os import path
from os.path import dirname
import numpy as np
import configparser
import warnings
import pygrib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import sqlite3
from datetime import datetime, timedelta


# database tables
hycom_tables = ['salinity', 'water_temp', 'water_u', 'water_v']
chs_table    = 'bathy'
era5_tables  = ['significant_height_of_combined_wind_waves_and_swell', 'mean_wave_direction', 'mean_wave_period']
wwiii_tables = ['hs', 'dp', 'tp']


def storage_cfg():
    """ return filepath containing storage configuration string

    first tries to check the config.ini file in kadlu root folder, if there's a 
    problem defaults to kadlu/storage and issues a warning
    """

    def default_storage(msg):
        """ helper function for storage_cfg() """
        storage_location = (path.abspath(path.dirname(dirname(dirname(dirname(__file__))))) + "/storage/")
        if not os.path.isdir(storage_location):
            os.mkdir(storage_location)
        warnings.warn(f"{msg} storage location will be set to {storage_location}")
        return storage_location

    cfg = configparser.ConfigParser()       # read .ini into dictionary object
    cfg.read(path.join(path.dirname(dirname(dirname(dirname(__file__)))), "config.ini"))
    try:
        storage_location = cfg["storage"]["StorageLocation"]
    except KeyError:                        # missing config.ini file
        return default_storage("missing kadlu/config.ini.")

    if storage_location is '':              # null value in config.ini
        return default_storage("null value in kadlu/config.ini.")

    if not path.isdir(storage_location):    # verify the location exists
        return default_storage("storage location doesn't exist.")

    return storage_location


def database_cfg():
    """ configure and connect to sqlite database """
    """
        time is stored as an INT in the database, where each integer
        is epoch hours since 2000-01-01 00:00
    """
    conn = sqlite3.connect(storage_cfg() + "geospatial.db")
    db = conn.cursor()

    # bathymetry table (CHS)
    db.execute(f"CREATE TABLE IF NOT EXISTS {chs_table}"
               "(   val     REAL,   "
               "    lat     REAL,   "
               "    lon     REAL,   "
               "    source  TEXT   )") 
    db.execute(f"CREATE UNIQUE INDEX IF NOT EXISTS "
               f"idx_{chs_table} on {chs_table}(lon, lat)")

    # hycom environmental data tables
    for fetchvar in hycom_tables:
        db.execute(f"CREATE TABLE IF NOT EXISTS {fetchvar}"
                    "(  val     INT,    "
                    "   lat     REAL,   "
                    "   lon     REAL,   "
                    "   time    INT,    "
                    "   depth   INT,    "
                    "   source  TEXT   )")
        db.execute(f"CREATE UNIQUE INDEX IF NOT EXISTS "
                   f"idx_{fetchvar} on {fetchvar}(time, lon, lat)")

    # wave data tables
    for fetchvar in era5_tables + wwiii_tables:
        db.execute(f"CREATE TABLE IF NOT EXISTS {fetchvar}"
                    "(  val     INT,    "
                    "   lat     REAL,   "
                    "   lon     REAL,   "
                    "   time    INT,    "
                    "   source  TEXT   )")
        db.execute(f"CREATE UNIQUE INDEX IF NOT EXISTS "
                   f"idx_{fetchvar} on {fetchvar}(time, lon, lat)")

    return conn, db


def dt_2_epoch(dt_arr):
    """ convert datetimes to epoch hours """
    t0 = datetime(2000, 1, 1, 0, 0, 0)
    delta = lambda dt : (dt - t0).total_seconds()/60/60
    dt_arr = np.array([dt_arr]) if np.array([dt_arr]).shape == (1,) else dt_arr
    return list(map(int, map(delta, dt_arr)))


def epoch_2_dt(ep_arr):
    """ convert epoch hours to datetimes """
    t0 = datetime(2000, 1, 1)
    return list(map(lambda ep : t0 + timedelta(hours=ep), ep_arr))


def index(val, sorted_arr):
    """ converts value in coordinate array to grid index """
    if val > sorted_arr[-1]: return len(sorted_arr) - 1
    return np.nonzero(sorted_arr >= val)[0][0]
  

class Boundary():
    """ compute intersecting boundaries using the separating axis theorem """

    def __init__(self, south, north, west, east, fetchvar=''):
        self.south, self.north, self.west, self.east, self.fetchvar = south, north, west, east, fetchvar

    def __str__(self): return self.fetchvar

    def intersects(self, other):  # separating axis theorem
        return not (self.east  < other.west or 
                    self.west  > other.east or 
                    self.north < other.south or 
                    self.south > other.north )


def ll_2_regionstr(south, north, west, east, regions, default=[]):
    """ convert input bounds to region strings using Boundary() class and separating axis theorem """

    if west > east:  # recursive function call if query intersects antimeridian
        return np.union1d(ll_2_regionstr(south, north, west,  180, regions, default), 
                          ll_2_regionstr(south, north, -180, east, regions, default))

    query = Boundary(south, north, west, east)
    matching = [str(reg) for reg in regions if query.intersects(reg)]

    if len(matching) == 0: 
        warnings.warn(f"No regions matched for query. Defaulting to {default} ({len(default)} regions)")
        return default

    return np.unique(matching)


def str_def(self, info, args):
    """ builds string definition for data source class objects """
    fcns = [fcn for fcn in dir(self) if callable(getattr(self, fcn)) and not fcn.startswith("__")]
    strlen = list(map(lambda f : len(f), fcns))
    whitespace = ''.join(map(lambda f : ' ', range(0, np.max(strlen) - np.min(strlen))))
    return f"{info}\n\nClass functions:\n\t" + "\n\t".join(map(lambda f : f"{f}{whitespace[len(f)-np.min(strlen):]}{args}", fcns ))


def plot_sample_grib(gribfiles, title_text="A sample plot"):

    #fig, axs = plt.subplots(2, int(len(gribfiles)/2))
    #
    #for x in range(0, len(gribfiles)):
    #    grb = pygrib.open(gribfiles[x])[1]
    #    ax = axs[x] if len(axs.shape) == 1 else axs
    #    if len(axs.shape) >= 2: ax = axs[int(x/axs.shape[0])][x%axs.shape[0]]

    for f in gribfiles:
        fig, ax = plt.subplots(1, 1)
        grb = pygrib.open(f)[1]

        data = grb.values
        lat, lon = grb.latlons()
        m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), resolution='l', ax=ax)
        x, y = m(lon,lat)

        # Paint map with parameter values under projected coordinates.
        #cs = ax.pcolormesh(x,y,data,shading='flat',cmap=plt.cm.jet)
        cs = ax.contourf(x, y, data, cmap=plt.cm.jet)

        # map filigree
        m.drawcoastlines()
        m.fillcontinents()
        m.drawmapboundary()
        m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0])
        m.drawmeridians(np.arange(-180.,180.,10.),labels=[0,0,0,1])

        # Plot legend, title.
        fig.colorbar(cs,orientation='vertical', ax=ax)
        #ax.set_title(title_text)
        plt.title(title_text)

        # Show plot.
        fig.tight_layout()
        plt.show()


def plot_coverage(lat, lon):
    fig = plt.figure()
    #m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), resolution='c')
    m=Basemap(projection='mill',lat_ts=10,
            llcrnrlon=-180, urcrnrlon=180,
            llcrnrlat=-90,urcrnrlat=90, 
            resolution='c')
    x, y = m(lon,lat)
    m.drawcoastlines()
    m.fillcontinents()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,10.),labels=[0,0,0,1])
    plt.scatter(x, y, 1, marker='.', color='xkcd:ocean blue', zorder=10)
    fig.tight_layout()
    plt.show()

"""
print(f"Range: lat {min(lat)}->{max(lat)}\tlon {min(lon)}->{max(lon)}")
plot_coverage(lat, lon)
"""

