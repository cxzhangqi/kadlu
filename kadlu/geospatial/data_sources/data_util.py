import os, sys
from os import path
from os.path import dirname
import numpy as np
import configparser
import warnings
import sqlite3
from datetime import datetime, timedelta
import json
import pickle
from functools import reduce
from hashlib import md5
from contextlib import contextmanager, redirect_stdout, redirect_stderr

#import pygrib
#import matplotlib.pyplot as plt


# database tables for data fetching and loading
chs_table    = 'bathy'
hycom_tables = ['salinity', 'water_temp', 'water_u', 'water_v']
wwiii_tables = ['hs', 'dp', 'tp', 'windU', 'windV']
era5_tables  = [
        'significant_height_of_combined_wind_waves_and_swell',
        'mean_wave_direction', 
        'mean_wave_period', 
        'u_component_of_wind', 
        'v_component_of_wind'
    ]

cfg = configparser.ConfigParser()       # read .ini into dictionary object
cfg.read(path.join(path.dirname(dirname(dirname(dirname(__file__)))), "config.ini"))


def storage_cfg():
    """ return filepath containing storage configuration string

    first checks the config.ini file in kadlu root folder, if there's a 
    problem defaults to kadlu/storage and issues a warning
    """

    def default_storage(msg):
        """ helper function for storage_cfg() """
        storage_location = (path.abspath(path.dirname(dirname(dirname(dirname(__file__))))) + "/storage/")
        if not os.path.isdir(storage_location):
            os.mkdir(storage_location)
        print(f"NOTICE: {msg} storage location will be set to {storage_location}")
        return storage_location

    try:
        storage_location = cfg["storage"]["storage_location"]
    except KeyError:                        # missing config.ini file
        return default_storage("missing kadlu/config.ini.")

    if storage_location == '':              # null value in config.ini
        return default_storage("null value in kadlu/config.ini.")

    if not path.isdir(storage_location):    # verify the location exists
        return default_storage("storage location doesn't exist.")

    return storage_location


def database_cfg():
    """ configure and connect to sqlite database

        time is stored as an integer in the database, where each value
        is epoch hours since 2000-01-01 00:00

        returns:
            conn:   
                database connection object
            db:
                connection cursor object
    """
    conn = sqlite3.connect(storage_cfg() + "geospatial.db")
    db = conn.cursor()

    # bathymetry table (CHS)
    db.execute(f"CREATE TABLE IF NOT EXISTS {chs_table}  "
               "(   val     REAL    NOT NULL, "
               "    lat     REAL    NOT NULL, "
               "    lon     REAL    NOT NULL, "
               "    source  TEXT    NOT NULL )") 
    db.execute(f"CREATE UNIQUE INDEX IF NOT EXISTS "
               f"idx_{chs_table} on {chs_table}(lon, lat, val, source)")

    # hycom environmental data tables
    for var in hycom_tables:
        db.execute(f"CREATE TABLE IF NOT EXISTS {var}"
                    "( val     REAL NOT NULL, "
                    "  lat     REAL NOT NULL, "
                    "  lon     REAL NOT NULL, "
                    "  time    INT  NOT NULL, "
                    "  depth   INT  NOT NULL, "
                    "  source  TEXT NOT NULL )")
        db.execute(f"CREATE UNIQUE INDEX IF NOT EXISTS "
                   f"idx_{var} on {var}(time, lon, lat, depth, val, source)")

    # wave data tables
    for var in era5_tables + wwiii_tables:
        db.execute(f"CREATE TABLE IF NOT EXISTS {var}"
                    "( val     REAL    NOT NULL, "
                    "  lat     REAL    NOT NULL, "
                    "  lon     REAL    NOT NULL, "
                    "  time    INT     NOT NULL, "
                    "  source  TEXT    NOT NULL )") 
        db.execute(f"CREATE UNIQUE INDEX IF NOT EXISTS "
                   f"idx_{var} on {var}(time, lon, lat, val, source)")

    db.execute("CREATE TABLE IF NOT EXISTS fetch_map"
                '(  hash    INT  NOT NULL,  '
                '   bytes   BLOB           )' )
    db.execute(f"CREATE UNIQUE INDEX IF NOT EXISTS "
                 f"idx_fetched on fetch_map(hash)")

    return conn, db


def bin_db():
    """ database for storing serialized objects in memory """
    conn = sqlite3.connect('file::memory:?cache=shared', uri=True)
    db = conn.cursor()
    db.execute('CREATE TABLE IF NOT EXISTS bin'
                '(  hash    INT  NOT NULL,  '
                '   bytes   BLOB NOT NULL  )' )
    db.execute(f"CREATE UNIQUE INDEX IF NOT EXISTS "
                 f"idx_bin on bin(hash)")

    raise ResourceWarning('fcn not used')
    #return conn, db


def hash_key(kwargs, seed, block=('lock',)):
    """ compute unique hash and convert to 8-byte int as serialization key """
    qry = kwargs.copy()
    for x in block:
        if x in qry.keys(): del qry[x]
    string = seed + json.dumps(qry, sort_keys=True, default=str)
    key = int(md5(string.encode('utf-8')).hexdigest(), base=16)
    return key >> 80  # bitshift value by 80 bits: SQLite max value is 64 bits


def serialized(kwargs, seed=''):
    """ returns true if fetch query hash exists in database """
    conn, db = database_cfg()
    key = hash_key(kwargs, seed)
    if 'lock' in kwargs.keys(): kwargs['lock'].acquire()
    db.execute('SELECT * FROM fetch_map WHERE hash == ? LIMIT 1', (key,))
    res = db.fetchone()
    if 'lock' in kwargs.keys(): kwargs['lock'].release()
    if res is None: return False
    if res[1] is not None: return res[1]
    return True


def insert_hash(kwargs, seed='', obj=None):
    """ create hash index in database to record query history 

        optionally include an object to be serialized and cached
    """
    qry = kwargs.copy()
    conn, db = database_cfg()
    if 'lock' in qry.keys(): del qry['lock']
    key = hash_key(qry, seed)
    db.execute('INSERT OR IGNORE INTO fetch_map VALUES (?,?)',
               (key, pickle.dumps(obj)))
    conn.commit()
    return


def deserialize(kwargs, persisting=True, seed=''):
    """ read binary from the database and load it as python object """
    conn, db = bin_db()
    key = hash_key(kwargs, seed)
    db.execute('SELECT * FROM fetch_map WHERE hash == ?', (key,))
    res = db.fetchone()
    if res is None: raise KeyError('no data found for query')
    if not persisting:
        db.execute('DELETE FROM fetch_map WHERE hash == ?', (key,))
        conn.commit()
    return pickle.loads(res[1])


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


def flatten(cols, frames):
    """ dimensional reduction by taking average of time frames """
    # assert that frames are of equal size
    assert reduce(lambda a, b: (a==b)*a, frames[1:] - frames[:-1])

    ix = range(len(frames) -1)
    fsplit= np.array([cols[3][frames[f] : frames[f +1]] for f in ix])
    vals = (reduce(np.add, fsplit) / len(fsplit))

    if len(cols) == 4:
        _, y, x, _ = cols[:, frames[0] : frames[1]]
    else:
        _, y, x, _, z = cols[:, frames[0] : frames[1]]
        return vals, y, x, z


def reshape_2D(cols):
    return dict(values=cols[0], lats=cols[1], lons=cols[2])


def reshape_3D(cols):
    """ prepare loaded data for interpolation """
    if isinstance(cols[0], (float, int)):
        return dict(values=cols[0])
    frames = np.append(np.nonzero(cols[3][1:] > cols[3][:-1])[0] + 1, len(cols[3]))
    if len(np.unique(frames)) > 1: vals, y, x, z = flatten(cols, frames) 
    else: vals, y, x, _, z  = cols
    rows = np.array((vals, y, x, z)).T

    # reshape row data to 3D array
    xgrid, ygrid, zgrid = np.unique(x), np.unique(y), np.unique(z)
    gridspace = np.full((len(ygrid), len(xgrid), len(zgrid)), fill_value=-30000)
    # this could potentially be optimized to avoid an index lookup cost
    for row in rows:
        x_ix = index(row[2], xgrid)
        y_ix = index(row[1], ygrid)
        z_ix = index(row[3], zgrid)
        gridspace[y_ix, x_ix, z_ix] = row[0]

    # remove -30000 values for interpolation:
    # fill missing depth values with last value in each column
    # this section could be cleaned up
    for xi in range(0, gridspace.shape[0]):
        for yi in range(0, gridspace.shape[1]):
            col = gridspace[xi, yi]
            if sum(col == -30000) > 0 and sum(col == -30000) < len(col):
                col[col == -30000] = col[col != -30000][-1]
                gridspace[xi, yi] = col

    # TODO:
    # create default values for columns that are entirely null

    return dict(values=gridspace, lats=ygrid, lons=xgrid, depths=zgrid)


def str_def(self, info, args):
    """ builds string definition for data source class objects """
    fcns = [fcn for fcn in dir(self) if callable(getattr(self, fcn)) and not fcn.startswith("__")]
    strlen = list(map(lambda f : len(f), fcns))
    whitespace = ''.join(map(lambda f : ' ', range(0, np.max(strlen) - np.min(strlen))))
    return f"{info}\n\nClass functions:\n\t" + "\n\t".join(map(lambda f : f"{f}{whitespace[len(f)-np.min(strlen):]}{args}", fcns ))


@contextmanager
def dev_null():
    """ context manager to redirect output to /dev/null """
    with open(os.devnull, 'w') as null:
        try:
            with redirect_stderr(null) as err, redirect_stdout(null) as out: 
                yield (err, out)
        finally:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__


class Boundary():
    """ compute intersecting boundaries with separating axis theorem """
    def __init__(self, south, north, west, east, fetchvar=''):
        self.south, self.north, self.west, self.east, self.fetchvar = \
                south, north, west, east, fetchvar

    def __str__(self): return self.fetchvar

    def intersects(self, other):  # separating axis theorem
        return not (self.east  < other.west or
                    self.west  > other.east or
                    self.north < other.south or
                    self.south > other.north)


def ll_2_regionstr(south, north, west, east, regions, default=[]):
    """ convert input bounds to region strings using Boundary class """

    if west > east:  # recursive function call if query intersects antimeridian
        return np.union1d(ll_2_regionstr(south, north, west,  180, regions, default), 
                          ll_2_regionstr(south, north, -180, east, regions, default))

    query = Boundary(south, north, west, east)
    matching = [str(reg) for reg in regions if query.intersects(reg)]

    if len(matching) == 0: 
        warnings.warn(f"No regions matched for query. Defaulting to {default} ({len(default)} regions)")
        return default

    return np.unique(matching)


"""
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

"""
"""

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

#def gen_kwargs():
#    """ some sample fetch/load keyword args for rapid testing """
#    """
#    kwargs = gen_kwargs()
#    self = Ocean(**kwargs)
#    """
#    return dict(
#        start=datetime(2015, 1, 9), end=datetime(2015, 1, 9, 3),
#        #time=datetime(2015, 1, 9),
#        south=44,                   west=-64.5, 
#        north=46,                   east=-62.5, 
#        top=0,                      bottom=5000
#    )

