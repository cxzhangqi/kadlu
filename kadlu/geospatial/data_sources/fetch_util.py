import os
from os import path
from os.path import dirname
import numpy as np
import configparser
import warnings
import pygrib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def default_storage(msg):
    storage_location = (path.abspath(path.dirname(dirname(dirname(dirname(__file__))))) + "/storage/")
    if not os.path.isdir(storage_location):
        os.mkdir(storage_location)
    warnings.warn("%s storage location will be set to %s" % (msg, storage_location))
    return storage_location


def storage_cfg():
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


def str_def(self, info, args):
    """ builds string definition for data source class objects """
    fcns = [fcn for fcn in dir(self) if callable(getattr(self, fcn)) and not fcn.startswith("__")]
    strlen = list(map(lambda f : len(f), fcns))
    whitespace = ''.join(map(lambda f : ' ', range(0, np.max(strlen) - np.min(strlen))))
    return f"{info}\n\nClass functions:\n\t" + "\n\t".join(map(lambda f : f"{f}{whitespace[len(f)-np.min(strlen):]}{args}", fcns ))


def loadgrib(filenames, plot=False):
    """
    This needs to be updated to return the entire list of grib
    messages. This is because WWIII data source returns data for
    the entire month in 3-hour chunks. This will allow the 
    load() function within the wwiii module to parse out the
    desired 3-hour chunk from the list of messages. 

    see the wwiii testing scripts in the playroom folder

    matt_s 2019-08
    """
    if plot is not False:
        plot_sample_grib(filenames, plot)

    for fname in filenames:
        grib = pygrib.open(fname)

    #data = [None for msg in range(grib.messages)]
    #for x in range(0, grib.messages):
    #    data[x] = grib[x+1].data()  # grib indexing starts at 1

    #return np.array(data)

        for msg in grib:
            lat, lon = np.array(msg.latlons())
            vals, lat, lon = msg.data()
            mask = vals.mask
            data = vals.data
            break

    warnings.warn("This function is deprecated. Instead, use module function")
    return grib[1].data()


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

