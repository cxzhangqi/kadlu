import os
from os import path
from os.path import dirname
import numpy as np
import configparser
import warnings
import pygrib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def init_default_storage_dir(msg):
    storage_location = (path.abspath(dirname(dirname(dirname(dirname(__file__))))) + "/storage/")
    if not os.path.isdir(storage_location):
        os.mkdir(storage_location)
    warnings.warn("%s storage location will be set to %s" % (msg, storage_location))
    return storage_location


def instantiate_storage_config():
    cfg = configparser.ConfigParser()
    cfg.read(path.join(dirname(dirname(dirname(dirname(__file__)))), "config.ini"))
    try:
        storage_location = cfg["storage"]["StorageLocation"]
    except KeyError:  # missing config.ini file
        return init_default_storage_dir("missing kadlu/config.ini.")

    if storage_location is '':  # null value in config.ini
        return init_default_storage_dir("null value in kadlu/config.ini.")

    if not path.isdir(storage_location):  # verify the location exists
        return init_default_storage_dir("storage location doesn't exist.")

    return storage_location


def validate_wavesource(filepath, sourceDict):
    # this may become obsolete...
    # seems kind of redundant
    # matt_s 2019-08
    gribdata = pygrib.fromstring(open(filepath, "rb").read())
    try:
        assert(gribdata['shortName'] in sourceDict.keys() or gribdata['shortName'] in sourceDict.values())
    except AssertionError as err:
        print("Specified source file is not a wave source")
        raise


def plotSampleGrib(grb, title_text):
    plt.figure()
    data = grb.values
    lat, lon = grb.latlons()

    m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), resolution='l')
    x, y = m(lon,lat)

    # Paint map with parameter values under projected coordinates.
    cs = m.pcolormesh(x,y,data,shading='flat',cmap=plt.cm.jet)

    # map filigree
    m.drawcoastlines()
    m.fillcontinents()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,10.),labels=[0,0,0,1])

    # Plot legend, title.
    plt.colorbar(cs,orientation='vertical')
    plt.title(title_text)

    # Show plot.
    plt.show()

