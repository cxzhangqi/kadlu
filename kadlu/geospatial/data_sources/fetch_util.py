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
    for fname in filenames:
        print(f"fname: {fname}")
        grib = pygrib.open(fname)
        if plot is not False: plot_sample_grib(grib[1], plot)

    #data = [None for msg in range(grib.messages)]
    #for x in range(0, grib.messages):
    #    data[x] = grib[x+1].data()  # grib indexing starts at 1

    #return np.array(data)

    """
    for msg in grib:
        lat, lon = np.array(msg.latlons())
        vals, lat, lon = msg.data()
        mask = vals.mask
        data = vals.data
        print(lat.shape)
        break
    """

    return grib[1].data()


def plot_sample_grib(grb, title_text="A sample plot"):
    data = grb.values
    lat, lon = grb.latlons()

    plt.figure()
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

def print_fcn_names(source):
    for fcn in source.fetch_functions:
        print(f"{fcn.__name__}{source.header}")


