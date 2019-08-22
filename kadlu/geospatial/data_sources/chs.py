""" Module within the kadlu package for handling Non-Navigational NONNA-100 
    bathymetric data from the The Canadian Hydrographic Service (CHS). 

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import os
import numpy as np
from kadlu.geospatial.geospatial import crop, read_geotiff


def fetch(storage_location, south=-90, north=90, west=-180, east=180):
    """ Fetch Non-Navigational NONNA-100 bathymetric data from the The 
        Canadian Hydrographic Service (CHS).

        TODO: Get rid of the storage_location argument and instead use the config.ini file
        TODO: Implement fetching

        Args: 
            storage_location: str
                Path to where data files are stored in the local file system.
            south: float
                Southern boundary of the region of interest.
            north: float
                Northern boundary of the region of interest.
            west: float
                Western boundary of the region of interest.
            east: float
                Eastern boundary of the region of interest.

        Returns:
            paths: list
                Paths to the data files that were retrieved.
    """
    # select relevant files
    fnames = select_files(south, north, west, east)
    paths = list()
    for fname in fnames:
        paths.append(os.path.join(storage_location, fname))

    # attempt to fetch those files we do not already have
    for path in paths:
        exists = os.path.exists(path)
        # if not exists:
            # ... implement fetch part here ...

    # check again
    fetched = list()
    for path in paths:
        exists = os.path.exists(path)
        if exists:
            fetched.append(path)

    if len(fetched) < len(paths):
        print("Only fetched {0} of {1} maps necessary to fully cover the specified region".format(len(fetched), len(paths)))

    return fetched


def load(storage_location, south=-90, north=90, west=-180, east=180):
    """ Load Non-Navigational NONNA-100 bathymetric data from the Canadian Hydrographic 
        Service (CHS) within specified geographical region.

        TODO: Get rid of the storage_location argument and instead use the config.ini file

        Args: 
            south: float
                Southern boundary of the region of interest.
            north: float
                Northern boundary of the region of interest.
            west: float
                Western boundary of the region of interest.
            east: float
                Eastern boundary of the region of interest.

        Returns:
            bathy: 1d numpy array
                Bathymetry values
            lats: 1d numpy array
                Latitude values
            lons: 1d numpy array
                Longitude values
    """
    # fetch relevant data files
    files = fetch(storage_location, south, north, west, east)

    bathy, lats, lons = list(), list(), list()        

    # loop over geotiff files
    for f in files:

        # load data from a single file
        z,y,x = load_from_file(f) 

        # crop the region of interest
        indices, y, x = crop(y, x, south, north, west, east)
        z = z[indices]

        # collect data
        bathy.append(z)
        lats.append(y)
        lons.append(x)

    # concatenate
    bathy = np.ma.concatenate(bathy)
    lats = np.concatenate(lats)
    lons = np.concatenate(lons)

    return (bathy,lats,lons)


def load_from_file(path):
    """ Load bathymetric data from a GeoTIFF file provided by the Canadian Hydrographic 
        Service (CHS) as part of the Non-Navigational NONNA-100 bathymetric data series.

        Args: 
            path: str
                Path to the GeoTIFF file.

        Returns:
            z: 1d numpy array
                Bathymetry values
            y: 1d numpy array
                Latitude values
            x: 1d numpy array
                Longitude values
    """

    # read data from geotiff file
    z = read(path)

    # create lat-lon arrays
    lat, lon = latlon(path)

    # make a grid
    x, y = np.meshgrid(lon, lat)

    # select non-masked entries
    x = x[~z.mask]
    y = y[~z.mask]
    z = z[~z.mask]

    return z,y,x


def select_files(south, north, west, east):
    """ Select the bathymetry data files that overlap with a specific 
        geographic region.

        Args: 
            south: float
                Southern boundary of the region of interest.
            north: float
                Northern boundary of the region of interest.
            west: float
                Western boundary of the region of interest.
            east: float
                Eastern boundary of the region of interest.

        Returns:
            fnames: list
                Names of the bathymetry data files.
    """
    # lat range
    lat_min = int(np.floor(south))
    lat_max = int(np.floor(north))
    if lat_max == north:
        lat_max -= 1

    # lon range
    lon_min = int(np.floor(west))
    lon_max = int(np.floor(east))
    if lon_max == east:
        lon_max -= 1

    # create lat,lon arrays
    lats = np.arange(start=lat_min, stop=lat_max+1, step=1)
    lons = np.arange(start=lon_min, stop=lon_max+1, step=1)

    # create list of filenames
    fnames = list()
    for lat in lats:
        for lon in lons:
            fname = filename(lat, lon)
            fnames.append(fname)

    return fnames


def read(path):
    """ Read bathymetry values from the data file.

        Args: 
            path: str
                File name

        Returns:
            val: 1d numpy array
                Data values
    """
    z = read_geotiff(path=path)
    z = np.flip(z, axis=0)
    return z


def latlon(path, num_lat=1001, num_lon=1001):
    """ Create latitude and longitude arrays for a CHS bathymetry data file.

        The number of latitude and longitude values is 1001 by default.

        The latitude binning is 0.001 degrees.

        The longitude binning depends on the latitude, as follows

            * south of 68 deg N: 0.001 deg
            * latitudes from 68 to 80 deg N: 0.002 deg
            * 80 deg N and north: 0.004 deg

        Args: 
            path: str
                File name
            num_lat: int
                Number of latitude grid points
            num_lon: int
                Number of longitude grid points

        Returns:
            lats: numpy array
                Uniformly spaced latitude values 
            lons: numpy array
                Uniformly spaced longitude values 
    """
    # parse SW corner of the map
    south, west = parse_sw_corner(path)

    # latitude step size in degrees
    dlat = 0.001

    # longitude step size in degrees
    if south < 68:
        dlon = 0.001
    elif south >=68 and south < 80:
        dlon = 0.002
    elif south >= 80:
        dlon = 0.004

    lats = np.arange(num_lat, dtype=np.float)
    lats *= dlat
    lats += south

    lons = np.arange(num_lon, dtype=np.float)
    lons *= dlon
    lons += west

    return lats, lons


def parse_sw_corner(path):
    """ Parse latitude and longitude data from filename.

        Args: 
            path: str
                Path to data file

        Returns:
            south: float
                Southern boundary of map
            west: float
                Western boundary of map
    """
    fname = os.path.basename(path)

    south = int(fname[4:8]) / 100
    west = -int(fname[9:14]) / 100

    assert west >= -180 and west <= 180, 'Invalid parsed longitude value'
    assert south >= -90 and south <= 90, 'Invalid parsed latitude value'

    return south, west


def filename(south, west):
    """ Construct filename from latitude and longitude of SW corner.

        Args: 
            south: float
                Southern boundary of map
            west: float
                Western boundary of map

        Returns:
            fname: 
                File name
    """
    fname = "CA2_{0:04d}N{1:05d}W.tif".format(int(south * 100), -int(west * 100))
    return fname