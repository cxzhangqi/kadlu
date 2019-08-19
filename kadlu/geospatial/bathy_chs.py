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
from kadlu.geospatial.read import read_geotiff_2d
from kadlu.geospatial.bathy_reader import LatLon


def fetch(storage_location, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
    """ Fetch Non-Navigational NONNA-100 bathymetric data from the The 
        Canadian Hydrographic Service (CHS).

        TODO: Get rid of the storage_location argument and instead use the config.ini file
        TODO: Implement fetching

        Args: 
            storage_location: str
                Path to where data files are stored in the local file system.
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.

        Returns:
            paths: list
                Paths to the data files that were retrieved.
    """
    # select relevant files
    fnames = select_files(latlon_SW, latlon_NE)
    paths = list()
    for fname in fnames:
        paths.append(os.path.join(storage_location, fname))

    # attempt to fetch those files we do not already have
    for i,path in enumerate(paths):
        exists = os.path.exists(path)
        # if not exists:
            # ... implement fetch part here ...

    # check again
    fetched = list()
    for i,path in enumerate(paths):
        exists = os.path.exists(path)
        if exists:
            fetched.append(path)

    if len(fetched) < len(paths):
        print("Only fetched {0} of {1} maps necessary to fully cover the specified region".format(len(fetched), len(paths)))

    return fetched


def select_files(latlon_SW, latlon_NE):
    """ Select the bathymetry data files that overlap with a specific 
        geographic region.

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.

        Returns:
            fnames: list
                Names of the bathymetry data files.
    """
    # lat range
    lat_min = int(np.floor(latlon_SW.latitude))
    lat_max = int(np.floor(latlon_NE.latitude))
    if lat_max == latlon_NE.latitude:
        lat_max -= 1

    # lon range
    lon_min = int(np.floor(latlon_SW.longitude))
    lon_max = int(np.floor(latlon_NE.longitude))
    if lon_max == latlon_NE.longitude:
        lon_max -= 1

    # create lat,lon arrays
    lats = np.arange(start=lat_min, stop=lat_max+1, step=1)
    lons = np.arange(start=lon_min, stop=lon_max+1, step=1)

    # create list of filenames
    fnames = list()
    for lat in lats:
        for lon in lons:
            sw_corner = LatLon(lat, lon)
            fname = filename(sw_corner)
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
    z, _, _ = read_geotiff_2d(path=path)
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
    sw_corner = parse_sw_corner(path)

    # latitude step size in degrees
    dlat = 0.001

    # longitude step size in degrees
    if sw_corner.latitude < 68:
        dlon = 0.001
    elif sw_corner.latitude >=68 and sw_corner.latitude < 80:
        dlon = 0.002
    elif sw_corner.latitude >= 80:
        dlon = 0.004

    lats = np.arange(num_lat, dtype=np.float)
    lats *= dlat
    lats += sw_corner.latitude

    lons = np.arange(num_lon, dtype=np.float)
    lons *= dlon
    lons += sw_corner.longitude

    return lats, lons


def parse_sw_corner(path):
    """ Parse latitude and longitude data from filename.

        Args: 
            path: str
                Path to data file

        Returns:
            sw_corner: LatLon
                Latitude and longitude of the SW corner of the data set
    """
    fname = os.path.basename(path)

    north = int(fname[4:8]) / 100
    west = int(fname[9:14]) / 100
    east = -west

    assert east >= -180 and east <= 180, 'Invalid parsed longitude value'
    assert north >= -90 and north <= 180, 'Invalid parsed latitude value'

    sw_corner = LatLon(north, east)

    return sw_corner


def filename(sw_corner):
    """ Construct filename from latitude and longitude of SW corner.

        Args: 
            sw_corner: LatLon
                Latitude and longitude of the SW corner of the data set

        Returns:
            fname: str
                File name
    """
    north = int(sw_corner.latitude * 100)
    east = int(sw_corner.longitude * 100)
    west = -east
    fname = "CA2_{0:04d}N{1:05d}W.tif".format(north, west)
    return fname