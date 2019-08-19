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
from kadlu.geospatial.read import read_geotiff_2d
from kadlu.geospatial.bathy_reader import LatLon


def fetch(latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
    """ Fetch Non-Navigational NONNA-100 bathymetric data from the The 
        Canadian Hydrographic Service (CHS).

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.

        Returns:
            paths: list
                Paths to the data files that were retrieved.
    """

    # Select relevant data files covering the specified geographic region

    fnames = ["CA2_4300N06000W.tif",  "CA2_4400N06000W.tif"]
    paths = "/home/oliskir/src/meridian/kadlu/kadlu/tests/assets/tif" + fnames

    return paths


def read(path):
    """ Read bathymetry values from the data file.

        Args: 
            path: str
                File name

        Returns:
            val: 1d numpy array
                Data values
    """
    v, _, _ = read_geotiff_2d(path=path)
    return v


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
                File name

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