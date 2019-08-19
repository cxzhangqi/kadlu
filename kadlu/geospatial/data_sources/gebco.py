""" Module within the kadlu package for handling GEBCO bathymetric data. 

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
from netCDF4 import Dataset
from kadlu.geospatial.geospatial import crop
from kadlu.geospatial.bathy_reader import LatLon


def fetch(storage_location, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
    """ Fetch GEBCO bathymetry data.

        TODO: Implement fetching

        Args: 
            storage_location: str
                Path to where data files are stored in the local file system.
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.

        Returns:
            path: list
                Paths to the data files that were retrieved.
    """
    path = storage_location
    return path



def load(storage_location, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
    """ Load GEBCO bathymetry data within specified geographical region.

        TODO: Get rid of the path argument and instead use the config.ini file

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.

        Returns:
            bathy: 2d numpy array
                Bathymetry values
            lats: 1d numpy array
                Latitude values
            lons: 1d numpy array
                Longitude values
    """
    # fetch relevant data files
    path = fetch(storage_location, latlon_SW, latlon_NE)

    # load data
    d = Dataset(path)

    # access lat-lon arrays
    lats = np.array(d.variables["lat"])        
    lons = np.array(d.variables["lon"])     

    # crop
    indices, lats, lons = crop(lats, lons, latlon_SW, latlon_NE, grid=True)

    # access bathymetry values
    bathy = np.array(d.variables["bathy"][indices])

    return bathy, lats, lons