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


def fetch(storage_location, south=-90, north=90, west=-180, east=180):
    """ Fetch GEBCO bathymetry data.

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
            path: list
                Paths to the data files that were retrieved.
    """
    path = storage_location
    return path



def load(storage_location, south=-90, north=90, west=-180, east=180):
    """ Load GEBCO bathymetry data within specified geographical region.

        TODO: Get rid of the path argument and instead use the config.ini file

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
            bathy: 2d numpy array
                Bathymetry values
            lats: 1d numpy array
                Latitude values
            lons: 1d numpy array
                Longitude values
    """
    # fetch relevant data files
    path = fetch(storage_location, south=-90, north=90, west=-180, east=180)

    # load data
    d = Dataset(path)

    # access lat-lon arrays
    lats = np.array(d.variables["lat"])        
    lons = np.array(d.variables["lon"])     

    # crop
    indices, lats, lons = crop(lats, lons, south, north, west, east, grid=True)

    # access bathymetry values
    bathy = np.array(d.variables["bathy"])[indices]

    return (bathy,lats,lons)