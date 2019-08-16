""" Geospatial data reader module within the kadlu package

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
from collections import namedtuple
from netCDF4 import Dataset
from osgeo import gdal
import scipy.io as sio
from enum import Enum
from kadlu.utils import get_files


def read_netcdf(path, lat_name, lon_name, val_name, depth_name=None):
    """ Read geospatial data from a NetCDF file.

        The data may be on two-dimensional (lat,lon) grid or a three-dimensional 
        (lat,lon,depth) grid.

        Args: 
            path: str
                File path
            lat_name: str
                Name of NetCDF variable containing the latitude values
            lon_name: str
                Name of NetCDF variable containing the longitude values
            val_name: str
                Name of NetCDF variable containing the data values   
            depth_name: str
                Name of NetCDF variable containing the depth values (optional)       

        Returns:
            lat: 1d numpy array
                Latitude values
            lon: 1d numpy array
                Longitude values
            val: numpy array
                Data values
            depth: 1d numpy array
                Depth values. None, if data is on a 2d (lat,lon) grid
    """
    # load data
    d = Dataset(path)

    # access arrays
    lat = np.array(d.variables[lat_name])
    lon = np.array(d.variables[lon_name])
    val = np.array(d.variables[val_name])
    if depth_name is not None:
        depth = np.array(d.variables[depth_name])
    else:
        depth = None

    return lat, lon, val, depth


def read_netcdf_2d(path, lat_name, lon_name, val_name):
    """ Read geospatial data on a two-dimensional (lat,lon) grid from a NetCDF file.

        Args: 
            path: str
                File path
            lat_name: str
                Name of NetCDF variable containing the latitude values
            lon_name: str
                Name of NetCDF variable containing the longitude values
            val_name: str
                Name of NetCDF variable containing the data values   

        Returns:
            lat: 1d numpy array
                Latitude values
            lon: 1d numpy array
                Longitude values
            val: numpy array
                Data values
    """
    lat, lon, val, _ = read_netcdf(path, lat_name, lon_name, val_name)
    return lat, lon, val


def read_matlab(path, lat_name, lon_name, val_name, depth_name=None):
    """ Read geospatial data from a NetCDF file.

        The data may be on two-dimensional (lat,lon) grid or a three-dimensional 
        (lat,lon,depth) grid.

        Args: 
            path: str
                File path
            lat_name: str
                Name of NetCDF variable containing the latitude values
            lon_name: str
                Name of NetCDF variable containing the longitude values
            val_name: str
                Name of NetCDF variable containing the data values   
            depth_name: str
                Name of NetCDF variable containing the depth values (optional)       

        Returns:
            lat: 1d numpy array
                Latitude values
            lon: 1d numpy array
                Longitude values
            val: numpy array
                Data values
            depth: 1d numpy array
                Depth values. None, if data is on a 2d (lat,lon) grid
    """
    # load data
    m = sio.loadmat(path)

    # access arrays
    lat = np.squeeze(np.array(m[lat_name]))
    lon = np.squeeze(np.array(m[lon_name]))
    val = np.array(m[val_name])
    if depth_name is not None:
        depth = np.array(m[depth_name])
    else:
        depth = None

    return lat, lon, val, depth


def read_matlab_2d(path, lat_name, lon_name, val_name):
    """ Read geospatial data on a two-dimensional (lat,lon) grid from a MATLAB file.

        Args: 
            path: str
                File path
            lat_name: str
                Name of NetCDF variable containing the latitude values
            lon_name: str
                Name of NetCDF variable containing the longitude values
            val_name: str
                Name of NetCDF variable containing the data values   

        Returns:
            lat: 1d numpy array
                Latitude values
            lon: 1d numpy array
                Longitude values
            val: numpy array
                Data values
    """
    lat, lon, val, _ = read_matlab(path, lat_name, lon_name, val_name)
    return lat, lon, val
