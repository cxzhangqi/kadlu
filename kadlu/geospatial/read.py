""" Geospatial data reader module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import numpy as np
from netCDF4 import Dataset
from osgeo import gdal
import scipy.io as sio


def read_netcdf(path, val_name, lat_name=None, lon_name=None, depth_name=None):
    """ Read geospatial data from a NetCDF file.

        The data may be on two-dimensional (lat,lon) grid or a three-dimensional 
        (lat,lon,depth) grid.

        Args: 
            path: str
                File path
            val_name: str
                Name of NetCDF variable containing the data values   
            lat_name: str
                Name of NetCDF variable containing the latitude values 
            lon_name: str
                Name of NetCDF variable containing the longitude values
            depth_name: str
                Name of NetCDF variable containing the depth values     

        Returns:
            val: numpy array
                Data values
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
            depth: numpy array
                Depth values. None, if data is on a 2d (lat,lon) grid
    """
    # load data
    d = Dataset(path)

    lat, lon, depth = None, None, None

    # access arrays
    val = np.array(d.variables[val_name])

    if lat_name is not None:
        lat = np.array(d.variables[lat_name])

    if lon_name is not None:
        lon = np.array(d.variables[lon_name])

    if depth_name is not None:
        depth = np.array(d.variables[depth_name])

    return val, lat, lon, depth


def read_netcdf_2d(path, val_name, lat_name=None, lon_name=None):
    """ Read geospatial data on a two-dimensional (lat,lon) grid from a NetCDF file.

        Args: 
            path: str
                File path
            val_name: str
                Name of NetCDF variable containing the data values   
            lat_name: str
                Name of NetCDF variable containing the latitude values
            lon_name: str
                Name of NetCDF variable containing the longitude values

        Returns:
            val: numpy array
                Data values
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
    """
    val, lat, lon, _ = read_netcdf(path, val_name, lat_name, lon_name)
    return val, lat, lon


def read_matlab(path, val_name, lat_name=None, lon_name=None, depth_name=None):
    """ Read geospatial data from a NetCDF file.

        The data may be on two-dimensional (lat,lon) grid or a three-dimensional 
        (lat,lon,depth) grid.

        Args: 
            path: str
                File path
            val_name: str
                Name of NetCDF variable containing the data values   
            lat_name: str
                Name of NetCDF variable containing the latitude values
            lon_name: str
                Name of NetCDF variable containing the longitude values
            depth_name: str
                Name of NetCDF variable containing the depth values

        Returns:
            val: numpy array
                Data values
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
            depth: numpy array
                Depth values. None, if data is on a 2d (lat,lon) grid
    """
    # load data
    m = sio.loadmat(path)

    lat, lon, depth = None, None, None

    # access arrays
    val = np.array(m[val_name])

    if lat_name is not None:
        lat = np.squeeze(np.array(m[lat_name]))

    if lon_name is not None:
        lon = np.squeeze(np.array(m[lon_name]))

    if depth_name is not None:
        depth = np.squeeze(np.array(m[depth_name]))

    return val, lat, lon, depth


def read_matlab_2d(path, val_name, lat_name=None, lon_name=None):
    """ Read geospatial data on a two-dimensional (lat,lon) grid from a MATLAB file.

        Args: 
            path: str
                File path
            val_name: str
                Name of NetCDF variable containing the data values   
            lat_name: str
                Name of NetCDF variable containing the latitude values
            lon_name: str
                Name of NetCDF variable containing the longitude values

        Returns:
            val: numpy array
                Data values
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
    """
    val, lat, lon, _ = read_matlab(path, val_name, lat_name, lon_name)
    return val, lat, lon


def read_geotiff(path, val_id=1, lat_id=None, lon_id=None, depth_id=None):
    """ Read geospatial data from a GeoTIFF file.

        The data may be on two-dimensional (lat,lon) grid or a three-dimensional 
        (lat,lon,depth) grid.

        Args: 
            path: str
                File path
            val_id: int
                Number of the GeoTIFF raster band containing the data values   
            lat_id: int
                Number of the GeoTIFF raster band containing the latitude values
            lon_id: int
                Number of the GeoTIFF raster band containing the longitude values
            depth_id: int
                Number of the GeoTIFF raster band containing the depth values

        Returns:
            val: masked numpy array
                Data values. The array has been masked where invalid values occur (NaNs or infs).
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
            depth: numpy array
                Depth values. None, if data is on a 2d (lat,lon) grid
    """
    # load data
    data_set = gdal.Open(path)

    lat, lon, depth = None, None, None

    # access arrays
    band = data_set.GetRasterBand(val_id)
    val = data_set.ReadAsArray()

    # replace no-data value with nan
    nodata = band.GetNoDataValue()
    val[val == nodata] = np.nan
    val = np.ma.masked_invalid(val)

    if lat_id is not None:
        band = data_set.GetRasterBand(lat_id)
        lat = data_set.ReadAsArray()

    if lon_id is not None:
        band = data_set.GetRasterBand(lon_id)
        lon = data_set.ReadAsArray()

    if depth_id is not None:
        band = data_set.GetRasterBand(depth_id)
        depth = data_set.ReadAsArray()

    return val, lat, lon, depth


def read_geotiff_2d(path, val_id=1, lat_id=None, lon_id=None):
    """ Read geospatial data on a two-dimensional (lat,lon) grid from a GeoTIFF file.

        Args: 
            path: str
                File path
            val_id: int
                Number of the GeoTIFF raster band containing the data values   
            lat_id: int
                Number of the GeoTIFF raster band containing the latitude values
            lon_id: int
                Number of the GeoTIFF raster band containing the longitude values

        Returns:
            val: numpy array
                Data values
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
    """
    val, lat, lon, _ = read_geotiff(path, val_id, lat_id, lon_id)
    return val, lat, lon
