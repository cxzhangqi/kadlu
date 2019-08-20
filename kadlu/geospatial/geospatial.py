""" Geospatial utilities module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import numpy as np
from osgeo import gdal
import scipy.io as sio


def read_matlab(path, name):
    """ Read data from a MatLab file.

        Args: 
            path: str
                File path
            name: str
                Name of MatLab field containing the data values   

        Returns:
            values: numpy array
                Data values
    """
    # load data
    m = sio.loadmat(path)

    # access array
    values = np.array(m[name])

    return values


def read_geotiff(path, band_id=1):
    """ Read data from a GeoTIFF file.

        Args: 
            path: str
                File path
            band_id: int
                Number of the GeoTIFF raster band containing the data values   

        Returns:
            values: masked numpy array
                Data values. The array has been masked where invalid values occur (NaNs or infs).
    """
    # load data
    data_set = gdal.Open(path)

    # access values
    band = data_set.GetRasterBand(band_id)
    values = data_set.ReadAsArray()

    # replace no-data value with nan
    nodata = band.GetNoDataValue()
    values[values == nodata] = np.nan
    values = np.ma.masked_invalid(values)

    return values


def crop(lat, lon, south, north, west, east, grid=False):
    """ Select rectangular region bounded by the geographical coordinates 
        latlon_SW to the south-west and latlon_NE to the north-east.

        If grid is False, lat and lon must have the same length.

        Args: 
            lat: 1d or 2d numpy array
                Latitude values
            lon: 1d or 2d numpy array
                Longitude values
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (NE) boundary of the region of interest.
            grid: bool
                Specify how to combine elements of lat and lon.

        Returns:
            ind: numpy array
                Selected indices. 1d if grid is False, 2d if grid is True.
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
    """
    ind_lat = np.argwhere((lat >= south) & (lat <= north))
    ind_lat = np.squeeze(ind_lat)

    ind_lon = np.argwhere((lon >= west) & (lon <= east))
    ind_lon = np.squeeze(ind_lon)

    if np.ndim(ind_lat) == 2:
        latc = ind_lat[:,0] + 1j * ind_lat[:,1]
        lonc = ind_lon[:,0] + 1j * ind_lon[:,1]
        x = np.intersect1d(latc, lonc)
        xr = np.real(x)
        xi = np.imag(x)
        xr = xr.astype(int)
        xi = xi.astype(int)
        ind_lat = np.arange(np.min(xr), np.max(xr)+1)
        ind_lon = np.arange(np.min(xi), np.max(xi)+1)
        lat = lat[:,0]
        lon = lon[0,:]
        ind = np.ix_(ind_lat, ind_lon)       

    if grid:
        ind = np.ix_(ind_lat, ind_lon)       
    else:
        ind = np.intersect1d(ind_lat, ind_lon)
        ind_lat = ind
        ind_lon = ind

    lat = lat[ind_lat]
    lon = lon[ind_lon]

    return ind, lat, lon