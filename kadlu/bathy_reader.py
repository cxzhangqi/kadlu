""" Bathymetry reader module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import numpy as np
from collections import namedtuple
from netCDF4 import Dataset
from osgeo import gdal
import scipy.io as sio
from enum import Enum


LatLon = namedtuple('LatLon', ['latitude', 'longitude'])
""" Latitude and longitude coordinates for a given location.

    Args: 
        latitude: float
            Latitude in degrees from -90 (South Pole) to +90 (North Pole).
        longitude: float
            Longitude in degrees from -180 to +180 (West to East) with 0 
            corresponding to the Greenwich Meridian.
"""


class Format(Enum):
    NETCDF = 1
    MATLAB = 2
    GEOTIFF = 3
    XYZ = 4


class BathyReader():
    """ Class for reading bathymetry data from NetCDF files (*.nc) or MATLAB files (*.mat).

        Attributes: 
            input: str or list(str)
                File name or list of file names. 
            lat_name: str
                Name of the variable that contains the latitue values.
            lon_name: str
                Name of the variable that contains the longitude values.
            bathy_name: str
                Name of the variable that contains the bathymetry values.
            format: str
                Format that the bathymetry data is stored in.
                Options are: NETCDF, MATLAB, GEOTIFF, XYZ
                If no format is provided (default), the BathyReader will 
                attempt to determine the format from the input path.
    """
    def __init__(self, input, lat_name='lat', lon_name='lon', bathy_name='bathy', lon_axis=1, format=None):

        self._check_input(input)

        if format is None:
            self.format = self._determine_format(input)

        self.input = input
        self.lat_name = lat_name
        self.lon_name = lon_name
        self.bathy_name = bathy_name
        self.lon_axis = lon_axis

    def _check_input(self, input):
        d = np.ndim(input)
        l = len(input)
        assert d==0 or (d==1 and l>0), "Input should be a string or a (non-empty) list of strings"

    def _determine_format(self, input):

        d = np.ndim(input)

        if d == 0: 
            f = input
        else:
            f = input[0]

        ext = f[f.rfind('.'):]

        if ext == '.nc':
            fmt = Format.NETCDF
        elif ext == '.mat':
            fmt = Format.MATLAB
        elif ext == '.tif':
            fmt = Format.GEOTIFF
        elif ext == '.xyz':
            fmt = Format.XYZ
        else:
            print('Unrecognized file format')
            exit(1)

        return fmt


    def read(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry from file.

            Args: 
                latlon_SW: LatLon
                    South-western (SW) boundary of the region of interest.
                latlon_NE: LatLon
                    North-eastern (SE) boundary of the region of interest.

            Returns:
                lat: 1d numpy array
                    Latitude values
                lon: 1d numpy array
                    Longitude values
                bathy: 2d numpy array
                    Bathymetry values
        """
        # select reader
        if self.format is Format.NETCDF:
            _r = self._read_netcdf
        elif self.format is Format.MATLAB:
            _r = self._read_matlab
        elif self.format is Format.GEOTIFF:
            _r = self._read_geotiff
        elif self.format is Format.XYZ:
            _r = self._read_xyz

        # read
        lat, lon, bathy = _r(latlon_SW, latlon_NE)

        # ensure that lat and lon are strictly increasing
        if np.all(np.diff(lat) < 0):
            lat = np.flip(lat, axis=0)
            bathy = np.flip(bathy, axis=0)
        if np.all(np.diff(lon) < 0):
            lon = np.flip(lat, axis=0)
            bathy = np.flip(bathy, axis=1)
        
        assert np.all(np.diff(lat) > 0), 'Latitudes must be strictly ascending'
        assert np.all(np.diff(lon) > 0), 'Longitudes must be strictly ascending'

        # ensure axis=0 is latitude and axis=1 is longitude
        if self.lon_axis == 0:
            bathy = np.swapaxes(bathy, 0, 1)

        # check that axes have consistent sizes
        if bathy.shape[0] != lat.shape[0]:
            bathy = np.swapaxes(bathy, 0, 1)

        assert bathy.shape[0] == lat.shape[0], 'axis #0 of bathymetry matrix must have same size as latitude array'
        assert bathy.shape[1] == lon.shape[0], 'axis #1 of bathymetry matrix must have same size as longitude array'

        return lat, lon, bathy

    def _read_netcdf(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry matrices from file.

            Args: 
                latlon_SW: LatLon
                    South-western (SW) boundary of the region of interest.
                latlon_NE: LatLon
                    North-eastern (NE) boundary of the region of interest.

            Returns:
                lat: 1d numpy array
                    Latitude values
                lon: 1d numpy array
                    Longitude values
                bathy: 2d numpy array
                    Bathymetry values
        """
        # load data
        d = Dataset(self.input)
        lat = np.array(d.variables[self.lat_name])
        lon = np.array(d.variables[self.lon_name])
        bathy = np.array(d.variables[self.bathy_name])

        # select region of interest
        lat, lon, bathy = self._select_region(lat, lon, bathy, latlon_SW, latlon_NE)

        return lat, lon, bathy

    def _read_matlab(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry matrices from file.

            Args: 
                latlon_NW: LatLon
                    South-western (SW) boundary of the region of interest.
                latlon_SE: LatLon
                    North-eastern (NE) boundary of the region of interest.

            Returns:
                lat: 1d numpy array
                    Latitude values
                lon: 1d numpy array
                    Longitude values
                bathy: 2d numpy array
                    Bathymetry values
        """
        # load data
        m = sio.loadmat(self.input)
        lat = np.squeeze(np.array(m[self.lat_name]))
        lon = np.squeeze(np.array(m[self.lon_name]))
        bathy = np.array(m[self.bathy_name])

        # select region of interest
        lat, lon, bathy = self._select_region(lat, lon, bathy, latlon_SW, latlon_NE)

        return lat, lon, bathy

    def _read_geotiff(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry matrices from file.

            Args: 
                latlon_NW: LatLon
                    South-western (SW) boundary of the region of interest.
                latlon_SE: LatLon
                    North-eastern (NE) boundary of the region of interest.

            Returns:
                lat: 1d numpy array
                    Latitude values
                lon: 1d numpy array
                    Longitude values
                bathy: 2d numpy array
                    Bathymetry values
        """
        # load data
        m = sio.loadmat(self.input)
        lat = np.squeeze(np.array(m[self.lat_name]))
        lon = np.squeeze(np.array(m[self.lon_name]))
        bathy = np.array(m[self.bathy_name])

        # select region of interest
        lat, lon, bathy = self._select_region(lat, lon, bathy, latlon_SW, latlon_NE)

        return lat, lon, bathy

    def _read_xyz(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        return None, None, None

    def _select_region(self, lat, lon, bathy, latlon_SW, latlon_NE):
        """ Select rectangular region bounded by the geographical coordinates 
            latlon_SW to the south-west and latlon_NE to the north-east.

            Latitudes and longitudes can be given as 2d numpy arrays, but
            must be in regularly spaced grid.

            Args: 
                lat: 1d or 2d numpy array
                    Latitude values
                lon: 1d or 2d numpy array
                    Longitude values
                bathy: 2d numpy array
                    Bathymetry values
                latlon_SW: LatLon
                    South-western (SW) boundary of the region of interest.
                latlon_NE: LatLon
                    North-eastern (NE) boundary of the region of interest.

            Returns:
                lat: 1d numpy array
                    Latitude values
                lon: 1d numpy array
                    Longitude values
                bathy: 2d numpy array
                    Bathymetry values
        """
        ind_lat = np.argwhere((lat >= latlon_SW.latitude) & (lat <= latlon_NE.latitude))
        ind_lat = np.squeeze(ind_lat)

        ind_lon = np.argwhere((lon >= latlon_SW.longitude) & (lon <= latlon_NE.longitude))
        ind_lon = np.squeeze(ind_lon)

        if np.ndim(ind_lat) == 1:
            ind = np.ix_(ind_lat, ind_lon)       

        else:
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

        lat = lat[ind_lat]
        lon = lon[ind_lon]
        bathy = bathy[ind]

        return lat, lon, bathy
