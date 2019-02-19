""" Bathymetry reader module within the kadlu package

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
from kadlu.util import get_files


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
    GEOTIFF_CHS = 4
    XYZ = 5

def get_member(cls, member_name):
    for name, member in cls.__members__.items():
        if member_name == name:
            return member

    s = ", ".join(name for name, _ in cls.__members__.items())
    raise ValueError("Unknown value \'{0}\'. Select between: {1}".format(member_name, s))


class BathyReader():
    """ Class for reading bathymetry data from NetCDF files (*.nc) or MATLAB files (*.mat).

        Attributes: 
            input: str
                Normally, this should be the full path to a single data file.
                However, in the special case where format='GEOTIFF_CHS', this 
                could also be the name of a folder containing multiple *.tif files.
            lat_name: str
                Name of the variable that contains the latitue values.
            lon_name: str
                Name of the variable that contains the longitude values.
            bathy_name: str
                Name of the variable that contains the bathymetry values.
            format: str
                Format that the bathymetry data is stored in.
                Options are: NETCDF, MATLAB, GEOTIFF, GEOTIFF_CHS, XYZ
                If no format is provided (default), the BathyReader will 
                attempt to determine the format from the input path.
    """
    def __init__(self, input, lat_name='lat', lon_name='lon', bathy_name='bathy', lon_axis=1, format=None):

        self.input = input

        if format is None:
            format = self._determine_format_from_filename(input)

        # file format
        self.format = get_member(Format, format)

        self.lat_name = lat_name
        self.lon_name = lon_name
        self.bathy_name = bathy_name
        self.lon_axis = lon_axis

    def _determine_format_from_filename(self, input):

        fmt = None

        f = input
        ext = f[f.rfind('.'):]

        if ext == '.nc':
            fmt = 'NETCDF'
        elif ext == '.mat':
            fmt = 'MATLAB'
        elif ext == '.tif':
            fmt = 'GEOTIFF'
        elif ext == '.xyz':
            fmt = 'XYZ'

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
        elif self.format is Format.GEOTIFF_CHS:
            _r = self._read_geotiff_chs
        elif self.format is Format.XYZ:
            _r = self._read_xyz

        # read
        lat, lon, bathy = _r(latlon_SW, latlon_NE)

        # extra checks for gridded data
        lat, lon, bathy = self._check_grid(lat, lon, bathy)

        return lat, lon, bathy

    def _check_grid(self, lat, lon, bathy):

        if np.ndim(bathy) == 1:
            return lat, lon, bathy

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
        return None, None, None

    def _read_geotiff_chs(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):

        band_no = 1

        # check if file or directory
        if os.path.isdir(self.input):
            files = get_files(path=self.input, substr='.tif')
        else:
            files = [self.input]

        # empty lists
        lat, lon, bathy = list(), list(), list()

        # loop over files
        for f in files:

            # parse south-west corner from file name
            sw = self._chs_parse_sw_corner(f)

            # generate lat-lon arrays
            y, x = self._chs_generate_latlon_arrays(sw)

            # check if file overlaps with requested region
            overlap = np.min(y) <= latlon_NE.latitude \
                and np.max(y) >= latlon_SW.latitude \
                and np.min(x) <= latlon_NE.longitude \
                and np.max(x) >= latlon_SW.longitude

            if not overlap:
                continue

            # open file
            data_set = gdal.Open(f)

            # get bathy matrix
            band = data_set.GetRasterBand(band_no)
            z = data_set.ReadAsArray()

            # get nodata value from the GDAL band object
            nodata = band.GetNoDataValue()

            # flip bathy matrix
            z = np.flip(z, axis=0)

            # empty lists
            _x, _y, _z = list(), list(), list()

            # loop over bathy matrix, keeping only those entries that 
            # have data and are inside requested region
            N = z.shape[0]
            for i in range(N):
                for j in range(N):
                    if z[i,j] != nodata and \
                        x[i] >= latlon_SW.longitude and \
                        x[i] <= latlon_NE.longitude and \
                        y[j] >= latlon_SW.latitude and \
                        y[j] <= latlon_NE.latitude:

                        _x.append(x[i])
                        _y.append(y[j])
                        _z.append(z[i,j])

            # fill lists
            lat = lat + _y
            lon = lon + _x
            bathy = bathy + _z

        return np.array(lat), np.array(lon), np.array(bathy)

    def _chs_parse_sw_corner(self, path):
        f = path[path.rfind('/')+1:]
        north = int(f[4:8]) / 100
        west = int(f[9:14]) / 100
        east = west - 180
        return LatLon(north, east)

    def _chs_generate_latlon_arrays(self, sw):
        N = 1001
        dlat = 0.001
        if sw.longitude < 68:
            dlon = 0.001
        elif sw.longitude >=68 and sw.longitude < 80:
            dlon = 0.002
        elif sw.longitude >= 80:
            dlon = 0.004
        else:
            print('Request longitude is outside tabulated range') 

        lats = np.arange(N, dtype=np.float)
        lats *= dlat
        lats += sw.latitude

        lons = np.arange(N, dtype=np.float)
        lons *= dlon
        lons += sw.longitude

        return lats, lons

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
