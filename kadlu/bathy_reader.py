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
    """ Enum class for bathymetry file formats
    """
    NETCDF = 1
    MATLAB = 2
    GEOTIFF = 3
    XYZ = 4

def get_member(cls, member_name):
    """ Search for a member of an Enum-derived class by name.

        Returns ValueError if the class does not contain any member by that name.

        Args:
            cls: Class
                Class derived from Enum

            member_name: str
                Name of the member

        Returns:
            member: Class member    
                Class member with the specified name
    """
    for name, member in cls.__members__.items():
        if member_name == name:
            return member

    s = ", ".join(name for name, _ in cls.__members__.items())
    raise ValueError("Unknown value \'{0}\'. Select between: {1}".format(member_name, s))


class BathyReader():
    """ Class for reading bathymetry data.
    
        BathyReader can handle four different file formats:
            
            * NetCDF (*.nc) 
            * MATLAB (*.mat)
            * GeoTIFF (*.tif)
            * XYZ (*.xyz)

        BathyReader expects NetCDF and MATLAB files to contain bathymetry 
        data on a regular grid. Specically, it expects to find two 1d arrays 
        with uniformly spaced latitude and longitude coordinates, and a 2d 
        array with the bathymetry values at each grid point. The names of these 
        arrays can be specified by the user. See the list of attributes below.

        BathyReader expects GeoTIFF and XYZ files to conform to the standard 
        adopted by the Canadian Hydrographic Service (CHS) for their Non-Navigational 
        (NONNA-100) Bathymetric Data. In particular, the files must adhere to the 
        same naming scheme used by the CHS. For more information, see 

        https://open.canada.ca/data/en/dataset/d3881c4c-650d-4070-bf9b-1e00aabf0a1d

        Attributes: 
            input: str
                Normally, this should be the full path to a single data file.
                However, for GeoTIFF and XYZ formats this can also be the name of 
                a folder containing multiple *.tif or *.xyz files.
            lat_name: str
                Name of the variable that contains the latitue values. 
                Only relevant for NetCDF and MATLAB files.
            lon_name: str
                Name of the variable that contains the longitude values.
                Only relevant for NetCDF and MATLAB files.
            bathy_name: str
                Name of the variable that contains the bathymetry values.
                Only relevant for NetCDF and MATLAB files.
            lon_axis: int
                Used to specify which of the axes of the 2d bathymetry array 
                corresponds to the longitude axis.
                Only relevant for NetCDF and MATLAB files.
    """
    def __init__(self, input, lat_name='lat', lon_name='lon', bathy_name='bathy', lon_axis=1):

        # get files
        self.files = self._get_files(input)
        assert len(self.files) > 0, "You must provide at least one data file" 

        # file format
        self.format = self._detect_format(self.files)

        self.lat_name = lat_name
        self.lon_name = lon_name
        self.bathy_name = bathy_name
        self.lon_axis = lon_axis

    def _get_files(self, input):

        # check if input is a folder
        if os.path.isdir(input):
            files = list()
            extensions = ['.tif', 'xyz']
            for ext in extensions:
                files = get_files(path=input, substr=ext)
                if len(files) > 0:
                    break
        else:
            files = [input]

        return files

    def _detect_format(self, files):
        """ Detect file format.

            Args:
                files: list(str)
                    List of files

                fmt: Format
                    File format
        """
        fmt = None
        f = files[0]
        ext = f[f.rfind('.'):]

        if ext == '.nc':
            fmt = 'NETCDF'
        elif ext == '.mat':
            fmt = 'MATLAB'
        elif ext == '.tif':
            fmt = 'GEOTIFF'
        elif ext == '.xyz':
            fmt = 'XYZ'

        fmt = get_member(Format, fmt)

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
        d = Dataset(self.files[0])
        lat = np.array(d.variables[self.lat_name])
        lon = np.array(d.variables[self.lon_name])
        bathy = np.array(d.variables[self.bathy_name])

        # select region of interest
        lat, lon, bathy = self._select_region(lat, lon, bathy, latlon_SW, latlon_NE)

        return lat, lon, bathy

    def _read_matlab(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry from MATLAB file.

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
        m = sio.loadmat(self.files[0])
        lat = np.squeeze(np.array(m[self.lat_name]))
        lon = np.squeeze(np.array(m[self.lon_name]))
        bathy = np.array(m[self.bathy_name])

        # select region of interest
        lat, lon, bathy = self._select_region(lat, lon, bathy, latlon_SW, latlon_NE)

        return lat, lon, bathy

    def _read_geotiff(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry from a collection of GeoTIFF files.

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
                bathy: 1d numpy array
                    Bathymetry values
        """
        # empty lists
        lat, lon, bathy = list(), list(), list()

        # loop over files
        for f in self.files:
            x, y, z = self._read_single_geotiff(path=f, latlon_SW=latlon_SW, latlon_NE=latlon_NE)
            lon = lon + x
            lat = lat + y
            bathy = bathy + z

        return np.array(lat), np.array(lon), np.array(bathy)

    def _read_single_geotiff(self, path, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry from a single GeoTIFF file.

            Args: 
                path: str
                    Full path of the data file.
                latlon_NW: LatLon
                    South-western (SW) boundary of the region of interest.
                latlon_SE: LatLon
                    North-eastern (NE) boundary of the region of interest.

            Returns:
                x: list(float)
                    Longitude values
                y: list(float)
                    Latitude values
                z: list(float)
                    Bathymetry values
        """

        band_no = 1

        # parse south-west corner from file name
        sw = self._parse_sw_corner(path)

        # generate lat-lon arrays
        lats, lons = self._generate_latlon_arrays(sw)

        # get overlap with requested region
        lat_overlap = self._get_latitude_overlap(lats=lats, south=latlon_SW, north=latlon_NE)
        lon_overlap = self._get_longitude_overlap(lons=lons, west=latlon_SW, east=latlon_NE)
        if len(lat_overlap) == 0 or len(lon_overlap) == 0:
            return list(), list(), list()

        # open file
        data_set = gdal.Open(path)

        # get bathy matrix
        band = data_set.GetRasterBand(band_no)
        z = data_set.ReadAsArray()

        # get nodata value from the GDAL band object
        nodata = band.GetNoDataValue()

        # flip bathy matrix
        z = np.flip(z, axis=0)
        z = np.swapaxes(z, 0, 1)

        # grid
        x = lons
        y = lats
        x, y = np.meshgrid(x, y)

        # select entries with data
        has_data = np.where(z != nodata)
        x = x[has_data]
        y = y[has_data]
        z = z[has_data]

        # select data within the region of interest
        overlap = self._get_overlap(lats=y, lons=x, sw=latlon_SW, ne=latlon_NE)
        x = x[overlap]
        y = y[overlap]
        z = z[overlap]

        x = x.tolist()
        y = y.tolist()
        z = z.tolist()

        return x, y, z

    def _parse_sw_corner(self, path):
        """ Parse latitude and longitude data from filename assuming 
            the naming convention of the Canadian Hydrographic Service.

            Args: 
                path: str
                    File name

            Returns:
                sw_corner: LatLon
                    Latitude and longitude of the SW corner of the data set
        """
        f = path[path.rfind('/')+1:]
        north = int(f[4:8]) / 100
        west = int(f[9:14]) / 100
        east = west - 180

        assert east >= -180 and east <= 180, 'Invalid parsed longitude value'
        assert north >= -90 and north <= 180, 'Invalid parsed latitude value'

        sw_corner = LatLon(north, east)
        return sw_corner

    def _generate_latlon_arrays(self, sw_corner):
        """ Compute latitude and longitude values for a single 
            bathymetry file from the Canadian Hydrographic Service.

            The number of lat/lon values is 1001.

            The latitude binning is 0.001 degrees.

            The longitude binning depends on the latitude:

                * south of 68 deg N: 0.001 deg
                * latitudes from 68 to 80 deg N: 0.002 deg
                * 80 deg N and north: 0.004 deg

            Args: 
                sw_corner: LatLon
                    Latitude and longitude of the SW corner of the data set

            Returns:
                lats: numpy array
                    Uniformly spaced latitude values 
                lons: numpy array
                    Uniformly spaced longitude values 
        """
        # number of data points        
        N = 1001

        # latitude step size in degrees
        dlat = 0.001

        # longitude step size in degrees
        if sw_corner.longitude < 68:
            dlon = 0.001
        elif sw_corner.longitude >=68 and sw_corner.longitude < 80:
            dlon = 0.002
        elif sw_corner.longitude >= 80:
            dlon = 0.004

        lats = np.arange(N, dtype=np.float)
        lats *= dlat
        lats += sw_corner.latitude

        lons = np.arange(N, dtype=np.float)
        lons *= dlon
        lons += sw_corner.longitude

        return lats, lons

    def _get_latitude_overlap(self, lats, south, north):
        """ From an array of latitudes, selects those that 
            are within the prescribed boundaries.

            Args: 
                lats: numpy array
                    latitude values
                south: LatLon
                    Southern boundary
                north: LatLon
                    Northern boundary

            Returns:
                indeces: numpy array
                    Indeces of the entries that are within the prescribed boundaries
        """
        indeces = np.where(np.logical_and(lats <= north.latitude, lats >= south.latitude))
        return indeces

    def _get_longitude_overlap(self, lons, west, east):
        """ From an array of longitudes, selects those that 
            are within the prescribed boundaries.

            Args: 
                lons: numpy array
                    Longitude values
                west: LatLon
                    Western boundary
                east: LatLon
                    Eastern boundary

            Returns:
                indeces: numpy array
                    Indeces of the entries that are within the prescribed boundaries
        """
        indeces = np.where(np.logical_and(lons <= east.longitude, lons >= west.longitude))
        return indeces

    def _get_overlap(self, lats, lons, sw, ne):
        """ From arrays of latitude and longitude coordinates, selects those that 
            are within the prescribed boundaries.

            The latitude and longitude arrays must have the same length.

            Args: 
                lats: numpy array
                    latitude values
                lons: numpy array
                    Longitude values
                west: LatLon
                    Western boundary
                east: LatLon
                    Eastern boundary

            Returns:
                indeces: numpy array
                    Indeces of the entries that are within the prescribed boundaries
        """
        assert len(lats) == len(lons), 'lats and lons must have the same length'        

        lat_indeces = self._get_latitude_overlap(lats, sw, ne) 
        lon_indeces = self._get_longitude_overlap(lons, sw, ne) 

        indeces = np.intersect1d(lat_indeces, lon_indeces)

        return indeces

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
