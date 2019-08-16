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
from kadlu.utils import get_files


class BathyCHSReader():
    """ Class for loading Canadian Hydrographic Service Non-Navigational (NONNA-100) Bathymetric Data.
    
        Args:
            path: str
                Full path to the directory where the bathymetry data files are stored.
    """
    def __init__(self, path):

        # get files
        self.files = find_files(path, ext='tiff')
        assert len(self.files) > 0, "Could not find any data files" 

    def _geographic_overlap():

    def load(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Load available bathymetry data within specified geographical region.

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
                bathy: 1d numpy array
                    Bathymetry values
        """
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

        # lat-lon range
        lat_min = np.min(lats)
        lat_max = np.max(lats)
        lon_min = np.min(lons)
        lon_max = np.max(lons)

        # check if overlap
        lat_overlap = (lat_min <= latlon_NE.latitude and lat_max >= latlon_SW.latitude)
        lon_overlap = (lon_min <= latlon_NE.longitude and lon_max >= latlon_SW.longitude)

        if not (lat_overlap and lon_overlap):
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

        # check that shape fits with size of lat and lon arrays
        # if not, issue a warning and re-compute the lat and lon arrays
        # to ensure that they fit
        num_lat = z.shape[0]
        num_lon = z.shape[1]
        if num_lat != len(lats) or num_lon != len(lons):
            print('Warning: Bathymetry data from {0} has shape {1} x {2} whereas {3} x {4} was expected'.format(path, num_lat, num_lon, len(lats), len(lons)))
            lats, lons = self._generate_latlon_arrays(sw, num_lat=num_lat, num_lon=num_lon)

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
        east = -west

        assert east >= -180 and east <= 180, 'Invalid parsed longitude value'
        assert north >= -90 and north <= 180, 'Invalid parsed latitude value'

        sw_corner = LatLon(north, east)
        return sw_corner

    def _generate_latlon_arrays(self, sw_corner, num_lat=1001, num_lon=1001):
        """ Compute latitude and longitude values for a single 
            bathymetry file from the Canadian Hydrographic Service.

            The number of lat/lon values is 1001 by default.

            The latitude binning is 0.001 degrees.

            The longitude binning depends on the latitude:

                * south of 68 deg N: 0.001 deg
                * latitudes from 68 to 80 deg N: 0.002 deg
                * 80 deg N and north: 0.004 deg

            Args: 
                sw_corner: LatLon
                    Latitude and longitude of the SW corner of the data set
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
        indeces = np.squeeze(indeces)
        if len(indeces.shape) == 0:
            indeces = indeces[np.newaxis]

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
        indeces = np.squeeze(indeces)
        if len(indeces.shape) == 0:
            indeces = indeces[np.newaxis]

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
