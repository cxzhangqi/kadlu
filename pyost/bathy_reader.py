""" Bathymetry reader module within the pyost package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/pyost
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import numpy as np
from collections import namedtuple


LatLon = namedtuple('LatLon', ['latitude', 'longitude'])
""" Latitude and longitude coordinates for a given location.

    Args: 
        latitude: float
            Latitude in degrees from -90 (South Pole) to +90 (North Pole).
        longitude: float
            Longitude in degrees from -180 to +180 (West to East) with 0 
            corresponding to the Greenwich Meridian.
"""


class BathyReader():
    """ Parent class for bathymetry readers.

        Attributes: 
            path: str
                File name including path.
            lat_name: str
                Name of the variable that contains the latitue values.
            lon_name: str
                Name of the variable that contains the longitude values.
            bathy_name: str
                Name of the variable that contains the bathymetry values.
    """
    def __init__(self, path, lat_name='lat', lon_name='lon', bathy_name='bathy'):
        self.path = path
        self.lat_name = lat_name
        self.lon_name = lon_name
        self.bathy_name = bathy_name


from netCDF4 import Dataset

class BathyNetCDFReader(BathyReader):
    """ Class for reading bathymetry data from NetCDF files (*.nc).

        Attributes: 
            path: str
                File name including path.
            lat_name: str
                Name of the variable in the netCDF file that contains the latitue values.
            lon_name: str
                Name of the variable in the netCDF file that contains the longitude values.
            bathy_name: str
                Name of the variable in the netCDF file that contains the bathymetry values.
    """
    def __init__(self, path, lat_name='lat', lon_name='lon', bathy_name='elevation'):
        super().__init__(path, lat_name, lon_name, bathy_name) # initialize parent class

    def read(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry matrices from file.

            Args: 
                latlon_NW: LatLon
                    North-western (NW) boundary of the region of interest.
                latlon_SE: LatLon
                    South-eastern (SE) boundary of the region of interest.

            Returns:
                lat: 1d numpy array
                    Latitude values
                lon: 1d numpy array
                    Longitude values
                bathy: 2d numpy array
                    Bathymetry values
        """

        d = Dataset(self.path)

        lat = np.array(d.variables[self.lat_name])
        ind_lat = np.argwhere((lat >= latlon_SW.latitude) & (lat <= latlon_NE.latitude))
        ind_lat = np.squeeze(ind_lat)
        lat = lat[ind_lat]

        lon = np.array(d.variables[self.lon_name])
        ind_lon = np.argwhere((lon >= latlon_SW.longitude) & (lon <= latlon_NE.longitude))
        ind_lon = np.squeeze(ind_lon)
        lon = lon[ind_lon]

        ind = np.ix_(ind_lat, ind_lon)       

        bathy = d.variables[self.bathy_name]
        bathy = np.array(bathy)
        bathy = bathy[ind]

        return lat, lon, bathy


import scipy.io as sio

class BathyMatReader(BathyReader):
    """ Class for reading bathymetry data from MATLAB files (*.mat).

        Attributes: 
            path: str
                File name including path.
            lat_name: str
                Name of the variable in the MATLAB file that contains the latitue values.
            lon_name: str
                Name of the variable in the MATLAB file that contains the longitude values.
            bathy_name: str
                Name of the variable in the MATLAB file that contains the bathymetry values.
    """
    def __init__(self, path, lat_name='lat', lon_name='lon', bathy_name='elevation'):
        super().__init__(path, lat_name, lon_name, bathy_name) # initialize parent class

    def read(self, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and bathymetry matrices from file.

            Args: 
                latlon_NW: LatLon
                    North-western (NW) boundary of the region of interest.
                latlon_SE: LatLon
                    South-eastern (SE) boundary of the region of interest.

            Returns:
                lat: 1d or 2d numpy array
                    Latitude values
                lon: 1d or 2d numpy array
                    Longitude values
                bathy: 2d numpy array
                    Bathymetry values
        """

        m = sio.loadmat(self.path)

        lat = np.array(m[self.lat_name])
        lat = np.squeeze(lat)
        ind_lat = np.argwhere((lat >= latlon_SW.latitude) & (lat <= latlon_NE.latitude))
        ind_lat = np.squeeze(ind_lat)

        lon = np.array(m[self.lon_name])
        lon = np.squeeze(lon)
        ind_lon = np.argwhere((lon >= latlon_SW.longitude) & (lon <= latlon_NE.longitude))
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
        else:
            ind = np.ix_(ind_lat, ind_lon)       
    

        lat = lat[ind_lat]
        lon = lon[ind_lon]

        bathy = m[self.bathy_name][ind]
        bathy = np.array(bathy)

        print(bathy)

        return lat, lon, bathy
