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


def crop(lat, lon, latlon_SW, latlon_NE, grid=False):
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
    ind_lat = np.argwhere((lat >= latlon_SW.latitude) & (lat <= latlon_NE.latitude))
    ind_lat = np.squeeze(ind_lat)

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

    if grid:
        ind = np.ix_(ind_lat, ind_lon)       
    else:
        ind = np.intersect1d(ind_lat, ind_lon)
        ind_lat = ind
        ind_lon = ind

    lat = lat[ind_lat]
    lon = lon[ind_lon]

    return ind, lat, lon