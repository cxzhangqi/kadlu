""" Bathymetry interpolation module within the pyost package

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
from scipy.interpolate import RectBivariateSpline, RectSphereBivariateSpline
from bathy_reader import BathyReader, LatLon
from util import deg2rad, XYtoLL, LLtoXY, regXYgrid


class BathyInterpolator():
    """ Class for interpolating bathymetry data.

        Attributes: 
            bathy_reader: BathyReader
                Bathymetry data file reader
            latlon_SW: LatLon
                South-western (SW) boundary of the interpolation region.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the interpolation region.
            latlon_ref: LatLon
                Reference location (origo of XY coordinate system).
    """
    def __init__(self, bathy_reader, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180), origin=None):
        
        # read bathymetry data from file
        lat, lon, bathy = bathy_reader.read(latlon_SW, latlon_NE)

        # compute coordinates of origin, if not provided
        if origin is None:
            lat_ref = (lat[0] + lat[-1]) / 2
            lon_ref = (lon[0] + lon[-1]) / 2
            origin = LatLon(lat_ref, lon_ref)

        self.origin = origin

        # initialize lat-lon interpolator
        lat_rad, lon_rad = self._torad(lat, lon)
        self.interp_ll = RectSphereBivariateSpline(u=lat_rad, v=lon_rad, r=bathy)

        # define regular x-y grid
        x, y = regXYgrid(lat=lat, lon=lon, lat_ref=self.origin.latitude, lon_ref=self.origin.longitude)

        # transform to lat-lon
        lat_xy, lon_xy = XYtoLL(x=x, y=y, lat_ref=self.origin.latitude, lon_ref=self.origin.longitude, grid=True)

        # evaluate bathy on new grid using lat-lon interpolator
        Nx = len(x)
        Ny = len(y)
        bathy_xy = np.zeros(shape=(Nx,Ny))
        for i in range(Nx):
            bathy_xy[i,:] = self.eval_ll(lat=lat_xy[:,i], lon=lon_xy[:,i])

        # initialize x-y interpolator
        self.interp_xy = RectBivariateSpline(x=x, y=y, z=bathy_xy)

        # store grids
        self.lat_nodes = lat
        self.lon_nodes = lon
        self.x_nodes = x
        self.y_nodes = y

    def eval_xy(self, x, y, grid=False):
        """ Evaluate interpolated bathymetry in position coordinates (XY).

            x and y can be floats or arrays.

            If grid is set to False, the bathymetry will be evaluated at 
            the positions (x_i, y_i), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_N). Note that in this case, x and y must have 
            the same length.

            If grid is set to True, the bathymetry will be evaluated at 
            all combinations (x_i, y_j), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_M). Note that in this case, the lengths of x 
            and y do not have to be the same.

            Args: 
                x: float or array
                   x-coordinate of the positions(s) where the bathymetry is to be evaluated
                y: float or array
                   y-coordinate of the positions(s) where the bathymetry is to be evaluated
                grid: bool
                   Specify how to combine elements of x and y.

            Returns:
                zi: Interpolated bathymetry values
        """
        zi = self.interp_xy.__call__(x=x, y=y, grid=grid)
        return zi

    def eval_ll(self, lat, lon, grid=False):
        """ Interpolate bathymetry grid in latitude and longitude coordinates (LL).

            lat and lot can be floats or arrays.

            If grid is set to False, the bathymetry will be evaluated at 
            the coordinates (lat_i, lon_i), where lat=(lat_1,...,lat_N) 
            and lon=(lon_1,...,lon_N). Note that in this case, lat and 
            lon must have the same length.

            If grid is set to True, the bathymetry will be evaluated at 
            all combinations (lat_i, lon_j), where lat=(lat_1,...,lat_N) 
            and lon=(lon_1,...,lon_M). Note that in this case, the lengths 
            of lat and lon do not have to be the same.

            Args: 
                lat: float or array
                   latitude of the positions(s) where the bathymetry is to be evaluated
                lon: float or array
                   longitude of the positions(s) where the bathymetry is to be evaluated
                grid: bool
                   Specify how to combine elements of lat and lon.

            Returns:
                zi: Interpolated bathymetry values
        """
        lat = np.squeeze(np.array(lat))
        lon = np.squeeze(np.array(lon))
        lat_rad, lon_rad = self._torad(lat, lon)
        zi = self.interp_ll.__call__(theta=lat_rad, phi=lon_rad, grid=grid)
        return zi

    def _torad(self, lat, lon):
        lat_rad = (lat + 90) * deg2rad
        lon_rad = lon * deg2rad
        return lat_rad, lon_rad