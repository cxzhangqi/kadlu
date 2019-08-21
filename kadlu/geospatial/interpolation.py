# ================================================================================ #
#   Authors: Casey Hillard and Oliver Kirsebom                                     #
#   Contact: oliver.kirsebom@dal.ca                                                #
#   Organization: MERIDIAN (https://meridian.cs.dal.ca/)                           #
#   Team: Data Analytics                                                           #
#   Project: kadlu                                                                 #
#   Project goal: The kadlu library provides functionalities for modeling          #
#   underwater noise due to environmental source such as waves.                    #
#                                                                                  #
#   License: GNU GPLv3                                                             #
#                                                                                  #
#       This program is free software: you can redistribute it and/or modify       #
#       it under the terms of the GNU General Public License as published by       #
#       the Free Software Foundation, either version 3 of the License, or          #
#       (at your option) any later version.                                        #
#                                                                                  #
#       This program is distributed in the hope that it will be useful,            #
#       but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#       GNU General Public License for more details.                               # 
#                                                                                  #
#       You should have received a copy of the GNU General Public License          #
#       along with this program.  If not, see <https://www.gnu.org/licenses/>.     #
# ================================================================================ #

""" Interpolation module within the kadlu library

    This module handles two- and three-dimensional interpolation of geospatial data in 
    spherical and planar geometry.

    In the two-dimensional case, the interpolation can be made on both regular and 
    irregular grids. 
    
    In the three-dimensional case, only interpolation on regular grids has been 
    implemented, although an extension to irregular grids (following the same 
    methodology as in the two-dimensional case) should be straightforward. 

    Contents:
        GridData2D class:
        Interpolator2D class
        Interpolator3D class
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline, RectSphereBivariateSpline
from scipy.interpolate import RegularGridInterpolator
from kadlu.geospatial.bathy_reader import LatLon
from kadlu.utils import deg2rad, XYtoLL, LLtoXY, torad, DLDL_over_DXDY
from scipy.interpolate import griddata

from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from matplotlib import pyplot as plt


class GridData2D():
    """ Wrapper function around scipy's interpolate.griddata

        https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.interpolate.griddata.html

        An alternative to griddata could be Rbf, as discussed here:

        https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy

        Attributes: 
            u: 1d numpy array
                data points 1st coordinate
            v: 1d numpy array
                data points 2nd coordinate
            r: 1d numpy array
                data values
            method : {‘linear’, ‘nearest’, ‘cubic’}, optional
    """
    def __init__(self, u, v, r, method='cubic'):
        self.uv = np.column_stack((u,v))
        self.r = r
        self.method = method
        self.u_step = (np.max(u) - np.min(u)) / 1E4
        self.v_step = (np.max(v) - np.min(v)) / 1E4

    def __call__(self, theta, phi, grid=False, dtheta=0, dphi=0):
        """ Interpolate data

            theta and phi can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the positions (theta_i, phi_i), where theta=(theta_1,...,theta_N) and 
            phi=(phi_1,...,phi_N). Note that in this case, theta and phi must have 
            the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (theta_i, phi_j), where theta=(theta_1,...,theta_N) and 
            phi=(phi_1,...,phi_M). Note that in this case, the lengths of theta 
            and phi do not have to be the same.

            Args: 
                theta: float or array
                   1st coordinate of the points where the interpolation is to be evaluated
                phi: float or array
                   2nd coordinate of the points where the interpolation is to be evaluated
                grid: bool
                   Specify how to combine elements of theta and phi.
                dtheta: int
                    Order of theta-derivative
                dphi: int
                    Order of phi-derivative

            Returns:
                ri: Interpolated values
        """        
        assert dtheta + dphi > 0, "Interpolation of higher-order derivatives not implemented for irregular grids"

        if dtheta == dphi == 0:
            pts = self._prep_input(theta, phi, grid)
            ri = griddata(self.uv, self.r, pts, method=self.method)

        # quick and dirty interpolation of 1st derivative in u
        elif dtheta == 1 and dphi == 0:
            theta1 = theta - 0.5 * self.u_step
            theta2 = theta + 0.5 * self.u_step
            pts1 = self._prep_input(theta1, phi, grid)
            pts2 = self._prep_input(theta2, phi, grid)
            ri1 = griddata(self.uv, self.r, pts1, method=self.method)
            ri2 = griddata(self.uv, self.r, pts2, method=self.method)
            ri = (ri2 - ri1) / (theta2 - theta1)

        # quick and dirty interpolation of 1st derivative in v
        elif dtheta == 0 and dphi == 1:
            phi1 = phi - 0.5 * self.v_step
            phi2 = phi + 0.5 * self.v_step
            pts1 = self._prep_input(theta, phi1, grid)
            pts2 = self._prep_input(theta, phi2, grid)
            ri1 = griddata(self.uv, self.r, pts1, method=self.method)
            ri2 = griddata(self.uv, self.r, pts2, method=self.method)
            ri = (ri2 - ri1) / (phi2 - phi1)

        if grid:
            ri = np.reshape(ri, newshape=(M,N))

        return ri

    def _prep_input(self, theta, phi, grid):
        """ Transform the input data to format appropriate for scipy.interpolate.griddata.

            Args: 
                theta: 1d numpy array
                   1st coordinate of the points where the interpolation is to be evaluated
                phi: 1d numpy array
                   2nd coordinate of the points where the interpolation is to be evaluated
                grid: bool
                   Specify how to combine elements of theta and phi.

            Returns:
                pts: numpy array 
                    Input data
        """        
        if grid:
            M = len(theta)
            N = len(phi)
            theta, phi = np.meshgrid(theta, phi)
            theta = np.reshape(theta, newshape=(M*N))
            phi = np.reshape(phi, newshape=(M*N))        

        pts = np.column_stack((theta, phi))

        return pts


class Interpolator2D():
    """ Class for interpolating 2D (lat,lon) geospatial data.

        For irregular grids, the data values must be passed as a 
        1d array and all three arrays (values, lats, lons) must have 
        the same length.

        For regular grids, the data values must be passed as a 
        2d array with shape (M,N) where M and N are the lengths 
        of the latitude and longitude array, respectively.

        Attributes: 
            values: 1d or 2d numpy array
                Values to be interpolated
            lats: 1d numpy array
                Latitude values
            lons: 1d numpy array
                Longitude values
            latlon_ref: LatLon
                Reference location (origo of XY coordinate system).
            method : {‘linear’, ‘nearest’, ‘cubic’}, optional
                Interpolation method used only for irregular grids.
    """
    def __init__(self, values, lats, lons, origin=None, method='cubic'):
        
        # compute coordinates of origin, if not provided
        if origin is None:
            lat_ref = (np.min(lats) + np.max(lats)) / 2
            lon_ref = (np.min(lons) + np.max(lons)) / 2
            origin = LatLon(lat_ref, lon_ref)

        self.origin = origin

        # check if bathymetry data are on a regular or irregular grid
        reggrid = (np.ndim(values) == 2)

        # convert to radians
        lats_rad, lons_rad = torad(lats, lons)

        # initialize lat-lon interpolator
        if reggrid:
            self.interp_ll = RectSphereBivariateSpline(u=lats_rad, v=lons_rad, r=values)
        else:
            self.interp_ll = GridData2D(u=lats_rad, v=lons_rad, r=values, method=method)

        # store grids
        self.lat_nodes = lats
        self.lon_nodes = lons
        self.values = values

    def eval_xy(self, x, y, grid=False, x_deriv_order=0, y_deriv_order=0):
        """ Interpolate using planar coordinate system (xy).

            x and y can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the positions (x_i, y_i), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_N). Note that in this case, x and y must have 
            the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (x_i, y_j), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_M). Note that in this case, the lengths of x 
            and y do not have to be the same.

            Args: 
                x: float or array
                   x-coordinate of the positions(s) where the interpolation is to be evaluated
                y: float or array
                   y-coordinate of the positions(s) where the interpolation is to be evaluated
                grid: bool
                   Specify how to combine elements of x and y.
                x_deriv_order: int
                    Order of x-derivative
                y_deriv_order: int
                    Order of y-derivative

            Returns:
                zi: Interpolated interpolation values
        """
        lat, lon = self._xy_to_ll(x, y, grid=grid)

        if grid:
            M = lat.shape[0]
            N = lat.shape[1]
            lat = np.reshape(lat, newshape=(M*N))
            lon = np.reshape(lon, newshape=(M*N))

        zi = self.eval_ll(lat=lat, lon=lon, lat_deriv_order=y_deriv_order, lon_deriv_order=x_deriv_order)

        if x_deriv_order + y_deriv_order > 0:
            r = DLDL_over_DXDY(lat=lat, lat_deriv_order=y_deriv_order, lon_deriv_order=x_deriv_order)
            zi *= r

        if grid:
            zi = np.reshape(zi, newshape=(M,N))

        if np.ndim(zi) == 2:
            zi = np.swapaxes(zi, 0, 1)

        if np.ndim(zi) == 0 or (np.ndim(zi) == 1 and len(zi) == 1):
            zi = float(zi)

        return zi

    def eval_ll(self, lat, lon, grid=False, lat_deriv_order=0, lon_deriv_order=0):
        """ Interpolate using spherical coordinate system (latitude-longitude).

            lat and lot can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the coordinates (lat_i, lon_i), where lat=(lat_1,...,lat_N) 
            and lon=(lon_1,...,lon_N). Note that in this case, lat and 
            lon must have the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (lat_i, lon_j), where lat=(lat_1,...,lat_N) 
            and lon=(lon_1,...,lon_M). Note that in this case, the lengths 
            of lat and lon do not have to be the same.

            Derivates are given per radians^n, where n is the overall 
            derivative order.

            Args: 
                lat: float or array
                    latitude of the positions(s) where the interpolation is to be evaluated
                lon: float or array
                    longitude of the positions(s) where the interpolation is to be evaluated
                grid: bool
                    Specify how to combine elements of lat and lon.
                lat_deriv_order: int
                    Order of latitude-derivative
                lon_deriv_order: int
                    Order of longitude-derivative

            Returns:
                zi: Interpolated values (or derivates)
        """
        lat = np.squeeze(np.array(lat))
        lon = np.squeeze(np.array(lon))
        lat_rad, lon_rad = torad(lat, lon)

        zi = self.interp_ll.__call__(theta=lat_rad, phi=lon_rad, grid=grid, dtheta=lat_deriv_order, dphi=lon_deriv_order)

        if np.ndim(self.values) == 1 and np.ndim(zi) == 2:
            zi = np.swapaxes(zi, 0, 1)

        if np.ndim(zi) == 0 or (np.ndim(zi) == 1 and len(zi) == 1):
            zi = float(zi)

        return zi

    def slice(self, angle=0, distance=1000, bins=100, num_slices=1, origin=None):

        # x,y coordinates of origin
        if origin is None:
            origin = self.origin
            xo, yo = 0, 0
        else:
            xo, yo = LLtoXY(lat=origin.latitude, lon=origin.longitude,\
                lat_ref=self.origin.latitude, lon_ref=self.origin.longitude)

        # distance array
        dr = distance / float(bins)
        r = np.arange(bins, dtype=np.float)
        r *= dr
        r += 0.5 * dr

        val = list()

        # loop over angles
        a = angle
        da = 360. / float(num_slices)
        for _ in range(num_slices):
            x = r * np.cos(a * np.pi / 180.)
            y = r * np.sin(a * np.pi / 180.)
            x += xo
            y += yo
            b = self.eval_xy(x=x, y=y)
            val.append(b)
            a += da

        if num_slices == 1:
            val = val[0]
        
        return val

    def _ll_to_xy(self, lat, lon, lat_ref=None, lon_ref=None, grid=False):
        """ Convert from spherical (lat-lon) to planar (x-y) coordinate system.

            Args: 
                lat: float or array
                    latitude(s)
                lon: float or array
                    longitude(s)
                lat_ref: float
                    reference latitude
                lon_ref: float
                    reference longitude
                grid: bool
                    Specify how to combine elements of lat and lon.

            Returns:
                x: float or array
                   x value(s)
                y: float or array
                   y value(s)
        """
        if lat_ref is None:
            lat_ref = self.origin.latitude
        if lon_ref is None:
            lon_ref = self.origin.longitude

        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=lat_ref, lon_ref=lon_ref, grid=grid)
        return x, y

    def _xy_to_ll(self, x, y, lat_ref=None, lon_ref=None, grid=False):
        """ Convert from planar (x-y) to spherical (lat-lon) coordinate system.

            Args: 
                x: float or array
                   x value(s)
                y: float or array
                   y value(s)
                lat_ref: float
                    reference latitude
                lon_ref: float
                    reference longitude
                grid: bool
                    Specify how to combine elements of lat and lon.

            Returns:
                lat: float or array
                    latitude(s)
                lon: float or array
                    longitude(s)
        """
        if lat_ref is None:
            lat_ref = self.origin.latitude
        if lon_ref is None:
            lon_ref = self.origin.longitude

        lat, lon = XYtoLL(x=x, y=y, lat_ref=lat_ref, lon_ref=lon_ref, grid=grid)
        return lat, lon


class Interpolator3D():
    """ Class for interpolating 3D (lat,lon,depth) geospatial data.

        For regular grids, the data values must be passed as a 
        3d array with shape (M,N,K) where M,N,K are the lengths 
        of the latitude, longitude, and depth arrays, respectively.

        The current implementation does not handle irregular grids, 
        although an extension to irregular grids (following the same 
        methodology as in the two-dimensional case) should be 
        straightforward. 

        Attributes: 
            values: 3d numpy array
                Values to be interpolated
            lats: 1d numpy array
                Latitude values
            lons: 1d numpy array
                Longitude values
            depths: 1d numpy array
                Depth values
            latlon_ref: LatLon
                Reference location (origo of XY coordinate system).
            method : {‘linear’, ‘nearest’}, optional
                Interpolation method. Default is linear
    """
    def __init__(self, values, lats, lons, depths, origin=None, method='linear'):
        
        # compute coordinates of origin, if not provided
        if origin is None:
            lat_ref = (np.min(lats) + np.max(lats)) / 2
            lon_ref = (np.min(lons) + np.max(lons)) / 2
            origin = LatLon(lat_ref, lon_ref)

        self.origin = origin

        # check if bathymetry data are on a regular or irregular grid
        assert np.ndim(values) == 3, 'values must be 3-dimensional'

        # convert to radians
        lats_rad, lons_rad = torad(lats, lons)

        # initialize lat-lon interpolator
        self.interp_ll = RegularGridInterpolator((lats_rad, lons_rad, depths), values)

        # store grids
        self.lat_nodes = lats
        self.lon_nodes = lons
        self.depth_nodes = depths
        self.values = values

    def eval_xyz(self, x, y, z, v, grid=False):
        """ Interpolate using planar coordinate system (xy).

            x,y,z can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the positions (x_i, y_i, z_i), where x=(x_1,...,x_N),  
            y=(y_1,...,y_N), and z=(z_1,...,z_N). Note that in this case, 
            x,y,z must have the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (x_i, y_j, z_k), where x=(x_1,...,x_N), 
            y=(y_1,...,y_M), and z=(z_1,...,z_K). Note that in this case, the 
            lengths of x,y,z do not have to be the same.

            Args: 
                x: float or array
                   x-coordinate of the positions(s) where the interpolation is to be evaluated
                y: float or array
                   y-coordinate of the positions(s) where the interpolation is to be evaluated
                z: float or array
                   z-coordinate of the positions(s) where the interpolation is to be evaluated
                grid: bool
                   Specify how to combine elements of x,y,z.

            Returns:
                vi: Interpolated values
        """
        lat, lon = self._xy_to_ll(x, y, grid=grid)

        if grid:
            M = lat.shape[0]
            N = lat.shape[1]
            K = z.shape[0]
            lat = np.reshape(lat, newshape=(M*N*K))
            lon = np.reshape(lon, newshape=(M*N*K))
            z = np.reshape(z, newshape=(M*N*K))

        zi = self.eval_ll(lat=lat, lon=lon, depths=z)

        if x_deriv_order + y_deriv_order > 0:
            r = DLDL_over_DXDY(lat=lat, lat_deriv_order=y_deriv_order, lon_deriv_order=x_deriv_order)
            zi *= r

        if grid:
            zi = np.reshape(zi, newshape=(M,N))

        if np.ndim(zi) == 2:
            zi = np.swapaxes(zi, 0, 1)

        if np.ndim(zi) == 0 or (np.ndim(zi) == 1 and len(zi) == 1):
            zi = float(zi)

        return zi

    def eval_ll(self, lat, lon, grid=False, lat_deriv_order=0, lon_deriv_order=0):
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

            Bathymetry values are given in meters and derivates are given in meters 
            per radians^n, where n is the overall derivative order.

            Args: 
                lat: float or array
                    latitude of the positions(s) where the bathymetry is to be evaluated
                lon: float or array
                    longitude of the positions(s) where the bathymetry is to be evaluated
                grid: bool
                    Specify how to combine elements of lat and lon.
                lat_deriv_order: int
                    Order of latitude-derivative
                lon_deriv_order: int
                    Order of longitude-derivative

            Returns:
                zi: Interpolated bathymetry values (or derivates)
        """
        lat = np.squeeze(np.array(lat))
        lon = np.squeeze(np.array(lon))
        lat_rad, lon_rad = torad(lat, lon)

        zi = self.interp_ll.__call__(theta=lat_rad, phi=lon_rad, grid=grid, dtheta=lat_deriv_order, dphi=lon_deriv_order)

        if np.ndim(self.values) == 1 and np.ndim(zi) == 2:
            zi = np.swapaxes(zi, 0, 1)

        if np.ndim(zi) == 0 or (np.ndim(zi) == 1 and len(zi) == 1):
            zi = float(zi)

        return zi