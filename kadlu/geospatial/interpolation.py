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
        Interpolator2D class:
        Interpolator3D class:
        Uniform2D class:
        Uniform3D class:
        DepthInterpolator3D class
"""

import numpy as np
from scipy.interpolate import RectSphereBivariateSpline, RegularGridInterpolator, interp1d, interp2d, griddata
from kadlu.utils import deg2rad, XYtoLL, LLtoXY, torad, DLDL_over_DXDY, LatLon

from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from matplotlib import pyplot as plt


def interp_2D(values, lats=None, lons=None, origin=None, 
            method_irreg='cubic',lats_reg=None, lons_reg=None):

        if isinstance(values, (float, int)):
            return Uniform2D(values)
        
        else:
            return Interpolator2D(values, lats, lons, origin, method_irreg, lats_reg, lons_reg)

def interp_3D(values, lats=None, lons=None, depths=None, origin=None, method='linear'):

        if isinstance(values, (float, int)):
            return Uniform3D(values)
        
        else:
            return Interpolator3D(values, lats, lons, depths, origin, method)


class GridData2D():
    """ Interpolation of data on a two-dimensional irregular grid.
    
        Essentially, a wrapper function around scipy's interpolate.griddata.

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

            Only first-order derivatives have been implemented. A request to interpolate 
            higher-order derivatives will give an Assertion error.

            TODO: Improve the algorithm used to compute the derivatives. 

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
        assert dtheta + dphi <= 1, "Interpolation of higher-order derivatives not implemented for irregular grids"

        if grid:
            M = N = 1
            if np.ndim(theta) == 1: 
                M = len(theta)
            if np.ndim(phi) == 1: 
                N = len(phi)

        if dtheta == dphi == 0:
            pts = self._prep_input(theta, phi, grid)
            ri = griddata(self.uv, self.r, pts, method=self.method)

        # 1st derivative in u
        elif dtheta == 1 and dphi == 0:
            theta1 = theta - 0.5 * self.u_step
            theta2 = theta + 0.5 * self.u_step
            pts1 = self._prep_input(theta1, phi, grid)
            pts2 = self._prep_input(theta2, phi, grid)
            ri1 = griddata(self.uv, self.r, pts1, method=self.method)
            ri2 = griddata(self.uv, self.r, pts2, method=self.method)
            ri = (ri2 - ri1) / (theta2 - theta1)

        # 1st derivative in v
        elif dtheta == 0 and dphi == 1:
            phi1 = phi - 0.5 * self.v_step
            phi2 = phi + 0.5 * self.v_step
            pts1 = self._prep_input(theta, phi1, grid)
            pts2 = self._prep_input(theta, phi2, grid)
            ri1 = griddata(self.uv, self.r, pts1, method=self.method)
            ri2 = griddata(self.uv, self.r, pts2, method=self.method)
            ri = (ri2 - ri1) / (phi2 - phi1)

        if grid:
            ri = np.reshape(ri, newshape=(N,M))
            ri = np.swapaxes(ri, 0, 1)

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
            method_irreg : {‘linear’, ‘nearest’, ‘cubic’, ‘regular’}, optional
                Interpolation method used for irregular grids.
                Note that 'nearest' is usually significantly faster than 
                the 'linear' and 'cubic'.
                If the 'regular' is selected, the data is first mapped onto 
                a regular grid by means of a cubic interpolation (for points outside 
                the area covered by the data, a nearest-point interpolation is used).
                The coordiantes of the regular grid onto which the data is mapped 
                is given by the arguments 'lats_reg' and 'lons_reg'.
            lats_reg: 1d numpy array
                Latitude values for regular interpolation grid.
                Must be specified for irregular grids if the interpolation method 
                'regular' is chosen. In all other cases, the argument is ignored.
            lons_reg: 1d numpy array
                Longitude values for regular interpolation grid.
                Must be specified for irregular grids if the interpolation method 
                'regular' is chosen. In all other cases, the argument is ignored.
    """
    def __init__(self, values, lats, lons,
            origin=None, method_irreg='cubic',
            lats_reg=None, lons_reg=None):
        
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

        # necessary to resolve a mismatch between scipy and underlying Fortran code
        # https://github.com/scipy/scipy/issues/6556
        if np.min(lons_rad) < 0:
            self._lon_corr = np.pi
        else:
            self._lon_corr = 0

        lons_rad += self._lon_corr

        # initialize lat-lon interpolator
        if reggrid:
            if len(lats) > 2 and len(lons) > 2:
                self.interp_ll = RectSphereBivariateSpline(u=lats_rad, v=lons_rad, r=values)
            elif len(lats) > 1 and len(lons) > 1:
                z = np.swapaxes(values, 0, 1)
                self.interp_ll = interp2d(x=lats_rad, y=lons_rad, z=z, kind='linear')
            elif len(lats) == 1:
                self.interp_ll = interp1d(x=lons_rad, y=np.squeeze(values), kind='linear')
            elif len(lons) == 1:
                self.interp_ll = interp1d(x=lats_rad, y=np.squeeze(values), kind='linear')

        else:
            if method_irreg == 'regular':
                assert lats_reg is not None and lons_reg is not None,\
                    'lats_reg and lons_reg must be specified for irregular grids when the interpolation method is `regular`'
    
                # interpolators on irregular grid
                gd_cubic = GridData2D(u=lats_rad, v=lons_rad, r=values, method='linear') #method='cubic')
                gd_nearest = GridData2D(u=lats_rad, v=lons_rad, r=values, method='nearest')

                # map to regular grid
                lats_reg_rad, lons_reg_rad = torad(lats_reg, lons_reg)
                lons_reg_rad += self._lon_corr
                zi = gd_cubic.__call__(theta=lats_reg_rad, phi=lons_reg_rad, grid=True)
                zi_nearest = gd_nearest.__call__(theta=lats_reg_rad, phi=lons_reg_rad, grid=True)
                indices_nan = np.where(np.isnan(zi))
                zi[indices_nan] = zi_nearest[indices_nan] 

                # interpolator on regular grid
                self.interp_ll = RectSphereBivariateSpline(u=lats_reg_rad, v=lons_reg_rad, r=zi)

            else:
                self.interp_ll = GridData2D(u=lats_rad, v=lons_rad, r=values, method=method_irreg)

        # store data used for interpolation
        self.lat_nodes = lats
        self.lon_nodes = lons
        self.values = values

    def get_nodes(self):
        return (self.values, self.lat_nodes, self.lon_nodes)

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
        lat, lon = XYtoLL(x=x, y=y, lat_ref=self.origin.latitude, lon_ref=self.origin.longitude, grid=grid)

        if grid:
            M = lat.shape[0]
            N = lat.shape[1]
            lat = np.reshape(lat, newshape=(M*N))
            lon = np.reshape(lon, newshape=(M*N))

        zi = self.eval_ll(lat=lat, lon=lon, squeeze=False, lat_deriv_order=y_deriv_order, lon_deriv_order=x_deriv_order)

        if x_deriv_order + y_deriv_order > 0:
            r = DLDL_over_DXDY(lat=lat, lat_deriv_order=y_deriv_order, lon_deriv_order=x_deriv_order)
            zi *= r

        if grid:
            zi = np.reshape(zi, newshape=(M,N))

        if np.ndim(zi) == 2:
            zi = np.swapaxes(zi, 0, 1)

        zi = np.squeeze(zi)

        if np.ndim(zi) == 0 or (np.ndim(zi) == 1 and len(zi) == 1):
            zi = float(zi)

        return zi

    def eval_ll(self, lat, lon, grid=False, squeeze=True, lat_deriv_order=0, lon_deriv_order=0):
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
                    Specify how to combine elements of lat and lon. If lat and lon have different
                    lengths, specifying grid has no effect as it is automatically set to True.
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
        lon_rad += self._lon_corr

        if isinstance(self.interp_ll, interp2d):
            zi = self.interp_ll.__call__(x=lat_rad, y=lon_rad, dx=lat_deriv_order, dy=lon_deriv_order)
            if grid: zi = np.swapaxes(zi, 0, 1)
            if not grid and np.ndim(zi) == 2: zi = np.diagonal(zi)

        elif isinstance(self.interp_ll, interp1d):
            if len(self.lat_nodes) > 1:
                zi = self.interp_ll(x=lat_rad)
            elif len(self.lon_nodes) > 1:
                zi = self.interp_ll(x=lon_rad)

        else:
            zi = self.interp_ll.__call__(theta=lat_rad, phi=lon_rad, grid=grid, dtheta=lat_deriv_order, dphi=lon_deriv_order)

        if squeeze:
            zi = np.squeeze(zi)

        if np.ndim(zi) == 0 or (np.ndim(zi) == 1 and len(zi) == 1):
            zi = float(zi)

        return zi


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
        assert np.ndim(values) == 3, 'values must be a 3d array'

        # convert to radians
        lats_rad, lons_rad = torad(lats, lons)

        # necessary to resolve a mismatch between scipy and underlying Fortran code
        # https://github.com/scipy/scipy/issues/6556
        if np.min(lons_rad) < 0:
            self._lon_corr = np.pi
        else:
            self._lon_corr = 0

        # initialize lat-lon interpolator
        lons_rad += self._lon_corr
        self.interp_ll = RegularGridInterpolator((lats_rad, lons_rad, depths), values, method=method, bounds_error=False, fill_value=None)

        # store grids
        self.lat_nodes = lats
        self.lon_nodes = lons
        self.depth_nodes = depths
        self.values = values

    def get_nodes(self):
        return (self.values, self.lat_nodes, self.lon_nodes, self.depth_nodes)

    def eval_xy(self, x, y, z, grid=False):
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
                   Depth of the positions(s) where the interpolation is to be evaluated
                grid: bool
                   Specify how to combine elements of x,y,z.

            Returns:
                vi: Interpolated values
        """
        M = N = K = 1
        if np.ndim(y) == 1: 
            M = len(y)
        if np.ndim(x) == 1: 
            N = len(x)
        if np.ndim(z) == 1: 
            K = len(z)

        lat, lon, z = XYtoLL(x=x, y=y, lat_ref=self.origin.latitude, lon_ref=self.origin.longitude, grid=grid, z=z)

        if grid:
            lat = lat.flatten()
            lon = lon.flatten()
            z = z.flatten()

        vi = self.eval_ll(lat=lat, lon=lon, z=z, squeeze=False)

        if grid:
            vi = np.reshape(vi, newshape=(M,N,K))

        if np.ndim(vi) == 3:
            vi = np.swapaxes(vi, 0, 1)

        vi = np.squeeze(vi)

        if np.ndim(vi) == 0 or (np.ndim(vi) == 1 and len(vi) == 1):
            vi = float(vi)

        return vi

    def eval_ll(self, lat, lon, z, grid=False, squeeze=True):
        """ Interpolate using spherical coordinate system (lat-lon).

            lat,lot,z can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the coordinates (lat_i, lon_i, z_i), where lat=(lat_1,...,lat_N), 
            lon=(lon_1,...,lon_N) and z=(z_,...,z_K). Note that in this case, lat and 
            lon must have the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (lat_i, lon_j, z_k), where lat=(lat_1,...,lat_N), 
            lon=(lon_1,...,lon_M) and z=(z_1,...,z_K). Note that in this case, the lengths 
            of lat and lon do not have to be the same.

            Args: 
                lat: float or array
                    Latitude of the positions(s) where the bathymetry is to be evaluated
                lon: float or array
                    Longitude of the positions(s) where the bathymetry is to be evaluated
                z: float or array
                    Depth of the positions(s) where the interpolation is to be evaluated
                grid: bool
                    Specify how to combine elements of lat,lon,z.

            Returns:
                zi: Interpolated bathymetry values (or derivates)
        """
        M = N = K = 1
        if np.ndim(lat) == 1: 
            M = len(lat)
        if np.ndim(lon) == 1: 
            N = len(lon)
        if np.ndim(z) == 1: 
            K = len(z)

        lat = np.squeeze(np.array(lat))
        lon = np.squeeze(np.array(lon))
        lat_rad, lon_rad = torad(lat, lon)
        lon_rad += self._lon_corr

        z = np.squeeze(np.array(z))

        if grid:
            lat_rad, lon_rad, z = np.meshgrid(lat_rad, lon_rad, z)
            lat_rad = lat_rad.flatten()
            lon_rad = lon_rad.flatten()
            z = z.flatten()

        pts = np.column_stack((lat_rad, lon_rad, z))        
        vi = self.interp_ll(pts)

        if grid:
            vi = np.reshape(vi, newshape=(M,N,K))

        if squeeze:
            vi = np.squeeze(vi)

        if np.ndim(vi) == 0 or (np.ndim(vi) == 1 and len(vi) == 1):
            vi = float(vi)

        return vi


class Uniform2D():

    """
    def __init__(self, values, lats, lons):
        self.value = values[0]
    """
    def __init__(self, values, lats, lons):
        self.value = values

    def get_nodes(self):
        return self.value

    def eval_xy(self, x, y, grid=False, x_deriv_order=0, y_deriv_order=0):

        z = self.eval_ll(lat=y, lon=x, grid=grid, squeeze=False, lat_deriv_order=y_deriv_order, lon_deriv_order=x_deriv_order)

        if np.ndim(z) == 3:
            z = np.swapaxes(z, 0, 1)

        z = np.squeeze(z)

        return z

    def eval_ll(self, lat, lon, grid=False, squeeze=True, lat_deriv_order=0, lon_deriv_order=0):

        if np.ndim(lat) == 0: lat = np.array([lat])
        if np.ndim(lon) == 0: lon = np.array([lon])

        if grid:
            s = (len(lat), len(lon))

        else:
            assert len(lat) == len(lon), 'when grid is False, lat and lon must have the same length'

            s = len(lat)

        if lat_deriv_order + lon_deriv_order > 0:
            v = 0
        else:
            v = self.value

        z = np.ones(s) * v

        if squeeze:
            z = np.squeeze(z)

        return z


class Uniform3D():

    #def __init__(self, value):
    #    self.value = value
    """
    def __init__(self, values, lats, lons, depths):
        print(values[0][0][0])
        self.value = values.flatten()[0]
    """
    def __init__(self, values):
        self.value = values

    def get_nodes(self):
        return self.value

    def eval_xy(self, x, y, z, grid=False):

        v = self.eval_ll(lat=y, lon=x, z=z, grid=grid, squeeze=False)

        if np.ndim(v) == 3:
            v = np.swapaxes(v, 0, 1)

        v = np.squeeze(v)

        return v

    def eval_ll(self, lat, lon, z, grid=False, squeeze=True):

        if np.ndim(lat) == 0: lat = np.array([lat])
        if np.ndim(lon) == 0: lon = np.array([lon])
        if np.ndim(z) == 0: z = np.array([z])

        if grid:
            s = (len(lat), len(lon), len(z))

        else:
            assert len(lat) == len(lon) == len(z), 'when grid is False, lat,lon,z must have the same length'

            s = len(lat)

        v = np.ones(s) * self.value
        
        if squeeze:
            v = np.squeeze(v)

        return v


class DepthInterpolator3D():
    """ Class for interpolating 3D (lat,lon,depth) geospatial data
        that only varies with depth.

        Attributes: 
            values: 1d or 3d numpy array
                Values to be interpolated
            depths: 1d numpy array
                Depth values
            method : {‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’}, optional
                Interpolation method. Default is quadratic.
    """
    def __init__(self, values, depths, method='quadratic'):
        
        if np.ndim(values) == 3:
            values = values[0,0,:]

        # initialize 1d interpolator
        self.interp = interp1d(x=depths, y=values, kind=method, fill_value="extrapolate")

        # store interpolation data
        self.depth_nodes = depths
        self.values = values

    def get_nodes(self):
        return (self.values, self.depth_nodes)

    def eval_xy(self, x, y, z, grid=False):
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
                   Depth of the positions(s) where the interpolation is to be evaluated
                grid: bool
                   Specify how to combine elements of x,y,z.

            Returns:
                vi: Interpolated values
        """
        v = self.eval_ll(lat=y, lon=x, z=z, grid=grid, squeeze=False)

        if np.ndim(v) == 3:
            v = np.swapaxes(v, 0, 1)

        v = np.squeeze(v)
        return v

    def eval_ll(self, lat, lon, z, grid=False, squeeze=True):
        """ Interpolate using spherical coordinate system (lat-lon).

            lat,lot,z can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the coordinates (lat_i, lon_i, z_i), where lat=(lat_1,...,lat_N), 
            lon=(lon_1,...,lon_N) and z=(z_,...,z_K). Note that in this case, lat and 
            lon must have the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (lat_i, lon_j, z_k), where lat=(lat_1,...,lat_N), 
            lon=(lon_1,...,lon_M) and z=(z_1,...,z_K). Note that in this case, the lengths 
            of lat and lon do not have to be the same.

            Args: 
                lat: float or array
                    Latitude of the positions(s) where the bathymetry is to be evaluated
                lon: float or array
                    Longitude of the positions(s) where the bathymetry is to be evaluated
                z: float or array
                    Depth of the positions(s) where the interpolation is to be evaluated
                grid: bool
                    Specify how to combine elements of lat,lon,z.

            Returns:
                zi: Interpolated bathymetry values (or derivates)
        """
        if np.ndim(lat) == 0: lat = np.array([lat])
        if np.ndim(lon) == 0: lon = np.array([lon])
        if np.ndim(z) == 0: z = np.array([z])

        if grid:
            s = (len(lat), len(lon), len(z))

        else:
            assert len(lat) == len(lon) == len(z), 'when grid is False, lat,lon,z must have the same length'

            s = len(lat)

        v = self.interp(z)

        if grid:
            v = np.ones(s) * v[np.newaxis, np.newaxis, :]
        
        if squeeze:
            v = np.squeeze(v)

        return v
