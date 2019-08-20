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

""" Bathymetry interpolation module within the kadlu library

    This module handles two-dimensional interpolation of bathymetry data in 
    spherical and planar geometry, on regular and iregular grids.

    Contents:
        GridData class:
        BathyInterpolator class
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline, RectSphereBivariateSpline
from kadlu.geospatial.bathy_reader import BathyReader, LatLon
from kadlu.utils import deg2rad, XYtoLL, LLtoXY, torad, DLDL_over_DXDY
from scipy.interpolate import griddata

from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from matplotlib import pyplot as plt


class GridData():
    """ Wrapper function around scipy's interpolate.griddata

        https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.interpolate.griddata.html

        Attributes: 
            x: float or 1d array
                data points 1st coordinate
            y: float or 1d array
                data points 2nd coordinate
            z: float or 1d array
                data values
            method : {‘linear’, ‘nearest’, ‘cubic’}, optional
    """
    def __init__(self, x, y, z, method='cubic'):
        self.xy = np.column_stack((x,y))
        self.z = z
        self.method = method

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
                zi: Interpolated values
        """        
        if grid:
            M = len(theta)
            N = len(phi)
            theta, phi = np.meshgrid(theta, phi)
            theta = np.reshape(theta, newshape=(M*N))
            phi = np.reshape(phi, newshape=(M*N))        

        xi = np.column_stack((theta, phi))

        zi = griddata(self.xy, self.z, xi, method=self.method)

        if grid:
            zi = np.reshape(zi, newshape=(M,N))

        return zi


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
            method : {‘linear’, ‘nearest’, ‘cubic’}, optional
                Interpolation method used for unstructured data (GeoTIFF and XYZ)
    """
    def __init__(self, bathy_reader, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180), origin=None, method='cubic'):
        
        # read bathymetry data from file
        lat, lon, bathy = bathy_reader.read(latlon_SW, latlon_NE)

        # check that data was 
        assert len(lat) > 0, "Reader unable to retrieve any bathymetry data for selected region"

        # compute coordinates of origin, if not provided
        if origin is None:
            lat_ref = (np.min(lat) + np.max(lat)) / 2
            lon_ref = (np.min(lon) + np.max(lon)) / 2
            origin = LatLon(lat_ref, lon_ref)

        self.origin = origin

        # check if bathymetry data are on a regular or unstructured grid
        reggrid = (np.ndim(bathy) == 2)

        # convert to radians
        lat_rad, lon_rad = torad(lat, lon)

        # initialize lat-lon interpolator
        if reggrid:
            self.interp_ll = RectSphereBivariateSpline(u=lat_rad, v=lon_rad, r=bathy)
        else:
            self.interp_ll = GridData(x=lat_rad, y=lon_rad, z=bathy, method=method)

        # store grids
        self.lat_nodes = lat
        self.lon_nodes = lon
        self.bathy = bathy

    def eval_xy(self, x, y, grid=False, x_deriv_order=0, y_deriv_order=0):
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
                x_deriv_order: int
                    Order of x-derivative
                y_deriv_order: int
                    Order of y-derivative

            Returns:
                zi: Interpolated bathymetry values
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

        if np.ndim(self.bathy) == 1 and np.ndim(zi) == 2:
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
            xo, yo = self._ll_to_xy(origin.latitude, origin.longitude)

        # distance array
        dr = distance / float(bins)
        r = np.arange(bins, dtype=np.float)
        r *= dr
        r += 0.5 * dr

        bathy = list()

        # loop over angles
        a = angle
        da = 360. / float(num_slices)
        for _ in range(num_slices):
            x = r * np.cos(a * np.pi / 180.)
            y = r * np.sin(a * np.pi / 180.)
            x += xo
            y += yo
            b = self.eval_xy(x=x, y=y)
            bathy.append(b)
            a += da

        if num_slices == 1:
            bathy = bathy[0]
        
        return bathy

    def _ll_to_xy(self, lat, lon, lat_ref=None, lon_ref=None, grid=False):

        if lat_ref is None:
            lat_ref = self.origin.latitude
        if lon_ref is None:
            lon_ref = self.origin.longitude

        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=lat_ref, lon_ref=lon_ref, grid=grid)
        return x, y

    def _xy_to_ll(self, x, y, lat_ref=None, lon_ref=None, grid=False):

        if lat_ref is None:
            lat_ref = self.origin.latitude
        if lon_ref is None:
            lon_ref = self.origin.longitude

        lat, lon = XYtoLL(x=x, y=y, lat_ref=lat_ref, lon_ref=lon_ref, grid=grid)
        return lat, lon

    def plot_ll(self, lat_bins=100, lat_min=None, lat_max=None, lon_bins=100, lon_min=None, lon_max=None):        
        """ Draw a map of the elevation in polar coordinates.

            Returns:
                fig: matplotlib.figure.Figure
                    A figure object.
        """
        fig = self._plot(True, lat_bins, lat_min, lat_max, lon_bins, lon_min, lon_max)
        return fig

    def plot_xy(self, x_bins=100, x_min=None, x_max=None, y_bins=100, y_min=None, y_max=None):        
        """ Draw a map of the elevation in planar coordinates.

            Returns:
                fig: matplotlib.figure.Figure
                    A figure object.
        """
        fig = self._plot(False, x_bins, x_min, x_max, y_bins, y_min, y_max)
        return fig

    def _plot(self, ll=True, x_bins=100, x_min=None, x_max=None, y_bins=100, y_min=None, y_max=None):
        """ Draw a map of the elevation using either polar or
            planar coordinates.

            Args:
                ll: bool
                    If ll=True use polar coordinates; otherwise 
                        use planar coordinates.

            Returns:
                fig: matplotlib.figure.Figure
                    A figure object.
        """        
        # interpolation grid
        if ll:
            x0 = self.lon_nodes
            y0 = self.lat_nodes
        else:
            grid = (np.ndim(self.bathy) == 2)
            x0, y0 = self._ll_to_xy(lat=self.lat_nodes, lon=self.lon_nodes, grid=grid)

        # axes ranges
        if x_min is None:
            x_min = np.min(x0)
        if x_max is None:
            x_max = np.max(x0)
        if y_min is None:
            y_min = np.min(y0)
        if y_max is None:
            y_max = np.max(y0)

        # binning
        dx = (x_max - x_min) / x_bins
        dy = (y_max - y_min) / y_bins

        # create axes
        X = np.arange(x_bins, dtype=np.float)
        X *= dx
        X += x_min + 0.5*dx
        Y = np.arange(y_bins, dtype=np.float)
        Y *= dy
        Y += y_min + 0.5*dy

        # interpolate bathymetry
        if ll:
            Z = self.eval_ll(lat=Y, lon=X, grid=True)
            Z = np.swapaxes(Z, 0, 1)
        else:
            Z = self.eval_xy(x=X, y=Y, grid=True)

        # mask NaN entries
        Z = np.ma.masked_invalid(Z)

        # meshgrid
        X,Y = np.meshgrid(X,Y)#,indexing='ij')

        # x and y axis range
        xrange = x_max - x_min
        yrange = y_max - y_min
        if xrange > 1e3:
            X = X / 1.e3
            x_max = x_max / 1.e3
            x_min = x_min / 1.e3
            xlabel = 'x (km)'
        else:
            xlabel = 'x (m)'

        if yrange > 1e3:
            Y = Y / 1.e3
            y_max = y_max / 1.e3
            y_min = y_min / 1.e3
            ylabel = 'y (km)'
        else:
            ylabel = 'y (m)'

        # z axis binning and range
        zmin = np.min(Z)
        zmax = np.max(Z)
        zrange = zmax - zmin
        p = int(np.log10(zrange)) - 1
        p = max(0, p)
        dz = pow(10,p)
        zmax = np.ceil(zmax / dz) * pow(10,p)
        zmin = np.floor(zmin / dz) * pow(10,p)

        if zrange > 1e3:
            zmax = zmax / 1.e3
            zmin = zmin / 1.e3
            Z = Z / 1.e3

        # plot
        fig, ax = plt.subplots(figsize=(8,6))
        img = ax.imshow(Z.T, aspect='auto', origin='lower', extent=(x_min, x_max, y_min, y_max))

        # axes titles
        if ll:
            ax.set_xlabel('Longitude (degrees east)')
            ax.set_ylabel('Latitude (degrees north)')
        else:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        if zrange > 1e3:
            zlabel = 'Elevation (km)'
        else:
            zlabel = 'Elevation (m)'

        # Add a color bar which maps values to colors
        fig.colorbar(img, format='%.02f', label=zlabel)

        return fig
