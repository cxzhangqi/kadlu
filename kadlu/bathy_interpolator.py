""" Bathymetry interpolation module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import numpy as np
from scipy.interpolate import RectBivariateSpline, RectSphereBivariateSpline
from kadlu.bathy_reader import BathyReader, LatLon
from kadlu.util import deg2rad, XYtoLL, LLtoXY, torad
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MaxNLocator
from scipy.interpolate import griddata


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

    def __call__(self, theta, phi, grid=False):
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

            Returns:
                zi: Interpolated values
        """        
        if grid:
            theta, phi = np.meshgrid(theta, phi)

        xi = np.column_stack((theta, phi))
        zi = griddata(self.xy, self.z, xi, method=self.method)

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

        # check if bathymetry data are on a regular or unstructured grid
        reggrid = (np.ndim(bathy) == 2)

        # convert to radians
        lat_rad, lon_rad = torad(lat, lon)

        # initialize lat-lon interpolator
        if reggrid:
            self.interp_ll = RectSphereBivariateSpline(u=lat_rad, v=lon_rad, r=bathy)
        else:
            self.interp_ll = GridData(x=lat_rad, y=lon_rad, z=bathy)

        # store grids
        self.lat_nodes = lat
        self.lon_nodes = lon
        self.bathy_ll = bathy

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
        lat, lon = XYtoLL(x, y, lat_ref=self.origin.latitude, lon_ref=self.origin.longitude, grid=grid)
        zi = self.eval_ll(lat=lat, lon=lon)
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
        lat_rad, lon_rad = torad(lat, lon)

        zi = self.interp_ll.__call__(theta=lat_rad, phi=lon_rad, grid=grid)

        if np.ndim(zi) == 0:
            zi = float(zi)

        return zi

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
            x0, y0 = LLtoXY(lat=self.lat_nodes, lon=self.lon_nodes, lat_ref=self.origin.latitude, lon_ref=self.origin.longitude, grid=True)

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
        X = np.arange(x_bins+1, dtype=np.float)
        X *= dx
        X += x_min
        Y = np.arange(y_bins+1, dtype=np.float)
        Y *= dy
        Y += y_min

        # interpolate bathymetry
        if ll:
            Z = self.eval_ll(lat=Y, lon=X, grid=True)
        else:
            Z = self.eval_xy(x=X, y=Y, grid=True)

        Z = np.swapaxes(Z, 0, 1)

        # meshgrid
        X,Y = np.meshgrid(X,Y,indexing='ij')

        # x and y axis range
        xrange = np.max(X) - np.min(X)
        yrange = np.max(Y) - np.min(Y)
        if xrange > 1e3:
            X = X / 1.e3
        if yrange > 1e3:
            Y = Y / 1.e3

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
        fig, ax = plt.subplots(figsize=(10,10))
        img = ax.imshow(Z.T, aspect='auto', origin='lower', extent=(np.min(X), np.max(X), np.min(Y), np.max(Y)))

        # axes titles
        if ll:
            ax.set_xlabel('Longitude (degrees east)')
            ax.set_ylabel('Latitude (degrees north)')
        else:
            if xrange > 1e3:
                ax.set_xlabel('X (km)')
            else:
                ax.set_xlabel('X (m)')
            
            if yrange > 1e3:
                ax.set_ylabel('Y (km)')
            else:
                ax.set_ylabel('Y (m)')

        if zrange > 1e3:
            zlabel = 'Elevation (km)'
        else:
            zlabel = 'Elevation (m)'

        # Add a color bar which maps values to colors
        fig.colorbar(img, format='%.02f', label=zlabel)

        return fig
