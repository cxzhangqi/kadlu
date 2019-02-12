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
from scipy.interpolate import RectBivariateSpline, RectSphereBivariateSpline
from pyost.bathy_reader import BathyReader, LatLon
from pyost.util import deg2rad, XYtoLL, LLtoXY
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MaxNLocator
from scipy.interpolate import griddata


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
        lat, lon = XYtoLL(x, y, lat_ref=self.origin.latitude, lon_ref=self.origin.longitude)

        zi = self.eval_ll(lat=lat, lon=lon, grid=grid)

#        zi = self.interp_xy.__call__(x=x, y=y, grid=grid)

#        if np.ndim(zi) == 0:
#            zi = float(zi)

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

        if np.ndim(zi) == 0:
            zi = float(zi)

        return zi

    def _torad(self, lat, lon):
        """ Convert latitute and longitude values from degrees to radians.

            The method expects the latitude to be in the range (-90,90) and
            the longitude to be in the range (-180,180).

            The output latitude is in the range (0,pi) and the output 
            longitude is in the range (-pi,pi).

            Args: 
                lat: float or array
                   latitude(s) in degrees from -90 to +90.
                lon: float or array
                   longitude(s) in degrees from -180 to +180.

            Returns:
                lat_rad: float or array
                    latitude(s) in radians from 0 to pi.
                lon_rad: float or array
                    longitude(s) in radians from -pi to +pi.
        """
        lat_rad = (lat + 90) * deg2rad
        lon_rad = lon * deg2rad
        return lat_rad, lon_rad

    def plot_ll(self):        
        """ Draw a map of the elevation in polar coordinates.

            Returns:
                fig: matplotlib.figure.Figure
                    A figure object.
        """
        fig = self._plot(ll=True)
        return fig

    def plot_xy(self):        
        """ Draw a map of the elevation in planar coordinates.

            Returns:
                fig: matplotlib.figure.Figure
                    A figure object.
        """
        fig = self._plot(ll=False)
        return fig

    def _plot(self, ll=True):
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
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        if ll:
            X = self.lon_nodes
            Y = self.lat_nodes
            X,Y = np.meshgrid(X,Y,indexing='ij')
            Z = self.bathy_ll
            Z = np.swapaxes(Z, 0, 1)
        else:
            X, Y = LLtoXY(lat=self.lon_nodes, lon=self.lon_nodes)
            X,Y = np.meshgrid(X,Y,indexing='ij')
            Z = self.bathy_ll

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

        # plot the surface
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False)

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
            ax.set_zlabel('Elevation (km)')
        else:
            ax.set_zlabel('Elevation (m)')

        # Customize the z axis
        ax.set_zlim(zmin, zmax)
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.3, aspect=5)

        return fig

#    def plot_xy(self):
