""" Sound speed module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import gsw
import numpy as np
from kadlu.utils import LatLon, DLDL_over_DXDY, interp_grid_1d, deg2rad
from kadlu.geospatial.data_provider import DataProvider 
from kadlu.geospatial.interpolation import Interpolator2D, Interpolator3D


class SoundSpeed():
    """ Class for handling computation and interpolation of sound speed. 

        Args:

        Attributes: 

    """
    def __init__(self, env_data, xy_res=1000, num_depths=50, rel_err=1E-6):

        self.origin = env_data.origin

        # convert from meters to degrees
        lat_res = 1./deg2rad * xy_res * DLDL_over_DXDY(lat=self.origin.latitude, lat_deriv_order=1, lon_deriv_order=0)
        lon_res = 1./deg2rad * xy_res * DLDL_over_DXDY(lat=self.origin.latitude, lat_deriv_order=0, lon_deriv_order=1)

        # geographic boundaries
        S = env_data.SW.latitude
        N = env_data.NE.latitude
        W = env_data.SW.longitude
        E = env_data.NE.longitude

        # lat and lon coordinates
        num_lats = int(np.ceil((N - S) / lat_res)) + 1
        lats = np.linspace(S, N, num=num_lats)
        num_lons = int(np.ceil((E - W) / lon_res)) + 1
        lons = np.linspace(W, E, num=num_lons)

        # generate depth coordinates
        depths = self._depth_coordinates(env_data, lats, lons, num_depths=num_depths, rel_err=rel_err)

        # temperature and salinity
        t = env_data.temp(x=lons, y=lats, z=depths, geometry='spherical', grid=True)
        s = env_data.salinity(x=lons, y=lats, z=depths, geometry='spherical', grid=True)

        # sound speed
        grid_shape = t.shape
        la, lo, de = np.meshgrid(lats, lons, depths)
        la = la.flatten()
        lo = lo.flatten()
        de = de.flatten()
        t = t.flatten()
        s = s.flatten()
        c = self._sound_speed(lats=la, lons=lo, z=-de, t=t, SP=s)
        c = np.reshape(c, newshape=grid_shape)

        # create interpolator
        self.interp = Interpolator3D(values=c, lats=lats, lons=lons,\
                depths=depths, origin=self.origin, method='linear')

        self.lats = lats
        self.lons = lons
        self.depths = depths



    def _depth_coordinates(self, env_data, lats, lons, num_depths, rel_err):

        seafloor_depth = -env_data.bathy(x=lons, y=lats, grid=True, geometry='spherical')

        # find deepest point
        deepest_point = np.unravel_index(np.argmax(seafloor_depth), seafloor_depth.shape)

        # depth and lat,lon coordinates at deepest point
        max_depth = seafloor_depth[deepest_point]
        lat = lats[deepest_point[0]]
        lon = lons[deepest_point[1]]

        # compute temperature, salinity and sound speed for every 1 meter
        z = -np.arange(0, int(np.ceil(max_depth))+1)
        t = env_data.temp(x=lon, y=lat, z=z, geometry='spherical', grid=True)
        s = env_data.salinity(x=lon, y=lat, z=z, geometry='spherical', grid=True)        
        c = self._sound_speed(lats=lat, lons=lon, z=z, t=t, SP=s)

        # determine grid
        indices, _ = interp_grid_1d(y=c, x=z, num_pts=num_depths, rel_err=rel_err)
        depths = -z[indices]

        return depths


    def eval(self, x, y, z, grid=False, geometry='planar'):
        """ Evaluate interpolated sound speed in spherical (lat-lon) or  
            planar (x-y) geometry.

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
                    x-coordinate(s) or longitude(s)
                y: float or array
                    y-coordinate(s) or latitude(s)
                z: float or array
                    depth(s)
                grid: bool
                   Specify how to combine elements of x,y,z.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'

            Returns:
                s: Interpolated sound speed values
        """
        if geometry == 'planar':
            z = self.interp.eval_xy(x=x, y=y, z=z, grid=grid)

        elif geometry == 'spherical':
            z = self.interp.eval_ll(lat=y, lon=x, z=z, grid=grid)

        return z


    def _sound_speed(self, lats, lons, z, t, SP):
        """ Compute sound speed.

            Args:
                lats: numpy array
                    Latitudes (-90 to 90 degrees)
                lons: numpy array
                    Longitudes (-180 to 180 degrees)
                z: numpy array
                    Depths (meters)
                t: numpy array
                    In-situ temperature (Celsius)
                SP: numpy array
                    Practical Salinity (psu)

            Returns:
                c: numpy array
                    Sound speed (m/s) 
        """
        p = gsw.p_from_z(z=z, lat=lats)  # sea pressure
        SA = gsw.SA_from_SP(SP, p, lons, lats)  # absolute salinity
        CT = gsw.CT_from_t(SA, t, p)  # conservative temperature
        c = gsw.density.sound_speed(SA=SA, CT=CT, p=p)
        return c