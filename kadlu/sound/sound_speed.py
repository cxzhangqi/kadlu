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
from kadlu.geospatial.interpolation import Interpolator2D, Interpolator3D, Uniform3D, DepthInterpolator3D


class SoundSpeed():
    """ Class for handling computation and interpolation of sound speed. 

        The sound speed can be specified via the argument ssp (sound speed 
        profile) or computed from bathymetry, temperature and salinity data 
        passed via the ocean argument.

        ssp can be either a single value, in which case the sound speed is 
        the same everywhere, or a tuple (c,z) where c is an array of sound 
        speed values and z is an array of depths.

        The eval method is used to obtain the interpolated sound speed at 
        any 3D point(s) in space.

        Args:
            ocean: DataProvider
                Provider of environmental data, including bathymetry, 
                temperature and salinity.
            ssp: float or tuple
                Sound speed profile. May be specified either as a float, 
                in which case the sound speed is the same everywhere, or as 
                a tuple (c,z) where c is an array of sound speeds and z is 
                an array of depth values.
            xy_res: float
                Horizontal (xy) resolution of the interpolation grid in meters.
                The default value is 1000 meters.
            num_depths: int
                Number of depth values for the interpolation grid. If rel_err is specified, 
                num_depths becomes the maximum allowed number of depth values. 
                The default value is 50.
            rel_err: float
                Maximum deviation of the interpolation, expressed as a ratio of the 
                range of sound-speed values. The default value is 0.001.

        Attributes: 
            self.interp: instance of DepthInterpolator3D, Uniform3D, or Interpolator3D
                Sounds speed interpolation function
            self.origin: LatLon
                Origin of the x-y planar coordinate system.
                None, if ssp is specified.
    """
    def __init__(self, ocean=None, ssp=None, xy_res=1000, num_depths=50, rel_err=1E-3):

        assert ocean is not None or ssp is not None, "ocean or ssp must be specified"

        if ssp is not None:
            self.origin = None           
            self.data = ssp 

            if isinstance(ssp, tuple):
                self.interp = DepthInterpolator3D(values=ssp[0], depths=ssp[1])

            else:
                self.interp = Uniform3D(value=ssp)

        else:
            self.origin = ocean.origin

            # convert from meters to degrees
            lat_res = 1./deg2rad * xy_res * DLDL_over_DXDY(lat=self.origin.latitude, lat_deriv_order=1, lon_deriv_order=0)
            lon_res = 1./deg2rad * xy_res * DLDL_over_DXDY(lat=self.origin.latitude, lat_deriv_order=0, lon_deriv_order=1)

            # geographic boundaries
            S = ocean.SW.latitude
            N = ocean.NE.latitude
            W = ocean.SW.longitude
            E = ocean.NE.longitude

            # lat and lon coordinates
            num_lats = int(np.ceil((N - S) / lat_res)) + 1
            lats = np.linspace(S, N, num=num_lats)
            num_lons = int(np.ceil((E - W) / lon_res)) + 1
            lons = np.linspace(W, E, num=num_lons)

            # generate depth coordinates
            depths = self._depth_coordinates(ocean, lats, lons, num_depths=num_depths, rel_err=rel_err)

            # temperature and salinity
            t = ocean.temp(lat=lats, lon=lons, z=depths, grid=True)
            s = ocean.salinity(lat=lats, lon=lons, z=depths, grid=True)

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

            # store interpolation data
            self.data = (c, lats, lons, depths)


    def _depth_coordinates(self, ocean, lats, lons, num_depths, rel_err):

        seafloor_depth = -ocean.bathy(lat=lats, lon=lons, grid=True)

        # find deepest point
        deepest_point = np.unravel_index(np.argmax(seafloor_depth), seafloor_depth.shape)

        # depth and lat,lon coordinates at deepest point
        max_depth = seafloor_depth[deepest_point]
        lat = lats[deepest_point[0]]
        lon = lons[deepest_point[1]]

        # compute temperature, salinity and sound speed for every 1 meter
        z = -np.arange(0, int(np.ceil(max_depth))+1)
        t = ocean.temp(lat=lat, lon=lon, z=z, grid=True)
        s = ocean.salinity(lat=lat, lon=lon, z=z, grid=True)        
        c = self._sound_speed(lats=lat, lons=lon, z=z, t=t, SP=s)

        # determine grid
        indices, _ = interp_grid_1d(y=c, x=z, num_pts=num_depths, rel_err=rel_err)
        depths = -z[indices]

        return depths


    def eval(self, x=None, y=None, z=None, grid=False, geometry='planar'):
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

            If x,y,z are not specified, the method returns the underlying 
            sound-speed data on which the interpolation is performed, either 
            as a (c,lat,lon,z) tuple, or as a float if the sound speed is 
            the same everywhere.

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
                c: Interpolated sound speed values
        """
        if x is None and y is None and z is None:
            c = self.data

        else:
            if geometry == 'planar':
                c = self.interp.eval_xy(x=x, y=y, z=z, grid=grid)

            elif geometry == 'spherical':
                c = self.interp.eval_ll(lat=y, lon=x, z=z, grid=grid)

        return c


    def _sound_speed(self, lats, lons, z, t, SP):
        """ Compute sound speed from temperature, salinity, and depth 
            for a given latitude and longitude.

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