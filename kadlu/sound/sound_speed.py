""" Sound speed module within the kadlu package
"""
import gsw
import numpy as np
from kadlu.utils import interp_grid_1d, deg2rad
from kadlu.geospatial.interpolation import Interpolator2D, Interpolator3D, Uniform3D, DepthInterpolator3D


def sound_speed_teos10(lats, lons, z, t, SP):
    """ Compute sound speed from temperature, salinity, and depth 
        for a given latitude and longitude using the Thermodynamic 
        Equation of Seawater 2010.

        https://teos-10.github.io/GSW-Python/

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
    p = gsw.p_from_z(z=-z, lat=lats)  # sea pressure (gsw uses negative z below sea surface)
    SA = gsw.SA_from_SP(SP, p, lons, lats)  # absolute salinity
    CT = gsw.CT_from_t(SA, t, p)  # conservative temperature
    c = gsw.density.sound_speed(SA=SA, CT=CT, p=p)
    return c

class SoundSpeed():
    """ Class for handling computation and interpolation of sound speed. 

        The sound speed can be specified via the argument ssp (sound speed 
        profile) or computed from the ocean variables (temperature, salinity).

        ssp can be either a single value, in which case the sound speed is 
        the same everywhere, or a tuple (c,z) where c is an array of sound 
        speed values and z is an array of depths.

        The interp and interp_xyz method may be used to obtain the interpolated 
        sound speed at any set of coordinates.

        Args:
            ocean: instance of :class:`kadlu.geospatial.ocean.Ocean`
                Ocean variables
            ssp: float or tuple
                Sound speed profile. May be specified either as a float, 
                in which case the sound speed is the same everywhere, or as 
                a tuple (c,z) where c is an array of sound speeds and z is 
                an array of depth values.
            num_depths: int
                Number of depth values for the interpolation grid. The default value is 50.
            rel_err: float
                Maximum deviation of the interpolation, expressed as a ratio of the 
                range of sound-speed values. The default value is 0.001.
    """
    def __init__(self, ocean=None, ssp=None, num_depths=50, rel_err=1E-3):

        assert ocean is not None or ssp is not None, "ocean or ssp must be specified"

        if ssp is not None:
            if isinstance(ssp, tuple): self._interp = DepthInterpolator3D(values=ssp[0], depths=ssp[1])
            else: self._interp = Uniform3D(values=ssp)

        else:
            lat_res, lon_res = self._lat_lon_res(ocean, default_res=1.0) #default resolution is 1 degree, approx 100 km

            # geographic boundaries
            S,N,W,E = ocean.boundaries['south'], ocean.boundaries['north'], ocean.boundaries['west'], ocean.boundaries['east']

            # lat and lon coordinates
            num_lats = max(3, int(np.ceil((N - S) / lat_res)) + 1)
            lats = np.linspace(S, N, num=num_lats)
            num_lons = max(3, int(np.ceil((E - W) / lon_res)) + 1)
            lons = np.linspace(W, E, num=num_lons)

            # compute depth coordinates
            depths = self._depth_coordinates(ocean, lats, lons, num_depths=num_depths, rel_err=rel_err)

            # interpolate temperature and salinity on lat,lon,depth grid
            t = ocean.temp(lat=lats, lon=lons, depth=depths, grid=True)
            s = ocean.salinity(lat=lats, lon=lons, depth=depths, grid=True)

            # compute sound speed
            grid_shape = t.shape
            la, lo, de = np.meshgrid(lats, lons, depths)
            la = la.flatten()
            lo = lo.flatten()
            de = de.flatten()
            t = t.flatten()
            s = s.flatten()
            c = sound_speed_teos10(lats=la, lons=lo, z=de, t=t, SP=s)
            c = np.reshape(c, newshape=grid_shape)

            # create interpolator
            self._interp = Interpolator3D(values=c, lats=lats, lons=lons,\
                    depths=depths, origin=ocean.origin, method='linear')

    def _lat_lon_res(self, ocean, default_res):
        """ Determine lat,lon resolutions for interpolation grid

            Args:
                ocean: instance of :class:`kadlu.geospatial.ocean.Ocean`
                    Ocean variables
                default_res: float
                    Default resolution in degrees

            Returns:
                lat_res, lon_res: float,float
                    Resolutions in degrees.
        """
        temp_nodes = ocean.interps['temp'].get_nodes()
        temp_lat, temp_lon = temp_nodes[1], temp_nodes[1]
        salinity_nodes = ocean.interps['salinity'].get_nodes()
        salinity_lat, salinity_lon = salinity_nodes[1], salinity_nodes[1]
        temp_lat_res = default_res if temp_lat is None else np.abs(temp_lat[1] - temp_lat[0])
        temp_lon_res = default_res if temp_lon is None else np.abs(temp_lon[1] - temp_lon[0])
        salinity_lat_res = default_res if salinity_lat is None else np.abs(salinity_lat[1] - salinity_lat[0])
        salinity_lon_res = default_res if salinity_lon is None else np.abs(salinity_lon[1] - salinity_lon[0])
        lat_res = min(temp_lat_res, salinity_lat_res)
        lon_res = min(temp_lon_res, salinity_lon_res)
        return lat_res, lon_res

    def _depth_coordinates(self, ocean, lats, lons, num_depths, rel_err):
        """ Compute depth coordinates for lat,lon,depth interpolation grid.

            Args:
                ocean: instance of :class:`kadlu.geospatial.ocean`
                    Ocean variables
                lats,lons: numpy.array, numpy.array
                    Latitude and longitude coordinates
                num_depths: int
                    Number of depth values for the interpolation grid. The default value is 50.
                rel_err: float
                    Maximum deviation of the interpolation, expressed as a ratio of the 
                    range of sound-speed values. The default value is 0.001.

            Returns:
                depths: numpy.array
                    Depth coordinates
        """
        seafloor_depth = ocean.bathy(lat=lats, lon=lons, grid=True)

        # find deepest point
        deepest_point = np.unravel_index(np.argmax(seafloor_depth), seafloor_depth.shape)

        # depth and lat,lon coordinates at deepest point
        max_depth = seafloor_depth[deepest_point]
        lat = lats[deepest_point[0]]
        lon = lons[deepest_point[1]]

        # compute temperature, salinity and sound speed for every 1 meter
        z = np.arange(0, int(np.ceil(max_depth))+1)
        t = ocean.temp(lat=lat, lon=lon, depth=z, grid=True)
        s = ocean.salinity(lat=lat, lon=lon, depth=z, grid=True)        
        c = sound_speed_teos10(lats=lat, lons=lon, z=z, t=t, SP=s)

        # determine grid
        indices, _ = interp_grid_1d(y=c, x=z, num_pts=num_depths, rel_err=rel_err)
        depths = z[indices]

        return depths

    def interp(self, lat, lon, z, grid=False):
        """ Interpolate sound speed in spherical coordinates (lat-lon).

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
                : numpy.array
                    Interpolated sound speed values
        """
        return self._interp.interp(lat=lat, lon=lon, z=z, grid=grid)

    def interp_xy(self, x, y, z, grid=False):
        """ Interpolate sound speed in planar (x-y) geometry.

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
                    x-coordinates
                y: float or array
                    y-coordinates
                z: float or array
                    depth(s)
                grid: bool
                   Specify how to combine elements of x,y,z.

            Returns:
                : numpy.array
                    Interpolated sound speed values
        """
        return self._interp.interp_xy(x=x, y=y, z=z, grid=grid)


