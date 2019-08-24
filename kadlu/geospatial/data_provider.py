""" Geospatial data provider module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import numpy as np
from kadlu.utils import LatLon
import kadlu.geospatial.data_sources.chs as chs
import kadlu.geospatial.data_sources.gebco as gebco 
from kadlu.geospatial.interpolation import Interpolator2D, Interpolator3D


class DataProvider():
    """ Class for handling geospatial data requests. 

        TODO: Get rid of the storage_location argument and instead use the config.ini file.

        TODO: Implement loading of temp, salinity and wave data.

        Args:
            lat_bin_size: float
                Latitude bin size used for interpolation grid. 
                The default values is 0.001 degrees (110 meters)
            lon_bin_size: float
                Longitude bin size used for interpolation grid
                The default values is 0.001 degrees (110 meters at the Equator)

        Attributes: 

    """
    def __init__(self, storage_location, bathy_source=None, temp_source=None,\
                salinity_source=None, wave_source=None, wave_var=None,\
                south=-90, north=90, west=-180, east=180, time=None,\
                lat_ref=None, lon_ref=None, lat_bin_size=0.001, lon_bin_size=0.001):

        self.bathy_data = None
        self.temp_data = None
        self.salinity_data = None
        self.wave_data = None

        self.bathy_interpolator = None
        self.temp_interpolator = None
        self.salinity_interpolator = None
        self.wave_interpolator = None

        # load bathymetry data
        if bathy_source == "CHS":
            self.bathy_data = chs.load(storage_location, south, north, west, east)

        elif bathy_source == "GEBCO":
            self.bathy_data = gebco.load(storage_location, south, north, west, east)

        elif bathy_source is not None:
            print('Warning: Unknown bathymetry data source {0}.'.format(bathy_source))
            print('Proceeding without bathymetry data')

        # load temperature and salinity
        # TODO: replace this with calls to the actual load methods
        if self.bathy_data is None:
            max_depth = 10000
        else:
            max_depth = -np.min(self.bathy_data[0])
        print('\n*** generating dummy temperature data ***')
        self.temp_data = _generate_fake_data_3d(4, south, north, west, east, max_depth)
        print('*** generating dummy salinity data *** ')
        self.salinity_data = _generate_fake_data_3d(35, south, north, west, east, max_depth)
        print('*** generating dummy wave data *** ')
        self.wave_data = _generate_fake_data_2d(1, south, north, west, east)

        # reference coordinates for x-y coordinate system
        if lat_ref is None:
            lat_ref = (south + north) / 2

        if lon_ref is None:
            lon_ref = (west + east) / 2

        self.origin = LatLon(lat_ref, lon_ref)

        # coordinates of regular lat-lon grid for 2D interpolation
        num_lats = int(np.ceil((north - south) / lat_bin_size)) + 1
        num_lons = int(np.ceil((east - west) / lon_bin_size)) + 1
        lats = np.linspace(south, north, num=num_lats)
        lons = np.linspace(west, east, num=num_lons)

        print(lats)
        print(lons)
        print(lat_ref, lon_ref)

        # initialize bathymetry interpolation table
        if self.bathy_data is not None:    
            self.bathy_interpolator = Interpolator2D(values=self.bathy_data[0],\
                    lats=self.bathy_data[1], lons=self.bathy_data[2], origin=self.origin,\
                    method_irreg='regular', lats_reg=lats, lons_reg=lons)       

        # initialize temperature interpolation table
        if self.temp_data is not None:    
            self.temp_interpolator = Interpolator3D(values=self.temp_data[0],\
                    lats=self.temp_data[1], lons=self.temp_data[2], depths=self.temp_data[3],\
                    origin=self.origin, method='linear')       

        # initialize salinity interpolation table
        if self.salinity_data is not None:    
            self.salinity_interpolator = Interpolator3D(values=self.salinity_data[0],\
                    lats=self.salinity_data[1], lons=self.salinity_data[2], depths=self.salinity_data[3],\
                    origin=self.origin, method='linear')       

        # initialize wave interpolation table
        if self.wave_data is not None:    
            self.wave_interpolator = Interpolator2D(values=self.wave_data[0],\
                    lats=self.wave_data[1], lons=self.wave_data[2], origin=self.origin,\
                    method_irreg='regular', lats_reg=lats, lons_reg=lons)       


    def bathy(self, x, y, grid=False, geometry='planar'):
        """ Evaluate interpolated bathymetry in spherical (lat-lon) or  
            planar (x-y) geometry.

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
                    x-coordinate(s) or longitude(s)
                y: float or array
                    y-coordinate(s) or latitude(s)
                grid: bool
                    Specify how to combine elements of x and y. If x and y have different
                    lengths, specifying grid has no effect as it is automatically set to True.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'

            Returns:
                z: Interpolated bathymetry values
        """
        if geometry == 'planar':
            z = self.bathy_interpolator.eval_xy(x=x, y=y, grid=grid)

        elif geometry == 'spherical':
            z = self.bathy_interpolator.eval_ll(lat=y, lon=x, grid=grid)

        return z


    def bathy_gradient(self, x=None, y=None, axis='x', grid=False, geometry='planar'): 
        """ Evaluate interpolated bathymetry gradient in spherical (lat-lon) or  
            planar (x-y) geometry, along the specified axis.

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
                   x-coordinate(s) or longitude(s)
                y: float or array
                   y-coordinate(s) or latitude(s)
                grid: bool
                   Specify how to combine elements of x and y.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'
                axis: str
                    Axis along which gradient is computed. Can be either 'x' (default) or 'y'

            Returns:
                grad: Interpolated bathymetry gradient values
        """
        deriv_order = [(axis=='x'), (axis!='x')]

        if geometry == 'planar':                
            grad = self.bathy_interpolator.eval_xy(x=x, y=y, grid=grid, x_deriv_order=deriv_order[0], y_deriv_order=deriv_order[1])

        elif geometry == 'spherical':
            grad = self.bathy_interpolator.eval_ll(lat=y, lon=x, grid=grid, lat_deriv_order=deriv_order[1], lon_deriv_order=deriv_order[0])

        return grad


    def temp(self, x, y, z, grid=False, geometry='planar'):
        """ Evaluate interpolated temperature in spherical (lat-lon) or  
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
                t: Interpolated temperature values
        """
        if geometry == 'planar':
            t = self.temp_interpolator.eval_xy(x=x, y=y, z=z, grid=grid)

        elif geometry == 'spherical':
            t = self.temp_interpolator.eval_ll(lat=y, lon=x, z=z, grid=grid)

        return t


    def salinity(self, x, y, z, grid=False, geometry='planar'):
        """ Evaluate interpolated salinity in spherical (lat-lon) or  
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
                s: Interpolated salinity values
        """
        if geometry == 'planar':
            s = self.salinity_interpolator.eval_xy(x=x, y=y, z=z, grid=grid)

        elif geometry == 'spherical':
            s = self.salinity_interpolator.eval_ll(lat=y, lon=x, z=z, grid=grid)

        return s


    def wave(self, x, y, grid=False, geometry='planar'):
        """ Evaluate interpolated wave data in spherical (lat-lon) or  
            planar (x-y) geometry.

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
                    x-coordinate(s) or longitude(s)
                y: float or array
                    y-coordinate(s) or latitude(s)
                grid: bool
                    Specify how to combine elements of x and y. If x and y have different
                    lengths, specifying grid has no effect as it is automatically set to True.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'

            Returns:
                w: Interpolated wave data
        """
        if geometry == 'planar':
            w = self.wave_interpolator.eval_xy(x=x, y=y, grid=grid)

        elif geometry == 'spherical':
            w = self.wave_interpolator.eval_ll(lat=y, lon=x, grid=grid)

        return w


def _generate_fake_data_2d(value, south, north, west, east):
    N = 30
    lats = (np.arange(N) + 0.5) / N * (north - south) + south 
    lons = (np.arange(N) + 0.5) / N * (east - west) + west 
    shape = (len(lats), len(lons))
    values = value * np.ones(shape)
    return (values, lats, lons)
    
def _generate_fake_data_3d(value, south, north, west, east, max_depth):
    N = 30
    lats = (np.arange(N) + 0.5) / N * (north - south) + south 
    lons = (np.arange(N) + 0.5) / N * (east - west) + west 
    depths = 2 * np.arange(N) / (N-1) * max_depth
    shape = (len(lats), len(lons), len(depths))
    values = value * np.ones(shape)
    return (values, lats, lons, depths)