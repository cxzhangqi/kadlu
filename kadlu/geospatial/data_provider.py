""" Geospatial data provider module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import kadlu.geospatial.data_sources.chs as chs
import kadlu.geospatial.data_sources.gebco as gebco 
from kadlu.geospatial.interpolation import Interpolator2D


class DataProvider():
    """ Class for handling geospatial data requests. 

        TODO: Get rid of the storage_location argument and instead use the config.ini file

        Args:

        Attributes: 

    """
    def __init__(self, storage_location, bathy_source=None, temp_source=None,\
                salinity_source=None, wave_source=None, wave_var=None,\
                south=-90, north=90, west=-180, east=180, time=None,\
                origin=None, interpolation_method='cubic'):

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

        # initialize bathymetry interpolation table
        if self.bathy_data is not None:
            self.bathy_interpolator = Interpolator2D(values=self.bathy_data[0],\
                    lats=self.bathy_data[1], lons=self.bathy_data[2], origin=origin, method=interpolation_method)       
        

    def bathy(self, x=None, y=None, grid=False, geometry='planar'):
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
                   Specify how to combine elements of x and y.
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

    def bathy_gradient(self, axis='x', x=None, y=None, grid=False, geometry='planar'): 
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
