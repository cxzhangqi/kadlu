""" Ocean module within the kadlu package

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
from kadlu.geospatial.interpolation import Interpolator2D, Interpolator3D, Uniform2D, Uniform3D


class Ocean():
    """ Class for handling ocean data requests.

        TODO: Implement loading of temp, salinity and wave data.

        Args:
            lats: 1d numpy array
                Latitude coordinates used for interpolation. 
            lons: 1d numpy array
                Longitude coordinates used for interpolation.
            depths: 1d numpy array
                Depth coordinates used for interpolation.

        Attributes: 

    """
    def __init__(self, bathy=None, temp=None, salinity=None, wave=None, wave_var=None,\
                south=-90, north=90, west=-180, east=180, time=None,\
                lat_ref=None, lon_ref=None, water_density=1.0):

        self.water_density = water_density

        # origo of x-y coordinate system
        if lat_ref is None: lat_ref = 0.5 * (south + north)
        if lon_ref is None: lon_ref = 0.5 * (west + east)
        self.origin = LatLon(lat_ref, lon_ref)

        # save south-west and north-east corners as class attributes
        self.SW = LatLon(south, west)
        self.NE = LatLon(north, east)

        # load data and initialize interpolation tables
        self._init_bathy(bathy, south, north, west, east)
        self._init_temp(temp, south, north, west, east, time)
        self._init_salinity(salinity, south, north, west, east, time)
        self._init_wave(wave, south, north, west, east, time)


    def _init_bathy(self, bathy, south, north, west, east):

        if bathy is None:
            self.bathy_data = None
            self.bathy_interp = None

        elif isinstance(bathy, str):
    
            if bathy == "CHS":

                # load data
                self.bathy_data = chs.load(south, north, west, east)

                # lat coordinates
                num_lats = int(np.ceil((north - south) / 0.001)) + 1
                lats = np.linspace(south, north, num=num_lats)

                # lon coordinates
                num_lons = int(np.ceil((east - west) / 0.001)) + 1
                lons = np.linspace(west, east, num=num_lons)

                # interpolate
                self.bathy_interp = Interpolator2D(values=self.bathy_data[0],\
                        lats=self.bathy_data[1], lons=self.bathy_data[2], origin=self.origin,\
                        method_irreg='regular', lats_reg=lats, lons_reg=lons)       

            elif bathy == "GEBCO":

                # load data
                self.bathy_data = gebco.load(south, north, west, east)

                # interpolate
                self.bathy_interp = Interpolator2D(values=self.bathy_data[0],\
                        lats=self.bathy_data[1], lons=self.bathy_data[2], origin=self.origin)       

            else: 
                print('Error: Unknown bathymetry source {0}.'.format(bathy))
                exit(1)

        elif isinstance(bathy, tuple):
            self.bathy_data = bathy
            self.bathy_interp = Interpolator2D(values=self.bathy_data[0],\
                    lats=self.bathy_data[1], lons=self.bathy_data[2], origin=self.origin)       

        else:
            self.bathy_data = bathy
            self.bathy_interp = Uniform2D(bathy)


    def _init_temp(self, temp, south, north, west, east, time):

        if temp is None:
            self.temp_data = None
            self.temp_interp = None

        elif isinstance(temp, str):
    
            if temp == "HYCOM":
                print('Error: Data loading from {0} not yet implemented.'.format(temp))

            else: 
                print('Error: Unknown temperature source {0}.'.format(temp))
                exit(1)

        elif isinstance(temp, tuple):
            self.temp_data = temp
            self.temp_interp = Interpolator3D(values=self.temp_data[0], lats=self.temp_data[1],\
                    lons=self.temp_data[2], depths=self.temp_data[3],\
                    origin=self.origin, method='linear')

        else:
            self.temp_data = temp
            self.temp_interp = Uniform3D(temp)


    def _init_salinity(self, salinity, south, north, west, east, time):

        if salinity is None:
            self.salinity_data = None
            self.salinity_interp = None

        elif isinstance(salinity, str):
    
            if salinity == "HYCOM":
                print('Error: Data loading from {0} not yet implemented.'.format(salinity))

            else: 
                print('Error: Unknown salinity source {0}.'.format(salinity))
                exit(1)

        elif isinstance(salinity, tuple):
            self.salinity_data = salinity
            self.salinity_interp = Interpolator3D(values=self.salinity_data[0], lats=self.salinity_data[1],\
                    lons=self.salinity_data[2], depths=self.salinity_data[3],\
                    origin=self.origin, method='linear')

        else:
            self.salinity_data = salinity
            self.salinity_interp = Uniform3D(salinity)


    def _init_wave(self, wave, south, north, west, east, time):

        if wave is None:
            self.wave_data = None
            self.wave_interp = None

        elif isinstance(wave, str):
    
            if wave == "ERA5":
                print('Error: Data loading from {0} not yet implemented.'.format(wave))

            elif wave == "RDWPS":
                print('Error: Data loading from {0} not yet implemented.'.format(wave))

            elif wave == "WWIII":
                print('Error: Data loading from {0} not yet implemented.'.format(wave))

            else: 
                print('Error: Unknown wave source {0}.'.format(wave))
                exit(1)

        elif isinstance(wave, tuple):
            self.wave_data = wave
            self.wave_interp = Interpolator2D(values=self.wave_data[0],\
                    lats=self.wave_data[1], lons=self.wave_data[2], origin=self.origin)       

        else:
            self.wave_data = wave
            self.wave_interp = Uniform2D(wave)


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

            If x and y are not specified, the method returns the underlying 
            bathymetric data on which the interpolation is performed, either 
            as a (bathy,lat,lon) tuple, or as a float if the depth is the same 
            everywhere.

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
        if x is None and y is None:
            z = self.bathy_data
        
        else:
            if geometry == 'planar':
                z = self.bathy_interp.eval_xy(x=x, y=y, grid=grid)

            elif geometry == 'spherical':
                z = self.bathy_interp.eval_ll(lat=y, lon=x, grid=grid)

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
            grad = self.bathy_interp.eval_xy(x=x, y=y, grid=grid, x_deriv_order=deriv_order[0], y_deriv_order=deriv_order[1])

        elif geometry == 'spherical':
            grad = self.bathy_interp.eval_ll(lat=y, lon=x, grid=grid, lat_deriv_order=deriv_order[1], lon_deriv_order=deriv_order[0])

        return grad


    def temp(self, x=None, y=None, z=None, grid=False, geometry='planar'):
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

            If x,y,z are not specified, the method returns the underlying 
            temperature data on which the interpolation is performed, either 
            as a (temp,lat,lon,z) tuple, or as a float if the temperature is 
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
                t: Interpolated temperature values
        """
        if x is None and y is None and z is None:
            t = self.temp_data

        else:
            if geometry == 'planar':
                t = self.temp_interp.eval_xy(x=x, y=y, z=z, grid=grid)

            elif geometry == 'spherical':
                t = self.temp_interp.eval_ll(lat=y, lon=x, z=z, grid=grid)

        return t


    def salinity(self, x=None, y=None, z=None, grid=False, geometry='planar'):
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

            If x,y,z are not specified, the method returns the underlying 
            salinity data on which the interpolation is performed, either 
            as a (salinity,lat,lon,z) tuple, or as a float if the salinity is 
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
                s: Interpolated salinity values
        """
        if x is None and y is None and z is None:
            s = self.salinity_data

        if geometry == 'planar':
            s = self.salinity_interp.eval_xy(x=x, y=y, z=z, grid=grid)

        elif geometry == 'spherical':
            s = self.salinity_interp.eval_ll(lat=y, lon=x, z=z, grid=grid)

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

            If x and y are not specified, the method returns the underlying 
            wave data on which the interpolation is performed, either 
            as a (wave,lat,lon) tuple, or as a float if the wave data is the same 
            everywhere.

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
        if x is None and y is None:
            w = self.wave_data

        else:
            if geometry == 'planar':
                w = self.wave_interp.eval_xy(x=x, y=y, grid=grid)

            elif geometry == 'spherical':
                w = self.wave_interp.eval_ll(lat=y, lon=x, grid=grid)

        return w


