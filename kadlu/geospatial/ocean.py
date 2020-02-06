import numpy as np
import pickle
from multiprocessing import Process
from kadlu.utils import LatLon
from kadlu.geospatial.interpolation             import      \
        Interpolator2D,                                     \
        Interpolator3D,                                     \
        Uniform2D,                                          \
        Uniform3D,                                          \
        interp_2D,                                          \
        interp_3D
from kadlu.geospatial.data_sources.data_util    import      \
        reshape_2D,                                         \
        reshape_3D,                                         \
        hash_key,                                           \
        bin_db,                                             \
        deserialize
from kadlu.geospatial.data_sources.chs          import Chs
#from kadlu.geospatial.data_sources.gebco        import Gebco
from kadlu.geospatial.data_sources.hycom        import Hycom
from kadlu.geospatial.data_sources.era5         import Era5
from kadlu.geospatial.data_sources.wwiii        import Wwiii

# dictionary for mapping strings to callback functions
fetch_map = dict(
        bathy_chs           = Chs().fetch_bathymetry,
        #bathy_gebco         = Gebco().fetch_bathymetry, 
        temp_hycom          = Hycom().fetch_temp,
        salinity_hycom      = Hycom().fetch_salinity,
        wavedir_era5        = Era5().fetch_wavedirection,
        wavedir_wwiii       = Wwiii().fetch_wavedirection,
        waveheight_era5     = Era5().fetch_windwaveswellheight,
        waveheight_wwiii    = Wwiii().fetch_windwaveheight,
        waveperiod_era5     = Era5().fetch_waveperiod,
        waveperiod_wwiii    = Wwiii().fetch_waveperiod,
        windspeed_era5      = Era5().fetch_wind,
        windspeed_wwiii     = Wwiii().fetch_wind
    )

load_map = dict(
        bathy_chs           = Chs().load_bathymetry,
        #bathy_gebco         = Gebco().load_bathymetry, 
        temp_hycom          = Hycom().load_temp,
        salinity_hycom      = Hycom().load_salinity,
        wavedir_era5        = Era5().load_wavedirection,
        wavedir_wwiii       = Wwiii().load_wavedirection,
        waveheight_era5     = Era5().load_windwaveswellheight,
        waveheight_wwiii    = Wwiii().load_windwaveheight,
        waveperiod_era5     = Era5().load_waveperiod,
        waveperiod_wwiii    = Wwiii().load_waveperiod,
        windspeed_era5      = Era5().load_wind,
        windspeed_wwiii     = Wwiii().load_wind
    )


def serialize_interp(interpfcn, reshapefcn, loadfcn, kwargs, seed):
    """ serialize an interpolation and store binary data to database.
        used for caching results
    
        interpfcn:
            callback function for interpolation
        reshapefcn:
            callback function for reshaping row data into matrix format
            for interpolation
        loadfcn:
            callback function for loading data to reshape
        kwargs:
            dictionary containing keyword arguments used to generate
            the object. this will be hashed and used as a database key
        seed:
            seed the hash to differentiate interpolation variables with
            the same set of kwargs
    """
    conn, db = bin_db()
    key = hash_key(kwargs, seed)
    db.execute('SELECT * FROM bin WHERE hash == ? LIMIT 1', (key,))
    if db.fetchone() is not None: return
    kwargs['var'] = seed.split('_')[1]
    obj = interpfcn(**reshapefcn(loadfcn, **kwargs))
    db.execute('INSERT INTO bin VALUES (?, ?)', (key, pickle.dumps(obj)))
    conn.commit()
    return 


def load_callback(**kwargs):
    """ bootstrap array data into callable for serialization when loading """
    v = kwargs['var'] + '_'
    data = np.array((kwargs[f'{v}val'], kwargs[f'{v}lat'], kwargs[f'{v}lon']))
    if f'{v}depth' not in kwargs.keys(): return data
    return np.append(data, kwargs[f'{v}depth'])


class Ocean():
    """ class for handling ocean data requests 

        data will be loaded using the given data sources and boundaries
        from arguments. an interpolation for each variable will be computed in 
        parallel

        any of the below load_args may also accept a callback function instead
        of a string or array value if you wish to write your own data loading
        function. the boundary arguments supplied here will be passed to the 
        callable, i.e. north, south, west, east, top, bottom, start, end

        callables or array arguments must be ordered by [val, lat, lon] for 2D 
        data, or [val, lat, lon, depth] for 3D data.

        it is also possible to supply a single float value for 

        args:
            cache:
                boolean. if True, resulting interpolations will be stored as 
                binary to be reused later. caching is True by default
            fetch:
                boolean. Description here ...
            default:
                boolean. If True, the default data source will be used if none 
                is specified. If False, zero values will be used if no data 
                source is specified. Default is True.   
            load_bathymetry: 
                source of bathymetry data. can be 'chs' to load previously 
                fetched data, or array ordered by [val, lat, lon], or float
            load_temp:
                source of temperature data. can be 'hycom' to load previously
                fetched data, or array ordered by [val, lat, lon, depth], or float
            load_salinity:
                source of salinity data. can be 'hycom' to load previously
                fetched data, or array ordered by [val, lat, lon, depth], or float
            load_wavedir:
                source of wave direction data. can be 'era5' or 'wwiii' to load
                previously fetched data, or array ordered by [val, lat, lon], or float
            load_waveheight:
                source of wave height data. can be 'era5' or 'wwiii' to load
                previously fetched data, or array ordered by [val, lat, lon], or float
            load_waveperiod:
                source of wave period data. can be 'era5' or 'wwiii' to load
                previously fetched data, or array ordered by [val, lat, lon], or float
            load_wind:
                source of wind speed data. can be 'era5' or 'wwiii' to load
                previously fetched data, or array ordered by [val, lat, lon], or float
            north, south:
                latitude boundaries (float)
            west, east:
                longitude boundaries (float)
            top, bottom:
                depth range in metres (float)
                only applies to salinity and temperature
            start, end:
                time range for data load query (datetime)
                if multiple times exist within range, they will be averaged
                before computing interpolation
            time:
                specify a single datetime as an alternative to using 
                the start, end kwargs. the nearest fetched time data 
                will be loaded 
    """
    def __init__(self, 
            cache           = True,
            fetch           = False,
            default         = True,
            load_bathymetry = 0,
            load_temp       = 0,
            load_salinity   = 0,
            load_wavedir    = 0,
            load_waveheight = 0,
            load_waveperiod = 0,
            load_windspeed  = 0,
            **kwargs):

        if 'start' in kwargs.keys() and 'end' in kwargs.keys(): 
            print('WARNING: data will be averaged over time frames for interpolation')
                #.\nto avoid this behaviour, use the \'time\' '
                #'keyword argument instead of start/end')

        self.cache = cache
        self.fetch = fetch
        self.interp = {}
        
        # set default data sources
        if default:
            if load_bathymetry == 0: load_bathymetry = 'chs'
            if load_temp == 0:       load_temp = 'hycom'
            if load_salinity == 0:   load_salinity = 'hycom'
            if load_wavedir == 0:    load_wavedir = 'era5'
            if load_waveheight == 0: load_waveheight = 'era5'
            if load_waveperiod == 0: load_waveperiod = 'era5'
            if load_windspeed == 0:  load_windspeed = 'era5'

        # collect load args for later use
        self.load_args = [load_bathymetry,
                          load_temp,
                          load_salinity,
                          load_wavedir,
                          load_waveheight,
                          load_waveperiod,
                          load_windspeed]

        # load ocean data
        self.load(**kwargs)


    def load(self, **kwargs):
        """ Load data for the specified region and time period.

            Args:
                north, south:
                    latitude boundaries (float)
                west, east:
                    longitude boundaries (float)
                top, bottom:
                    depth range in metres (float)
                    only applies to salinity and temperature
                start, end:
                    time range for data load query (datetime)
                    if multiple times exist within range, they will be averaged
                    before computing interpolation
                time:
                    specify a single datetime as an alternative to using 
                    the start, end kwargs. the nearest fetched time data 
                    will be loaded 
        """
        # assume water density to be 1.0 g/cm^3 everywhere
        self.water_density = 1.0

        # south-west and north-east corners of the region considered
        if 'south' not in kwargs.keys(): kwargs['south'] =  -90
        if 'north' not in kwargs.keys(): kwargs['north'] =   90
        if 'west'  not in kwargs.keys(): kwargs['west']  = -180
        if 'east'  not in kwargs.keys(): kwargs['east']  =  180
        self.SW = LatLon(kwargs['south'], kwargs['west'])
        self.NE = LatLon(kwargs['north'], kwargs['east'])

        # place origin of planar (x,y) coordinate system at the center
        lat_ref = 0.5 * (self.SW.latitude + self.NE.latitude)
        lon_ref = 0.5 * (self.SW.longitude + self.NE.longitude)
        self.origin = LatLon(lat_ref, lon_ref)

        # load data and create interpolation objects for each variable
        self.add_var_2D('bathy',      self.load_args[0], kwargs)
        self.add_var_3D('temp',       self.load_args[1], kwargs)
        self.add_var_3D('salinity',   self.load_args[2], kwargs)
        self.add_var_2D('wavedir',    self.load_args[3], kwargs)
        self.add_var_2D('waveheight', self.load_args[4], kwargs)
        self.add_var_2D('waveperiod', self.load_args[5], kwargs)
        self.add_var_2D('windspeed',  self.load_args[6], kwargs)


    def add_var_2D(self, var, load_arg, kwargs):
        """ Add 2D variable to the ocean.

            Args:
                var: str
                    Variable name
                load_arg: 
                    Data source. Can be str to load previously 
                    fetched data, or array ordered by [val, lat, lon], 
                    or float
        """
        self._add_var(var, load_arg, interp_2D, reshape_2D, kwargs)


    def add_var_3D(self, var, load_arg, kwargs):
        """ Add 3D variable to the ocean.

            Args:
                var: str
                    Variable name
                load_arg: 
                    Data source. Can be str to load previously 
                    fetched data, or array ordered by [val, lat, lon, depth], 
                    or float
        """
        self._add_var(var, load_arg, interp_3D, reshape_3D, kwargs)


    def _add_var(self, var, load_arg, interp_fcn, reshape_fcn, kwargs):
        """ Add variable to the ocean.

            Args:
                var: str
                    Variable name
                load_arg: 
                    Data source. Can be str to load previously 
                    fetched data, or array ordered by [val, lat, lon, depth], 
                    or float
                interp_fcn: 
                    Interpolation object.
                reshape_fcn:
                    Reshaping function.
        """
        if callable(load_arg): pass

        elif isinstance(load_arg, str):
            key = f'{var}_{load_arg.lower()}'
            if self.fetch == True: fetch_map[key](**kwargs)
            load_arg = load_map[key]

        # TODO: handling of tuples and flots/ints with dummy_load 
        #       and dummy_reshape functions is VERY clumsy. A more 
        #       elegant and simple implementation would be desirable.

        elif isinstance(load_arg, (list, tuple, np.ndarray)):
            if len(load_arg) not in (3, 4):
                raise ValueError(f'invalid array shape for load_{var}. '
                'arrays must be ordered by [val, lat, lon] for 2D data, or '
                '[val, lat, lon, depth] for 3D data')
            kwargs[f'{var}_val'] = load_arg[0]
            kwargs[f'{var}_lat'] = load_arg[1]
            kwargs[f'{var}_lon'] = load_arg[2]
            if len(load_arg) == 4: kwargs[f'{var}_depth'] = load_arg[3]
            #load_arg = load_callback
            def dummy_load(**kwargs):
                if f'{var}_depth' in kwargs.keys():
                    return (kwargs[f'{var}_val'], kwargs[f'{var}_lat'], kwargs[f'{var}_lon'])
                else:
                    return (kwargs[f'{var}_val'], kwargs[f'{var}_lat'], kwargs[f'{var}_lon'])
            def dummy_reshape(load_fcn, **kwargs):
                t = load_fcn(**kwargs)
                if len(t) == 3:
                    return dict(values=t[0], lats=t[1], lons=t[2])
                elif len(t) == 4:
                    return dict(values=t[0], lats=t[1], lons=t[2], depths=t[3])
            load_arg = dummy_load
            reshape_fcn = dummy_reshape

        elif isinstance(load_arg, (int, float)):
            kwargs[f'{var}_val'] = load_arg
            def dummy_load(**kwargs):
                return kwargs[f'{var}_val']
            def dummy_reshape(load_fcn, **kwargs):
                return dict(values=load_fcn(**kwargs))
            load_arg = dummy_load
            reshape_fcn = dummy_reshape

        else: raise TypeError(f'invalid type for load_{var}. '
            'valid types include string, array, and callable')

        # compute interpolations in parallel processes
        # child processes will serialize the result for parent to deserialize
        # if cache=False, the serialized binary will be removed from database
        interp = Process(target=serialize_interp, 
                        args=(interp_fcn, reshape_fcn, load_arg, kwargs, f'interp_{var}'))

        interp.start()
        interp.join()

        self.interp[var] = deserialize(kwargs, self.cache, f'interp_{var}')
        self.interp[var].origin = self.origin


    def get_var(self, var, grid=False, **kwargs):
        """ Get data associated with a given variable.

            Coordinates may be specified via the keyword arguments x,y,z 
            (planar geometry) or lat,lon,z (spherical geometry).

            If coordinates are specified, the methods returns the 
            interpolated data values at the specified coordinates.

            If no coordinates are specified, the method returns the 
            full underlying data array. 

            Args:
                var: str
                    Variable name
                x,y,z: array-like
                    Coordinate in planar geometry with values given in meters.
                lat,lon: array-like
                    Coordinate in spherical geometry with values given in degrees.
                grid: bool
                    Specify how to combine elements of coordinate arrays.

            Returns:
                : array-like
                Interpolated data at specified coordinates, or underlying data array
        """
        assert var in self.interp.keys(), f'Requested variable ({var}) not found.'

        is_3d = ('z' in kwargs.keys())
        is_xy = ('x' in kwargs.keys() and 'y' in kwargs.keys())
        is_ll = ('lat' in kwargs.keys() and 'lon' in kwargs.keys())

        if is_3d: z = kwargs['z']

        if is_xy:
            x, y = kwargs['x'], kwargs['y']
            if is_3d: return self.interp[var].eval_xy(x, y, z, grid)
            else:     return self.interp[var].eval_xy(x, y, grid)
        elif is_ll:
            lat, lon = kwargs['lat'], kwargs['lon']
            if is_3d: return self.interp[var].eval_ll(lat, lon, z, grid)
            else:     return self.interp[var].eval_ll(lat, lon, grid)
        else:
            return self.interp[var].get_nodes()

    def get_deriv(self, var, axis, grid=False, **kwargs):
        """ Get the derivative of a given variable.

            Coordinates must be specified via the keyword arguments x,y,z 
            (planar geometry) or lat,lon,z (spherical geometry).

            Args:
                var: str
                    Variable name
                axis: str
                    Axis along which to compute the derivative. 
                x,y,z: array-like
                    Coordinate in planar geometry with values given in meters.
                lat,lon: array-like
                    Coordinate in spherical geometry with values given in degrees.
                grid: bool
                    Specify how to combine elements of coordinate arrays.

            Returns:
                : array-like
                    Derivative
        """
        assert var in self.interp.keys(), f'Requested variable ({var}) not found.'

        is_xy = ('x' in kwargs.keys() and 'y' in kwargs.keys())
        is_ll = ('lat' in kwargs.keys() and 'lon' in kwargs.keys())

        if is_xy:
            assert axis in ('x','y'), 'if x,y are specified, axis must be either x or y'

            x, y = kwargs['x'], kwargs['y']
            return self.interp[var].eval_xy(x, y, grid,
                x_deriv_order=(axis=='x'), y_deriv_order=(axis=='y'))

        elif is_ll:
            assert axis in ('lat','lon'), 'if lat,lon are specified, axis must be either lat or lon'

            lat, lon = kwargs['lat'], kwargs['lon']
            return self.interp[var].eval_ll(lat, lon, grid,
                lat_deriv_order=(axis=='lat'), lon_deriv_order=(axis=='lon'))

        else:
            print('x,y or lat,lon must be specified')
            exit(1)


    def bathy(self, grid=False, **kwargs):
        return self.get_var('bathy', grid, **kwargs)

    def bathy_deriv(self, axis, grid=False, **kwargs):
        return self.get_deriv('bathy', axis, grid, **kwargs)

    def temp(self, grid=False, **kwargs):
        return self.get_var('temp', grid, **kwargs)

    def salinity(self, grid=False, **kwargs):
        return self.get_var('salinity', grid, **kwargs)

    def wavedir(self, grid=False, **kwargs):
        return self.get_var('wavedir', grid, **kwargs)

    def waveheight(self, grid=False, **kwargs):
        return self.get_var('waveheight', grid, **kwargs)

    def waveperiod(self, grid=False, **kwargs):
        return self.get_var('waveperiod', grid, **kwargs)

    def windspeed(self, grid=False, **kwargs):
        return self.get_var('windspeed', grid, **kwargs)
