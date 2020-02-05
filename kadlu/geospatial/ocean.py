import numpy as np
import pickle
from multiprocessing import Process
from kadlu.geospatial.interpolation             import      \
        Interpolator2D,                                     \
        Interpolator3D,                                     \
        Uniform2D,                                          \
        Uniform3D
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
    return np.append(data, kwargs['f{v}depth'])


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
            load_bathymetry = None,
            load_temp       = None,
            load_salinity   = None,
            load_wavedir    = None,
            load_waveheight = None,
            load_waveperiod = None,
            load_windspeed  = None,
            **kwargs):

        # default data sources
        default_sources = {'bathy':'chs', 'temp':'hycom', 'salinity':'hycom',
                'wavedir':'era5', 'waveheight':'era5', 'waveperiod':'era5', 'windspeed':'era5'}

        if 'start' in kwargs.keys() and 'end' in kwargs.keys(): 
            print('WARNING: data will be averaged over time frames for interpolation')
                #.\nto avoid this behaviour, use the \'time\' '
                #'keyword argument instead of start/end')

        vardims = [2, 3, 3, 2, 2, 2, 2]
        vartypes = ['bathy', 'temp', 'salinity', 
                'wavedir', 'waveheight', 'waveperiod', 'windspeed']
        load_args = [load_bathymetry, load_temp, load_salinity, 
                load_wavedir, load_waveheight, load_waveperiod, load_windspeed]

        interp_opt = {2: Interpolator2D, 3: Interpolator3D} #interpolation options
        reshape_opt = {2: reshape_2D, 3: reshape_3D}        #reshape options
        interp_uni_opt = {2: Uniform2D, 3: Uniform3D}       #interpolation options for uniform data

        interp_fcn = [interp_opt[dim] for dim in vardims]
        reshape_fcn = [interp_opt[dim] for dim in vardims]

        # if default=False, replace None's with 0
        # otherwise, replace None's with default data sources
        for i in range(len(load_args)):
            if load_args[i] is None:
                if default == False: load_args[i] = 0
                else: load_args[i] = default_sources[vartypes[i]]

        # if load_args are not callable, convert string or array to callable
        for var, dim, load_arg, ix in zip(vartypes, vardims, load_args, range(len(vartypes))):
            if callable(load_arg): continue

            elif isinstance(load_arg, str):
                key = f'{var}_{load_arg.lower()}'
                if fetch == True: fetch_map[key](**kwargs)
                load_args[ix] = load_map[key]

            elif isinstance(load_arg, (list, tuple, np.ndarray)):
                if len(load_arg) not in (3, 4):
                    raise ValueError(f'invalid array shape for load_{var}. '
                    'arrays must be ordered by [val, lat, lon] for 2D data, or '
                    '[val, lat, lon, depth] for 3D data')
                kwargs[f'{var}_val'] = load_arg[0]
                kwargs[f'{var}_lat'] = load_arg[1]
                kwargs[f'{var}_lon'] = load_arg[2]
                if len(load_arg) == 4: kwargs[f'{var}_depth'] = load_arg[3]
                load_args[ix] = load_callback

            elif isinstance(load_arg, (int, float)):
                kwargs[f'{var}_val'] = load_arg
                load_args[ix] = lambda x: x
                interp_fcn[ix] = interp_opt[dim]
                reshape_fcn[ix] = lambda x: x

            else: raise TypeError(f'invalid type for load_{var}. '
                  'valid types include string, array, and callable')

        # compute interpolations in parallel processes
        # child processes will serialize the result for parent to deserialize
        # if cache=False, the serialized binary will be removed from database
        interpolations = []
        for i in range(len(interp_fcn)):
            proc = Process(target=serialize_interp, 
                           args=(interp_fcn[i], reshape_fcn[i], load_args[i], kwargs, 'interp_{0}'.format(vartypes[i])))
            interpolations.append(proc)

#                Process(target=serialize_interp, 
#                        args=(Interpolator2D, reshape_2D, load_args[0], kwargs, 'interp_bathy')),
#                Process(target=serialize_interp, 
#                        args=(Interpolator3D, reshape_3D, load_args[1], kwargs, 'interp_temp')),
#                Process(target=serialize_interp, 
#                        args=(Interpolator3D, reshape_3D, load_args[2], kwargs, 'interp_salinity')),
#                Process(target=serialize_interp, 
#                        args=(Interpolator2D, reshape_2D, load_args[3], kwargs, 'interp_wavedir')), 
#                Process(target=serialize_interp, 
#                        args=(Interpolator2D, reshape_2D, load_args[4], kwargs, 'interp_waveheight')),
#                Process(target=serialize_interp, 
#                        args=(Interpolator2D, reshape_2D, load_args[5], kwargs, 'interp_waveperiod')),
#                Process(target=serialize_interp, 
#                        args=(Interpolator2D, reshape_2D, load_args[6], kwargs, 'interp_windspeed'))
#            ]
        for i in interpolations: i.start()
        for i in interpolations: i.join()

        self.interp_bathy      = deserialize(kwargs, cache, 'interp_bathy')
        self.interp_temp       = deserialize(kwargs, cache, 'interp_temp')
        self.interp_salinity   = deserialize(kwargs, cache, 'interp_salinity')
        self.interp_wavedir    = deserialize(kwargs, cache, 'interp_wavedir')
        self.interp_waveheight = deserialize(kwargs, cache, 'interp_waveheight')
        self.interp_waveperiod = deserialize(kwargs, cache, 'interp_waveperiod')
        self.interp_wind       = deserialize(kwargs, cache, 'interp_windspeed')


    def bathy(self, lat, lon, grid=False):
        return self.interp_bathy.eval_ll(lat=lat, lon=lon, grid=grid)

    def bathy_gradient(self, lat, lon, axis='x', grid=False):
        assert axis in ('x', 'y'), 'axis must be \'x\' or \'y\''
        return self.interp_bathy.eval_ll(
                lat=lat, lon=lon, grid=grid,
                lat_deriv_order=(axis != 'x'), lon_deriv_order=(axis == 'x'))

    def temp(self, lat, lon, z, grid=False):
        return self.interp_temp.eval_ll(lat=lat, lon=lon, z=z, grid=grid)

    def salinity(self, lat, lon, z, grid=False):
        return self.interp_salinity.eval_ll(lat=lat, lon=lon, z=z, grid=grid)

    def wavedir(self, lat, lon, grid=False):
        return self.interp_wavedir.eval_ll(lat=lat, lon=lon, grid=grid)

    def waveheight(self, lat, lon, grid=False):
        return self.interp_waveheight.eval_ll(lat=lat, lon=lon, grid=grid)

    def waveperiod(self, lat, lon, grid=False):
        return self.interp_waveperiod.eval_ll(lat=lat, lon=lon, grid=grid)

    def windspeed(self, lat, lon, grid=False):
        return self.interp_wind.eval_ll(lat=lat, lon=lon, grid=grid)

