import numpy as np
import pickle
from multiprocessing import Process
from kadlu.geospatial.interpolation             import      \
        Interpolator2D,                                     \
        Interpolator3D
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


def serialize_interp(interpfcn, reshapefcn, loadfcn, kwargs, seed=''):
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
    """ bootstrap array data into a callable for serialization when loading """
    v = kwargs['var'] + '_'
    data = np.array((kwargs[f'{v}val'], kwargs[f'{v}lat'], kwargs[f'{v}lon']))
    if f'{v}depth' not in kwargs.keys(): return data
    return np.append(data, kwargs['f{v}depth'])


class Ocean():
    """ class for handling ocean data requests 
        
        initialization args are the keyword argmuments to be passed to
        the load functions. the loaded data will then be used for
        interpolation

        args:
            north, south:
                latitude boundaries (float)
            east, west:
                longitude boundaries (float)
            top, bottom:
                depth range. only applies to salinity and temperature (float)
            start, end:
                time range for data load query (datetime)
                note that if these kwargs are used, data will be
                averaged across the time dimension. to avoid this
                behaviour, use the 'time' arg instead
            time:
                specify a single datetime as an alternative to using 
                the start, end kwargs. the nearest fetched time data 
                will be loaded 
    """
    def __init__(self, 
            load_bathymetry = 'chs',
            load_temp       = 'hycom',
            load_salinity   = 'hycom',
            load_wavedir    = 'era5',
            load_waveheight = 'era5',
            load_waveperiod = 'era5',
            load_windspeed  = 'era5',
            cache           = True,
            **kwargs):

        vartypes = ['bathy', 'temp', 'salinity', 
                'wavedir', 'waveheight', 'waveperiod', 'windspeed']
        load_args = [load_bathymetry, load_temp, load_salinity, 
                load_wavedir, load_waveheight, load_waveperiod, load_windspeed]

        # if load_args are not callable, convert string or array to callable
        for var, load_arg, ix in zip(vartypes, load_args, range(0, len(vartypes))):
            if callable(load_arg): continue

            elif isinstance(load_arg, str):
                load_args[ix] = load_map[f'{var}_{load_arg.lower()}']

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

            else: raise TypeError(f'unknown type for load_{var}. '
                  'valid types include string, array, and callable')

        # compute interpolations in parallel processes
        # child processes will serialize the result for parent to deserialize
        # if cache=False, the serialized binary will be removed from database
        interpolations = [
                Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_args[0], kwargs, 'interp_bathy')),
                Process(target=serialize_interp, 
                        args=(Interpolator3D, reshape_3D, load_args[1], kwargs, 'interp_temp')),
                Process(target=serialize_interp, 
                        args=(Interpolator3D, reshape_3D, load_args[2], kwargs, 'interp_salinity')),
                Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_args[3], kwargs, 'interp_wavedir')), 
                Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_args[4], kwargs, 'interp_waveheight')),
                Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_args[5], kwargs, 'interp_waveperiod')),
                Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_args[6], kwargs, 'interp_wind'))
            ]

        for i in interpolations: i.start()
        for i in interpolations: i.join()

        self.interp_bathy      = deserialize(kwargs, cache, 'interp_bathy')
        self.interp_temp       = deserialize(kwargs, cache, 'interp_temp')
        self.interp_salinity   = deserialize(kwargs, cache, 'interp_salinity')
        self.interp_wavedir    = deserialize(kwargs, cache, 'interp_wavedir')
        self.interp_waveheight = deserialize(kwargs, cache, 'interp_waveheight')
        self.interp_waveperiod = deserialize(kwargs, cache, 'interp_waveperiod')
        self.interp_wind       = deserialize(kwargs, cache, 'interp_wind')


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

