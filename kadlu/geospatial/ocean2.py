import numpy as np
from multiprocessing import Process
from datetime import datetime
import pickle
from kadlu.geospatial.interpolation             import      \
        Interpolator2D,                                     \
        Interpolator3D,                                     \
        Uniform2D,                                          \
        Uniform3D
from kadlu.geospatial.data_sources.data_util    import      \
        index,                                              \
        flatten,                                            \
        reshape_2D,                                         \
        reshape_3D,                                         \
        hash_key,                                           \
        bin_db,                                             \
        deserialize
from kadlu.geospatial.data_sources.chs          import Chs
from kadlu.geospatial.data_sources.hycom        import Hycom
from kadlu.geospatial.data_sources.era5         import Era5
from kadlu.geospatial.data_sources.wwiii        import Wwiii

#import kadlu.geospatial.data_sources.gebco as gebco 

# dictionary for mapping strings to callback functions
load_map = dict(
    # bathymetry
        bathy_chs           = Chs.load_bathymetry,
        #bathy_gebco        = Gebco.load_bathymetry, 
    # temperature
        temp_hycom          = Hycom().load_temp,
    # salinity
        salinity_hycom      = Hycom().load_salinity,
    # wave direction
        wavedir_era5        = Era5().load_wavedirection,
        wavedir_wwiii       = Wwiii().load_wavedirection,
    # wave height
        waveheight_era5     = Era5().load_windwaveswellheight,
        waveheight_wwiii    = Wwiii().load_windwaveheight,
    # wave period
        waveperiod_era5     = Era5().load_waveperiod,
        waveperiod_wwiii    = Wwiii().load_waveperiod,
    # wind speed
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
    kwargs['vartype'] = seed[7:]
    obj = interpfcn(**reshapefcn(loadfcn, **kwargs))
    db.execute('INSERT OR IGNORE INTO bin VALUES (?, ?)', (key, pickle.dumps(obj)))
    conn.commit()
    return 


def load_callback(**kwargs):
    """ used to bootstrap array data into a callable for serialization """
    v = kwargs['vartype'] + '_'
    if 'depth' not in kwargs.keys():
        return np.array((
                kwargs[f'{v}val'], 
                kwargs[f'{v}lat'], 
                kwargs[f'{v}lon']
            ))
    return np.array((
            kwargs[f'{v}val'], 
            kwargs[f'{v}lat'], 
            kwargs[f'{v}lon'], 
            kwargs[f'{v}depth']
        ))


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
            cache_results=True,
            load_bathymetry = 'chs',
            load_temp       = 'hycom',
            load_salinity   = 'hycom',
            load_wavedir    = 'era5',
            load_waveheight = 'era5',
            load_waveperiod = 'era5',
            load_windspeed  = 'era5',
            **kwargs):

        for vartype, loadvar in zip(
                ['bathy', 'temp', 'salinity', 'wavedir', 
                    'waveheight', 'waveperiod', 'windspeed'], 
                [load_bathymetry, load_temp, load_salinity, load_wavedir, 
                    load_waveheight, load_waveperiod, load_windspeed]):
            if callable(loadvar): continue
            elif isinstance(loadvar, str):
                loadvar = load_map[f'{vartype}_{loadvar.lower()}']
            elif isinstance(loadvar, (list, np.ndarray)):
                kwargs[f'{vartype}_val'] = loadvar[0]
                kwargs[f'{vartype}_lat'] = loadvar[1]
                kwargs[f'{vartype}_lon'] = loadvar[2]
                if len(loadvar) == 4: kwargs[f'{vartype}_depth'] = loadvar[3]
                loadvar = load_callback
                if len(loadvar) not in (3, 4): 
                    raise ValueError(f'invalid array shape for load_{vartype}. '
                    'arrays must be ordered by [val, lat, lon] for 2D data, or '
                    '[val, lat, lon, depth] for 3D data')
            else:   raise ValueError(f'unknown type for load_{vartype}. '
                    'valid types include string, array, and callable')

        if not cache_results:
            kwargs['vartype']      = 'bathy'
            self.interp_bathy      = Interpolator2D(**reshape_2D(load_bathymetry, **kwargs))
            kwargs['vartype']      = 'temp'
            self.interp_temp       = Interpolator3D(**reshape_3D(load_temp,       **kwargs))
            kwargs['vartype']      = 'salinity'
            self.interp_salinity   = Interpolator3D(**reshape_3D(load_salinity,   **kwargs))
            kwargs['vartype']      = 'wavedir'
            self.interp_wavedir    = Interpolator2D(**reshape_2D(load_wavedir,    **kwargs))
            kwargs['vartype']      = 'waveheight'
            self.interp_waveheight = Interpolator2D(**reshape_2D(load_waveheight, **kwargs))
            kwargs['vartype']      = 'waveperiod'
            self.interp_waveperiod = Interpolator2D(**reshape_2D(load_waveperiod, **kwargs))
            kwargs['vartype']      = 'wind'
            self.interp_wind       = Interpolator2D(**reshape_2D(load_windspeed,  **kwargs))
        else:
            # compute interpolations in parallel processes
            # the resulting binary will be serialized to database for caching
            # if it already exists in the database, interpolation will be skipped
            processes = [
                    Process(
                        target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_bathymetry, kwargs, 'interp_bathy')
                    ),
                    Process(
                        target=serialize_interp, 
                        args=(Interpolator3D, reshape_3D, load_temp,       kwargs, 'interp_temp')
                    ),
                    Process(
                        target=serialize_interp, 
                        args=(Interpolator3D, reshape_3D, load_salinity,   kwargs, 'interp_salinity')
                    ),
                    Process(
                        target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_wavedir,    kwargs, 'interp_wavedir')
                    ),
                    Process(
                        target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_waveheight, kwargs, 'interp_waveheight')
                    ),
                    Process(
                        target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_waveperiod, kwargs, 'interp_waveperiod')
                    ),
                    Process(
                        target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, load_windspeed,  kwargs, 'interp_wind')
                    )
                ]
            for p in processes:
                p.start()
            for p in processes:
                p.join()

            self.interp_bathy      = deserialize(kwargs, 'interp_bathy')
            self.interp_temp       = deserialize(kwargs, 'interp_temp')
            self.interp_salinity   = deserialize(kwargs, 'interp_salinity')
            self.interp_wavedir    = deserialize(kwargs, 'interp_wavedir')
            self.interp_waveheight = deserialize(kwargs, 'interp_waveheight')
            self.interp_waveperiod = deserialize(kwargs, 'interp_waveperiod')
            self.interp_wind       = deserialize(kwargs, 'interp_wind')


    def bathy(self, lat, lon, grid=False):
        return self.interp_bathy.eval_ll(lat=lat, lon=lon, grid=grid)

    def bathy_gradient(self, lat, lon, axis='x', grid=False):
        assert axis in ('x', 'y'), 'axis must be \'x\' or \'y\''
        return self.interp_bathy.eval_ll(
                lat=lat, lon=lon, grid=grid,
                lat_deriv_order=(axis != 'x'), lon_deriv_order=(axis == 'x')
            )

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

