import numpy as np
import pickle
from multiprocessing import Process, Queue
import warnings
import time
from kadlu.geospatial.interpolation             import      \
        Interpolator2D,                                     \
        Interpolator3D
from kadlu.geospatial.data_sources.data_util    import      \
        hash_key,                                           \
        bin_db,                                             \
        reshape_2D,                                         \
        reshape_3D,                                         \
        deserialize
from kadlu.geospatial.data_sources.chs          import Chs
#from kadlu.geospatial.data_sources.gebco        import Gebco
from kadlu.geospatial.data_sources.hycom        import Hycom
from kadlu.geospatial.data_sources.era5         import Era5
from kadlu.geospatial.data_sources.wwiii        import Wwiii


# binary database for caching interpolation results
conn, db = bin_db()


# dicts for mapping strings to callback functions
fetch_map = dict(
        bathy_chs           = Chs().fetch_bathymetry,
        #bathy_gebco         = Gebco().fetch_bathymetry, 
        temp_hycom          = Hycom().fetch_temp,
        salinity_hycom      = Hycom().fetch_salinity,
        wavedir_era5        = Era5().fetch_wavedirection,
        waveheight_era5     = Era5().fetch_windwaveswellheight,
        waveperiod_era5     = Era5().fetch_waveperiod,
        windspeed_era5      = Era5().fetch_wind,
        wavedir_wwiii       = Wwiii().fetch_wavedirection,
        waveheight_wwiii    = Wwiii().fetch_windwaveheight,
        waveperiod_wwiii    = Wwiii().fetch_waveperiod,
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


def worker(interpfcn, reshapefcn, cols, q, seed):
    """ compute interpolation in parallel worker process
    
        interpfcn:
            callback function for interpolation
        reshapefcn:
            callback function for reshaping row data into matrix format
            for interpolation
        cols:
            columns of data as returned from load function
        q:
            shared queue object to pass binary back to parent
        seed:
            seed the hash to differentiate interpolation variables with
            the same set of kwargs
    """
    obj = interpfcn(**reshapefcn(cols))
    q.put((seed, obj))
    return


def load_callback(var, **kwargs):
    """ bootstrap array data into callable when loading """
    v = var + '_'
    data = np.array((kwargs[f'{v}val'], kwargs[f'{v}lat'], kwargs[f'{v}lon']))
    if f'{v}depth' not in kwargs.keys(): return data
    return np.append(data, kwargs['f{v}depth'])


def test_2d_reshape(cols):
    # testing doing interpolation by passing 1d arrays rather than reshaping it
    # for demonstration purposes only, we can make this cleaner later
    return dict(
            values=cols[0],
            lats=cols[1],
            lons=cols[2]
        )


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
        data, or [val, lat, lon, depth] for 3D data

        args:
            cache:
                boolean. if True, resulting interpolations will be stored as 
                binary to be reused later. caching is True by default
            load_bathymetry: 
                source of bathymetry data. can be 'chs' to load previously 
                fetched data, or array ordered by [val, lat, lon]
            load_temp:
                source of temperature data. can be 'hycom' to load previously
                fetched data, or array ordered by [val, lat, lon, depth]
            load_salinity:
                source of salinity data. can be 'hycom' to load previously
                fetched data, or array ordered by [val, lat, lon, depth]
            load_wavedir:
                source of wave direction data. can be 'era5' or 'wwiii' to load
                previously fetched data, or array ordered by [val, lat, lon]
            load_waveheight:
                source of wave height data. can be 'era5' or 'wwiii' to load
                previously fetched data, or array ordered by [val, lat, lon]
            load_waveperiod:
                source of wave period data. can be 'era5' or 'wwiii' to load
                previously fetched data, or array ordered by [val, lat, lon]
            load_wind:
                source of wind speed data. can be 'era5' or 'wwiii' to load
                previously fetched data, or array ordered by [val, lat, lon]
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
            load_bathymetry = 'chs',
            load_temp       = 'hycom',
            load_salinity   = 'hycom',
            load_wavedir    = 'era5',
            load_waveheight = 'era5',
            load_waveperiod = 'era5',
            load_windspeed  = 'era5',
            **kwargs):
        """
            cache           = True
            fetch           = False
            load_bathymetry = 'chs'
            load_temp       = 'hycom'
            load_salinity   = 'hycom'
            load_wavedir    = 'era5'
            load_waveheight = 'era5'
            load_waveperiod = 'era5'
            load_windspeed  = 'era5'
        """

        print('NOTICE: data will be averaged over time frames for interpolation'
              '. for finer temporal resolution, define smaller time bounds')

        vartypes = ['bathy', 'temp', 'salinity', 
                'wavedir', 'waveheight', 'waveperiod', 'windspeed']
        load_args = [load_bathymetry, load_temp, load_salinity, 
                load_wavedir, load_waveheight, load_waveperiod, load_windspeed]

        # if load_args are not callable, convert string or array to callable
        for var, load_arg, ix in zip(vartypes, load_args, range(len(vartypes))):
            if callable(load_arg): continue

            elif isinstance(load_arg, str):
                key = f'{var}_{load_arg.lower()}'
                if fetch == True: fetch_map[key](**kwargs)
                load_args[ix] = load_map[key]

            elif isinstance(load_arg, (int, float)):
                kwargs[f'{var}_val'] = np.array([load_arg for n in range(10)])
                kwargs[f'{var}_lat'] = np.linspace(kwargs['south'], kwargs['north'], 10)
                kwargs[f'{var}_lon'] = np.linspace(kwargs['west'], kwargs['east'], 10)
                if var in ('temp', 'salinity'):
                    kwargs[f'{var}_depth'] = np.linspace(kwargs['top'], kwargs['bottom'], 25)
                load_args[ix] = load_callback

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

            else: raise TypeError(f'invalid type for load_{var}. '
                  'valid types include string, array, and callable')


        key = hash_key(kwargs, 'interp_ocean')
        db.execute('SELECT * FROM bin WHERE hash == ? LIMIT 1', (key,))
        res = db.fetchone()
        if res is not None: self.interp = pickle.loads(res[1])
        else:
            # compute interpolations in parallel processes
            print('preparing interpolations...')
            q = Queue()
            cols = [callback(var=var, **kwargs) for var, callback in zip(vartypes, load_args)]

            interpolations = [
                    Process(target=worker,
                            args=(Interpolator2D, test_2d_reshape, cols[0], q, 'interp_bathy')),
                    Process(target=worker,
                            args=(Interpolator3D, reshape_3D, cols[1], q, 'interp_temp')),
                    Process(target=worker,
                            args=(Interpolator3D, reshape_3D, cols[2], q, 'interp_salinity')),
                    Process(target=worker,
                            args=(Interpolator2D, reshape_2D, cols[3], q, 'interp_wavedir')), 
                    Process(target=worker,
                            args=(Interpolator2D, reshape_2D, cols[4], q, 'interp_waveheight')),
                    Process(target=worker,
                            args=(Interpolator2D, reshape_2D, cols[5], q, 'interp_waveperiod')),
                    Process(target=worker,
                            args=(Interpolator2D, reshape_2D, cols[6], q, 'interp_wind'))
                ]

            for i in interpolations: i.start()
            
            self.interp = {}
            joined = 0;
            #for i in range(len(interpolations)):
            while joined < 7:
                binary = q.get()
                self.interp[binary[0].split('_')[1]] = binary[1]
                joined += 1
            for i in interpolations: i.join()
            assert q.empty()  # all done

            if cache: 
                db.execute('INSERT INTO bin VALUES (?, ?)', (key, pickle.dumps(self.interp)))
                conn.commit()


    def bathy(self, lat, lon, grid=False):
        return self.interp['bathy'].eval_ll(lat=lat, lon=lon, grid=grid)

    def bathy_xy():
        pass

    def bathy_gradient(self, lat, lon, axis='x', grid=False):
        assert axis in ('x', 'y'), 'axis must be \'x\' or \'y\''
        return self.interp['bathy'].eval_ll(lat=lat, lon=lon, grid=grid,
                lat_deriv_order=(axis != 'x'), lon_deriv_order=(axis == 'x'))

    def bathy_gradient_xy():
        pass

    def temp(self, lat, lon, depth, grid=False):
        return self.interp['temp'].eval_ll(lat=lat, lon=lon, z=depth, grid=grid)

    def temp_xy():
        pass

    def salinity(self, lat, lon, depth, grid=False):
        return self.interp['salinity'].eval_ll(lat=lat, lon=lon, z=depth, grid=grid)

    def salinity_xy():
        pass

    def wavedir(self, lat, lon, grid=False):
        return self.interp['wavedir'].eval_ll(lat=lat, lon=lon, grid=grid)

    def wavedir_xy():
        pass

    def waveheight(self, lat, lon, grid=False):
        return self.interp['waveheight'].eval_ll(lat=lat, lon=lon, grid=grid)

    def waveheight_xy():
        pass

    def waveperiod(self, lat, lon, grid=False):
        return self.interp['waveperiod'].eval_ll(lat=lat, lon=lon, grid=grid)

    def waveperiod_xy():
        pass

    def windspeed(self, lat, lon, grid=False):
        return self.interp['wind'].eval_ll(lat=lat, lon=lon, grid=grid)

    def windspeed_xy():
        pass

"""
from kadlu.geospatial.data_sources.data_util import gen_kwargs
kwargs = gen_kwargs()
self = Ocean(**kwargs)
self.temp([42.5], [-62.5], [100])
"""

