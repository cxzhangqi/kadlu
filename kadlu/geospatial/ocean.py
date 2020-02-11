import numpy as np
import pickle
from multiprocessing import Process, Queue
import warnings
import time
from kadlu.geospatial.interpolation             import      \
        Interpolator2D,                                     \
        Interpolator3D,                                     \
        Uniform2D,                                          \
        Uniform3D
from kadlu.geospatial.data_sources.data_util    import      \
        hash_key,                                           \
        bin_db,                                             \
        reshape_2D,                                         \
        reshape_3D,                                         \
        dt_2_epoch
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


def worker(interpfcn, reshapefcn, data, var, q):
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
    obj = interpfcn(**reshapefcn(var=var, data=data))
    q.put((var, obj))
    return


def load_callback(*, data, var, **kwargs):
    """ bootstrap data into callable when loading """
    v = var + '_'
    if f'{v}time' not in data.keys(): 
        return (data[f'{v}val'], data[f'{v}lat'], data[f'{v}lon'])
    if f'{v}depth' not in data.keys(): 
        return (data[f'{v}val'], data[f'{v}lat'], data[f'{v}lon'], data[f'{v}time'])
    return (data[f'{v}val'], data[f'{v}lat'], data[f'{v}lon'], data[f'{v}time'], data[f'{v}depth'])


class Ocean():
    """ class for handling ocean data requests 

        data will be loaded using the given data sources and boundaries
        from arguments. an interpolation for each variable will be computed in 
        parallel

        any of the below load_args may also accept a callback function instead
        of a string or array value if you wish to write your own data loading
        function. the boundary arguments supplied here will be passed to the 
        callable, i.e. north, south, west, east, top, bottom, start, end

        bles or array arguments must be ordered by [val, lat, lon] for 2D 
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
            cache           = False,
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
            cache           = False 
            fetch           = False
            load_bathymetry = 0#'chs'
            load_temp       = 0#'hycom',
            load_salinity   = 0#'hycom'
            load_wavedir    = 0#'era5'
            load_waveheight = 0#'era5'
            load_waveperiod = 0#'era5'
            load_windspeed  = 0#'era5'
        """

        print('data will be averaged over time frames for interpolation. '
              'for finer temporal resolution, define smaller time bounds')

        vartypes = ['bathy', 'temp', 'salinity', 
                'wavedir', 'waveheight', 'waveperiod', 'windspeed']
        load_args = [load_bathymetry, load_temp, load_salinity, 
                load_wavedir, load_waveheight, load_waveperiod, load_windspeed]
        callbacks = []
        data = {}

        # if load_args are not callable, convert string or array to callable
        for var, load_arg, ix in zip(vartypes, load_args, range(len(vartypes))):
            if callable(load_arg): continue

            elif isinstance(load_arg, str):
                key = f'{var}_{load_arg.lower()}'
                if fetch == True: fetch_map[key](**kwargs)
                callbacks.append(load_map[key])

            elif isinstance(load_arg, (int, float)):
                data[f'{var}_val'] = load_arg
                data[f'{var}_lat'] = kwargs['south'] 
                data[f'{var}_lon'] = kwargs['west']
                data[f'{var}_time'] = dt_2_epoch(kwargs['start'])[0] 
                if var in ('temp', 'salinity'):
                    data[f'{var}_depth'] = np.linspace(kwargs['top'], kwargs['bottom'], 10)
                callbacks.append(load_callback)

            elif isinstance(load_arg, (list, tuple, np.ndarray)):
                if len(load_arg) not in (3, 4):
                    raise ValueError(f'invalid array shape for load_{var}. '
                    'arrays must be ordered by [val, lat, lon] for 2D data, or '
                    '[val, lat, lon, depth] for 3D data')
                data[f'{var}_val'] = load_arg[0]
                data[f'{var}_lat'] = load_arg[1]
                data[f'{var}_lon'] = load_arg[2]
                if len(load_arg) == 4: data[f'{var}_depth'] = load_arg[3]
                callbacks.append(load_callback)

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
            interpolators = []
            var3D = ('temp', 'salinity')
            for v, arg in zip(vartypes, load_args):
                if isinstance(arg, (float, int)):
                    if v in var3D: interpolators.append(Uniform3D)
                    else: interpolators.append(Uniform2D)
                else:
                    if v in var3D: interpolators.append(Interpolator3D)
                    else: interpolators.append(Interpolator2D)

            reshapers = [reshape_2D if v not in ('temp', 'salinity') else reshape_3D for v in vartypes]
            columns = [callback(data=data, var=v, **kwargs) for v, callback in zip(vartypes, callbacks)]
            interpolations = map(
                lambda i,r,c,v: Process(target=worker, args=(i,r,c,v,q)),
                interpolators, reshapers, columns, vartypes
            )

            joined = 0;
            self.interp = {}
            for i in interpolations: i.start()
            while joined < len(vartypes):
                obj = q.get()
                self.interp[obj[0]] = obj[1]
                joined += 1
            for i in interpolations: i.join()
            q.close()

            if cache: 
                db.execute('INSERT INTO bin VALUES (?, ?)', (key, pickle.dumps(self.interp)))
                conn.commit()

            #self.lat_default = np.linspace(kwargs['south'], kwargs['north'],100)
            #self.lon_default = np.linspace(kwargs['west'],  kwargs['east'], 100)
            #self.depth_default = np.linspace(kwargs['top'],  kwargs['bottom'], 100)

            self.lat_default = kwargs['north']
            self.lon_default = kwargs['west']
            self.depth_default = kwargs['top']


    def bathy(self, lat=None, lon=None, grid=False):
        if lat is None: lat=self.lat_default
        if lon is None: lon=self.lon_default
        return self.interp['bathy'].eval_ll(lat=lat, lon=lon, grid=grid)

    def bathy_xy(self, x, y, grid=False):
        pass

    def bathy_deriv(self, lat=None, lon=None, axis='x', grid=False):
        if lat == None: lat=self.lat_default
        if lon == None: lon=self.lon_default
        assert axis in ('x', 'y'), 'axis must be \'x\' or \'y\''
        return self.interp['bathy'].eval_ll(lat=lat, lon=lon, grid=grid,
                lat_deriv_order=(axis != 'x'), lon_deriv_order=(axis == 'x'))

    def bathy_deriv_xy(self, x, y, grid=False):
        pass

    def temp(self, lat=None, lon=None, depth=None, grid=False):
        if lat == None: lat=self.lat_default
        if lon == None: lon=self.lon_default
        if depth == None: depth=self.depth_default
        return self.interp['temp'].eval_ll(lat=lat, lon=lon, z=depth, grid=grid)

    def temp_xy(self, x, y, z, grid=False):
        pass

    def salinity(self,lat=None, lon=None, depth=None, grid=False):
        if lat == None: lat=self.lat_default
        if lon == None: lon=self.lon_default
        if depth == None: depth=self.depth_default
        return self.interp['salinity'].eval_ll(lat=lat, lon=lon, z=depth, grid=grid)

    def salinity_xy(self, x, y, z, grid=False):
        pass

    def wavedir(self, lat=None, lon=None, depth=None, grid=False):
        if lat == None: lat=self.lat_default
        if lon == None: lon=self.lon_default
        return self.interp['wavedir'].eval_ll(lat=lat, lon=lon, grid=grid)

    def wavedir_xy(self, x, y, grid=False):
        pass

    def waveheight(self, lat=None, lon=None, depth=None, grid=False):
        if lat == None: lat=self.lat_default
        if lon == None: lon=self.lon_default
        return self.interp['waveheight'].eval_ll(lat=lat, lon=lon, grid=grid)

    def waveheight_xy(self, x, y, grid=False):
        pass

    def waveperiod(self, lat=None, lon=None, depth=None, grid=False):
        if lat == None: lat=self.lat_default
        if lon == None: lon=self.lon_default
        return self.interp['waveperiod'].eval_ll(lat=lat, lon=lon, grid=grid)

    def waveperiod_xy(self, x, y, grid=False):
        pass

    def windspeed(self, lat=None, lon=None, depth=None, grid=False):
        if lat == None: lat=self.lat_default
        if lon == None: lon=self.lon_default
        return self.interp['windspeed'].eval_ll(lat=lat, lon=lon, grid=grid)

    def windspeed_xy(self, x, y, grid=False):
        pass

"""
from kadlu.geospatial.data_sources.data_util import gen_kwargs
kwargs = gen_kwargs()
self = Ocean(**kwargs)
self.temp([42.5], [-62.5], [100])
"""

