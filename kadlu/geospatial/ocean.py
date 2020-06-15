""" The ocean module provides an interface to fetching, loading 
    and interpolating ocean variables.
"""
import os
import logging
from datetime import timedelta
from multiprocessing import Process, Queue

import numpy as np

from kadlu.geospatial.interpolation             import      \
        Interpolator2D,                                     \
        Interpolator3D,                                     \
        Uniform2D,                                          \
        Uniform3D
from kadlu.geospatial.data_sources.data_util    import      \
        reshape_2D,                                         \
        reshape_3D,                                         \
        dt_2_epoch,                                         \
        fmt_coords
from kadlu.geospatial.data_sources.source_map   import      \
        default_val,                                        \
        load_map,                                           \
        var3d
from kadlu.geospatial.data_sources.chs          import Chs
from kadlu.geospatial.data_sources.hycom        import Hycom
from kadlu.geospatial.data_sources.era5         import Era5
from kadlu.geospatial.data_sources.wwiii        import Wwiii
from kadlu.geospatial.data_sources.fetch_handler import fetch_handler
from kadlu.utils import center_point


def worker(interpfcn, reshapefcn, cols, var, q):
    """ compute interpolation in parallel worker process

        interpfcn:
            callback function for interpolation
        reshapefcn:
            callback function for reshaping row data into matrix format
            for interpolation
        cols:
            data as returned from load function
        var:
            variable type. used as key in Ocean().interps dictionary
        q:
            shared queue object to pass interpolation back to parent
    """
    obj = interpfcn(**reshapefcn(cols))
    q.put((var, obj))
    return


def load_callback(*, v, data, **kwargs):
    """ bootstrap data into callable to prepare for parallelization """
    return [data[key] for key in 
               (f'{v}_val', f'{v}_lat', f'{v}_lon', f'{v}_time', f'{v}_depth')
               if key in data.keys()]


class Ocean():
    """ class for handling ocean data requests 

        data will be loaded using the given data sources and boundaries
        from arguments. an interpolation for each variable will be computed in 
        parallel

        data will be averaged over time frames for interpolation. for finer
        temporal resolution, define smaller time boundaries

        any of the below load_args may also accept a callback function instead
        of a string or array value if you wish to write your own data loading
        function. the boundary arguments supplied here will be passed to the 
        callable, i.e. north, south, west, east, top, bottom, start, end

        callables or array arguments must be ordered by [val, lat, lon] for 2D 
        data, or [val, lat, lon, depth] for 3D data

        args:
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

        attrs:
            interps: dict
                Dictionary of data interpolators
            origin: tuple(float, float)
                Latitude and longitude coordinates of the centre point of the 
                geographic bounding box. This point serves as the origin of the 
                planar x-y coordinate system.
            boundaries: dict
                Bounding box for the ocean volume in space and time
    """

    def __init__(self,
            load_bathymetry=0,  load_temp=0,        load_salinity=0, 
            load_wavedir=0,     load_waveheight=0,  load_waveperiod=0, 
            load_wind_uv=0,     load_wind_u=0,      load_wind_v=0,
            load_water_uv=0,    load_water_u=0,     load_water_v=0,
            fetch=4, **kwargs):


        for kw in [k for k in ('south', 'west', 'north', 'east', 'top', 'bottom', 
                'start', 'end') if k not in kwargs.keys()]:
            kwargs[kw] = default_val[kw]

        data = {}
        callbacks = []
        vartypes = ['bathy',            'temp',             'salinity', 
                    'wavedir',          'waveheight',       'waveperiod', 
                    'wind_uv',          'wind_u',           'wind_v', 
                    'water_uv',         'water_u',          'water_v',]
        load_args = [load_bathymetry,   load_temp,          load_salinity, 
                     load_wavedir,      load_waveheight,    load_waveperiod, 
                     load_wind_uv,      load_wind_u,        load_wind_v, 
                     load_water_uv,     load_water_u,       load_water_v,]

        # if load_args are not callable, convert it to a callable function
        for v, load_arg, ix in zip(vartypes, load_args, range(len(vartypes))):
            if callable(load_arg): callbacks.append(load_arg)

            elif isinstance(load_arg, str):
                key = f'{v}_{load_arg.lower()}'
                assert key in load_map.keys(), f'no map for {key} in\n{load_map=}'
                callbacks.append(load_map[key])
                if fetch is not False:
                    fetch_handler(v, load_arg.lower(), parallel=fetch, **kwargs)

            elif isinstance(load_arg, (int, float)):
                data[f'{v}_val'] = load_arg
                data[f'{v}_lat'] = kwargs['south']
                data[f'{v}_lon'] = kwargs['west']
                data[f'{v}_time'] = dt_2_epoch(kwargs['start'])
                if v in var3d: data[f'{v}_depth'] = kwargs['top']
                callbacks.append(load_callback)

            elif isinstance(load_arg, (list, tuple, np.ndarray)):
                if len(load_arg) not in (3, 4):
                    raise ValueError(f'invalid array shape for load_{v}. '
                    'arrays must be ordered by [val, lat, lon] for 2D data, or '
                    '[val, lat, lon, depth] for 3D data')
                data[f'{v}_val'] = load_arg[0]
                data[f'{v}_lat'] = load_arg[1]
                data[f'{v}_lon'] = load_arg[2]
                if len(load_arg) == 4: data[f'{v}_depth'] = load_arg[3]
                callbacks.append(load_callback)

            else: raise TypeError(f'invalid type for load_{v}. '
                  'valid types include string, float, array, and callable')

        q = Queue()

        # prepare data pipeline
        pipe = zip(callbacks, vartypes)
        is_3D = [v in var3d for v in vartypes]
        is_arr = [not isinstance(arg, (int, float)) for arg in load_args]
        columns = [fcn(v=v, data=data, **kwargs) for fcn, v in pipe]
        intrpmap = [(Uniform2D, Uniform3D), (Interpolator2D, Interpolator3D)]
        reshapers = [reshape_3D if v else reshape_2D for v in is_3D]
        # map interpolations to dictionary
        self.interps = {}
        interpolators = map(lambda x, y: intrpmap[x][y], is_arr, is_3D)
        interpolations = map(
            lambda i,r,c,v,q=q: Process(target=worker, args=(i,r,c,v,q)),
            interpolators, reshapers, columns, vartypes
        )

        # assert that no empty arrays were returned by load function
        for col, var in zip(columns, vartypes):
            if isinstance(col, dict) or isinstance(col[0], (int, float)): continue
            assert len(col[0]) > 0, (
                    f'no data found for {var} in region {fmt_coords(kwargs)}. '
                    f'consider expanding the region')

        # compute interpolations in parallel and store in dict attribute
        if not os.environ.get('LOGLEVEL') == 'DEBUG':
            for i in interpolations: i.start()
            while len(self.interps.keys()) < len(vartypes):
                obj = q.get()
                self.interps[obj[0]] = obj[1]
            for i in interpolations: i.join()

        # debug mode: disable parallelization for nicer stack traces
        elif os.environ.get('LOGLEVEL') == 'DEBUG':
            logging.debug('OCEAN DEBUG MSG: parallelization disabled')
            for i,r,c,v in zip(interpolators, reshapers, columns, vartypes):
                logging.debug(f'interpolating {v}')
                logging.debug(f'{i = }\n{r = }\n{c = }\n{v = }')
                obj = i(**r(c))
                q.put((v, obj))

            while len(self.interps.keys()) < len(vartypes):
                obj = q.get()
                self.interps[obj[0]] = obj[1]
                logging.debug(f'done {obj[0]}... {len(self.interps.keys())}/{len(vartypes)}')

        q.close()

        # set ocean boundaries and interpolator origins
        self.boundaries = kwargs.copy()  
        self.origin = center_point(lat=[kwargs['south'], kwargs['north']], 
                                   lon=[kwargs['west'],  kwargs['east']])
        for v in vartypes: self.interps[v].origin = self.origin

        return

    def bathy(self, lat, lon, grid=False):
        return self.interps['bathy'].interp(lat, lon, grid)

    def bathy_xy(self, x, y, grid=False):
        return self.interps['bathy'].interp_xy(x, y, grid)

    def bathy_deriv(self, lat, lon, axis, grid=False):
        assert axis in ('lat', 'lon'), 'axis must be \'lat\' or \'lon\''
        return self.interps['bathy'].interp(lat, lon, grid,
              lat_deriv_order=(axis=='lat'), lon_deriv_order=(axis=='lon'))

    def bathy_deriv_xy(self, x, y, axis, grid=False):
        assert axis in ('x', 'y'), 'axis must be \'x\' or \'y\''
        return self.interps['bathy'].interp_xy(x, y, grid,
                x_deriv_order=(axis=='x'), y_deriv_order=(axis=='y'))

    def salinity(self,lat, lon, depth, grid=False):
        return self.interps['salinity'].interp(lat, lon, depth, grid)

    def salinity_xy(self, x, y, z, grid=False):
        return self.interps['salinity'].interp_xy(x, y, z, grid)

    def temp(self, lat, lon, depth, grid=False):
        return self.interps['temp'].interp(lat, lon, depth, grid)

    def temp_xy(self, x, y, z, grid=False):
        return self.interps['temp'].interp_xy(x, y, z, grid)

    def wavedir(self, lat, lon, grid=False):
        return self.interps['wavedir'].interp(lat, lon, grid)

    def wavedir_xy(self, x, y, grid=False):
        return self.interps['wavedir'].interp_xy(x, y, grid)

    def waveheight(self, lat, lon, grid=False):
        return self.interps['waveheight'].interp(lat, lon, grid)

    def waveheight_xy(self, x, y, grid=False):
        return self.interps['waveheight'].interp_xy(x, y, grid)

    def waveperiod(self, lat, lon, grid=False):
        return self.interps['waveperiod'].interp(lat, lon, grid)

    def waveperiod_xy(self, x, y, grid=False):
        return self.interps['waveperiod'].interp_xy(x, y, grid)

    def wind_uv(self, lat, lon, grid=False):
        return self.interps['wind_uv'].interp(lat, lon, grid)

    def wind_uv_xy(self, x, y, grid=False):
        return self.interps['wind_uv'].interp_xy(x, y, grid)

    def wind_u(self, lat, lon, grid=False):
        return self.interps['wind_u'].interp(lat, lon, grid)

    def wind_u_xy(self, x, y, grid=False):
        return self.interps['wind_u'].interp_xy(x, y, grid)

    def wind_v(self, lat, lon, grid=False):
        return self.interps['wind_v'].interp(lat, lon, grid)

    def wind_v_xy(self, x, y, grid=False):
        return self.interps['wind_v'].interp_xy(x, y, grid)

    def water_uv(self, lat, lon, grid=False):
        return self.interps['water_uv'].interp(lat, lon, grid)

    def water_uv_xy(self, x, y, grid=False):
        return self.interps['water_uv'].interp_xy(x, y, grid)

    def water_u(self, lat, lon, grid=False):
        return self.interps['water_u'].interp(lat, lon, grid)

    def water_u_xy(self, x, y, grid=False):
        return self.interps['water_u'].interp_xy(x, y, grid)

    def water_v(self, lat, lon, grid=False):
        return self.interps['water_v'].interp(lat, lon, grid)

    def water_v_xy(self, x, y, grid=False):
        return self.interps['water_v'].interp_xy(x, y, grid)

