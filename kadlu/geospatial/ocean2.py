import numpy as np
from multiprocessing import Process
from datetime import datetime
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
        serialize_interp,                                   \
        deserialize
from kadlu.geospatial.data_sources.chs          import Chs
from kadlu.geospatial.data_sources.hycom        import Hycom
from kadlu.geospatial.data_sources.era5         import Era5
from kadlu.geospatial.data_sources.wwiii        import Wwiii

#import kadlu.geospatial.data_sources.gebco as gebco 


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
    def __init__(self, cache_results=True, **kwargs):
        self.load_bathymetry        = Chs().load_bathymetry
        self.load_temp              = Hycom().load_temp
        self.load_salinity          = Hycom().load_salinity
        self.load_wavedirection     = Era5().load_wavedirection
        self.load_waveheight        = Era5().load_windwaveswellheight
        self.load_waveperiod        = Era5().load_waveperiod
        self.load_windspeed         = Era5().load_wind

        if not cache_results:
            self.bathy_interp         = Interpolator2D( **reshape_2D(self.load_bathymetry,    **kwargs))
            self.temp_interp          = Interpolator3D( **reshape_3D(self.load_temp,          **kwargs))
            self.salinity_interp      = Interpolator3D( **reshape_3D(self.load_salinity,      **kwargs))
            self.wavedirection_interp = Interpolator2D( **reshape_2D(self.load_wavedirection, **kwargs))
            self.waveheight_interp    = Interpolator2D( **reshape_2D(self.load_waveheight,    **kwargs))
            self.waveperiod_interp    = Interpolator2D( **reshape_2D(self.load_waveperiod,    **kwargs))
            self.windspeed_interp     = Interpolator2D( **reshape_2D(self.load_windspeed,     **kwargs))
        else:
            # run interpolations in parallel processes
            # the resulting binary will be serialized to database for caching
            # if it already exists in the database, interpolation will be skipped
            # and the cached result will be used 
            processes = [
                    Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, self.load_bathymetry, kwargs, 'interp_bathy')),
                    Process(target=serialize_interp, 
                        args=(Interpolator3D, reshape_3D, self.load_temp, kwargs, 'interp_temp')),
                    Process(target=serialize_interp, 
                        args=(Interpolator3D, reshape_3D, self.load_salinity, kwargs, 'interp_salinity')),
                    Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, self.load_wavedirection, kwargs, 'interp_wavedirection')),
                    Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, self.load_waveheight, kwargs, 'interp_waveheight')),
                    Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, self.load_waveperiod, kwargs, 'interp_waveperiod')),
                    Process(target=serialize_interp, 
                        args=(Interpolator2D, reshape_2D, self.load_windspeed, kwargs, 'interp_windspeed'))
                ]
            for p in processes:
                p.start()
            for p in processes:
                p.join()

            self.bathy_interp         = deserialize(kwargs, 'interp_bathy')
            self.temp_interp          = deserialize(kwargs, 'interp_temp')
            self.salinity_interp      = deserialize(kwargs, 'interp_salinity')
            self.wavedirection_interp = deserialize(kwargs, 'interp_wavedirection')
            self.waveheight_interp    = deserialize(kwargs, 'interp_waveheight')
            self.waveperiod_interp    = deserialize(kwargs, 'interp_waveperiod')
            self.windspeed_interp     = deserialize(kwargs, 'interp_windspeed')


    def bathy(self, lat, lon, grid=False):
        return self.bathy_interp.eval_ll(lat=lat, lon=lon, grid=grid)

    def bathy_gradient(self, lat, lon, axis='x', grid=False):
        assert axis in ('x', 'y'), 'axis must be \'x\' or \'y\''
        return self.bathy_interp.eval_ll(
                lat=lat, lon=lon, grid=grid,
                lat_deriv_order=(axis != 'x'), lon_deriv_order=(axis == 'x')
            )

    def temp(self, lat, lon, z, grid=False):
        return self.temp_interp.eval_ll(lat=lat, lon=lon, z=z, grid=grid)

    def salinity(self, lat, lon, z, grid=False):
        return self.salinity_interp.eval_ll(lat=lat, lon=lon, z=z, grid=grid)

    def wave_direction(self, lat, lon, grid=False):
        return self.wavedirection_interp.eval_ll(lat=lat, lon=lon, grid=grid)

    def wave_height(self, lat, lon, grid=False):
        return self.waveheight_interp.eval_ll(lat=lat, lon=lon, grid=grid)

    def wave_period(self, lat, lon, grid=False):
        return self.waveperiod_interp.eval_ll(lat=lat, lon=lon, grid=grid)

    def wind_speed(self, lat, lon, grid=False):
        return self.windspeed_interp.eval_ll(lat=lat, lon=lon, grid=grid)

