import numpy as np
import multiprocessing 
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
        serialize,                                          \
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
    def __init__(self, **kwargs):
        # later we can define the preferred data sources on class initialization
        # or config.ini. for now, just use these ones
        self.load_bathymetry = Chs().load_bathymetry
        self.load_temp = Hycom().load_temp
        self.load_salinity = Hycom().load_salinity
        self.load_wavedirection = Era5().load_wavedirection
        self.load_waveheight = Era5().load_windwaveswellheight
        self.load_waveperiod = Era5().load_waveperiod
        
        #self.set_origin(0, 0)
        #self.SW = LatLon(-90, -180)
        #self.NE = LatLon(90, 180)

        # temporary placeholder values for testing
        kwargs['water_density'] = 1.0
        kwargs['origin']        = (0, 0)  # lat_ref, lon_ref
        kwargs['method']        = 'linear'
        #kwargs['geometry']      = 'spherical'

        """
        interpolate data and store result
        note that the reshape_2D function is incomplete
        the data keys will be updated accordingly when finished
        """

        #data = reshape_2D(self.load_bathymetry, **kwargs)
        #serialize(kwargs, (Interpolator2D(**data), data), 'interp_bathy')
        #self.bathy_data = reshape_2D(self.load_bathymetry, kwargs)
        #self.bathy_interp = Interpolator2D(**self.bathy_data)

        #data = reshape_3D(self.load_temp, **kwargs)
        #serialize(kwargs, (Interpolator3D(**data), data), 'interp_temperature')
        self.temp_data = reshape_3D(self.load_temp, kwargs)
        self.temp_interp = Interpolator3D(**self.temp_data)
        
        self.salinity_data = reshape_3D(self.load_salinity, kwargs)
        self.salinity_interp = Interpolator3D(**self.salinity_data)

    def bathy(self, lat, lon, grid=False):
        #interp, data = deserialize(kwargs, 'interp_bathy')
        return self.interp_bathy.eval_ll(lat=lat, lon=lon, grid=grid)

    def bathy_gradient(self, lat, lon, axis='x', grid=False):
        #interp, data = deserialize(kwargs, 'interp_bathy')
        return self.bathy_interp.eval_ll(
                lat=lat, 
                lon=lon, 
                grid=grid,
                lat_deriv_order=(axis != 'x'), 
                lon_deriv_order=(axis == 'x')
            )
        
    def temp(self, lat, lon, z, grid=False):
        #interp, data = deserialize(kwargs, 'interp_temperature')
        return self.temp_interp.eval_ll(lat=lat, lon=lon, z=z, grid=grid)

    def salinity(self, lat, lon, z, grid=False):
        #interp, data = deserialize(kwargs, 'interp_salinity')
        return self.salinity_interp.eval_ll(lat=lat, lon=lon, z=z, grid=grid)

    def wave_direction(self, lat, lon, grid=False):
        #data = reshape_2D(self.load_wavedirection, **kwargs)
        #return Interpolator2D(**data).eval_ll(lat=data['lats'], lon=data['lons'], grid=data['values'])
        pass

    def wave_height():
        #kwargs = default(kwargs)  # temporary fix
        #data = reshape_2D(self.load_waveheight, **kwargs)
        #return Interpolator2D(**data).eval_ll(lat=data['lats'], lon=data['lons'], grid=data['values'])
        pass

    def wave_period():
        #kwargs = default(kwargs)  # temporary fix
        #data = reshape_2D(self.load_waveperiod, **kwargs)
        #return Interpolator2D(**data).eval_ll(lat=data['lats'], lon=data['lons'], grid=data['values'])
        pass

    def wind_speed():
        # need to think about how to reshape windspeed tuple data...
        # take square mean before?
        pass

