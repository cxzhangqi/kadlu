from kadlu.geospatial.data_sources.chs import Chs
from kadlu.geospatial.data_sources.hycom import Hycom
from kadlu.geospatial.data_sources.era5 import Era5
from kadlu.geospatial.data_sources.wwiii import Wwiii
from datetime import datetime, timedelta

# dicts for mapping strings to callback functions
# helpful for doing stuff like passing a string 'chs' to the ocean module,
# and having the module determine which function to use for loading

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

# some reasonable default kwargs
default_val = dict(
        south=45.0, west=-68.5,
        north=51.5, east=-56.5,
        top=0, bottom=5000,
        start=datetime.now()-timedelta(hours=3), end=datetime.now(),
        water_density=1, seafloor_density=1
    )

