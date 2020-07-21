"""
    some data fetching utils, function maps, and constant variables
"""

from datetime import datetime, timedelta

import kadlu.geospatial.data_sources.chs as chs
import kadlu.geospatial.data_sources.era5 as era5
import kadlu.geospatial.data_sources.gebco as gebco
import kadlu.geospatial.data_sources.hycom as hycom
import kadlu.geospatial.data_sources.wwiii as wwiii


# dicts for mapping strings to callback functions
# helpful for passing source strings to the ocean module,
# and having the module determine which function to use for loading
fetch_map = dict(
        bathy_chs           = chs  .Chs()  .fetch_bathymetry,
        temp_hycom          = hycom.Hycom().fetch_temp,
        salinity_hycom      = hycom.Hycom().fetch_salinity,
        water_uv_hycom      = hycom.Hycom().fetch_water_v,
        water_u_hycom       = hycom.Hycom().fetch_water_u,
        water_v_hycom       = hycom.Hycom().fetch_water_v,
        wavedir_era5        = era5 .Era5() .fetch_wavedirection,
        waveheight_era5     = era5 .Era5() .fetch_windwaveswellheight,
        waveperiod_era5     = era5 .Era5() .fetch_waveperiod,
        wind_uv_era5        = era5 .Era5() .fetch_wind_uv,
        wind_u_era5         = era5 .Era5() .fetch_wind_u,
        wind_v_era5         = era5 .Era5() .fetch_wind_v,
        wavedir_wwiii       = wwiii.Wwiii().fetch_wavedirection,
        waveheight_wwiii    = wwiii.Wwiii().fetch_windwaveheight,
        waveperiod_wwiii    = wwiii.Wwiii().fetch_waveperiod,
        wind_uv_wwiii       = wwiii.Wwiii().fetch_wind_uv,
        wind_u_wwiii        = wwiii.Wwiii().fetch_wind_u,
        wind_v_wwiii        = wwiii.Wwiii().fetch_wind_v,
        bathy_gebco         = gebco.Gebco().fetch_bathymetry,
    )

load_map = dict(
        bathy_chs           = chs  .Chs()  .load_bathymetry,
        temp_hycom          = hycom.Hycom().load_temp,
        salinity_hycom      = hycom.Hycom().load_salinity,
        water_uv_hycom      = hycom.Hycom().load_water_uv,
        water_u_hycom       = hycom.Hycom().load_water_u,
        water_v_hycom       = hycom.Hycom().load_water_v,
        wavedir_era5        = era5 .Era5() .load_wavedirection,
        waveheight_era5     = era5 .Era5() .load_windwaveswellheight,
        waveperiod_era5     = era5 .Era5() .load_waveperiod,
        wind_uv_era5        = era5 .Era5() .load_wind_uv,
        wind_u_era5         = era5 .Era5() .load_wind_u,
        wind_v_era5         = era5 .Era5() .load_wind_v,
        wavedir_wwiii       = wwiii.Wwiii().load_wavedirection,
        waveheight_wwiii    = wwiii.Wwiii().load_windwaveheight,
        waveperiod_wwiii    = wwiii.Wwiii().load_waveperiod,
        wind_uv_wwiii       = wwiii.Wwiii().load_wind_uv,
        wind_u_wwiii        = wwiii.Wwiii().load_wind_u,
        wind_v_wwiii        = wwiii.Wwiii().load_wind_v,
        bathy_gebco         = gebco.Gebco().load_bathymetry,
    )

# some reasonable default kwargs
default_val = dict(
        south=44.25, west=-64.5,
        north=44.70, east=-63.33,
        top=0, bottom=5000,
        start=datetime(2015, 3, 1), end=datetime(2015, 3, 1, 12)
    )

var3d = ('temp', 'salinity', 'water_u', 'water_v', 'water_uv',)

source_map = (
    """
    CHS   (Canadian Hydrography Service)
          bathymetry:       bathymetric data in Canada's waterways. metres, variable resolution \n
    ERA5  (Global environmental dataset from Copernicus Climate Data Store)
          wavedir:          mean wave direction, degrees
          waveheight:       combined height of wind, waves, and swell. metres
          waveperiod:       mean wave period, seconds
          wind_uv:          wind speed computed as sqrt(u^2 + v^2) / 2, where u, v are direction vectors
          wind_u:           wind speed coordinate U-vector, m/s
          wind_v:           wind speed coordinate V-vector, m/s \n
    GEBCO (General Bathymetric Chart of the Oceans)
          bathymetry:       global bathymetric and topographic data. metres below sea level \n
    HYCOM (Hybrid Coordinate Ocean Model)
          salinity:         g/kg salt in water
          temp:             degrees celsius
          water_uv:         ocean current computed as sqrt(u^2 + v^2) / 2, where u, v are direction vectors
          water_u:          ocean current coordinate U-vector, m/s
          water_v:          ocean current coordinate V-vector, m/s \n
    WWIII (WaveWatch Ocean Model Gen 3)
          wavedir:          primary wave direction, degrees
          waveheight:       combined height of wind and waves, metres
          waveperiod:       primary mean wave period, seconds
          wind_uv:          wind speed computed as sqrt(u^2 + v^2) / 2, where u, v are direction vectors
          wind_u:           wind speed coordinate U-vector, m/s
          wind_v:           wind speed coordinate V-vector, m/s
    """)

