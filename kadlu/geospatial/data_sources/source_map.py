"""
    Utilities for data fetching specific to the 
    data sources used by Kadlu.
"""

from datetime import datetime, timedelta
from multiprocessing import Queue, Lock, Process

import numpy as np

from kadlu.geospatial.data_sources.chs      import Chs
from kadlu.geospatial.data_sources.hycom    import Hycom
from kadlu.geospatial.data_sources.era5     import Era5
from kadlu.geospatial.data_sources.wwiii    import Wwiii
from kadlu.geospatial.data_sources.data_util import                 \
        database_cfg,                                               \
        serialized


# dicts for mapping strings to callback functions
# helpful for doing stuff like passing a string 'chs' to the ocean module,
# and having the module determine which function to use for loading
fetch_map = dict(
        bathy_chs           = Chs().fetch_bathymetry,
        temp_hycom          = Hycom().fetch_temp,
        salinity_hycom      = Hycom().fetch_salinity,
        water_u_hycom       = Hycom().fetch_water_u,
        water_v_hycom       = Hycom().fetch_water_v,
        wavedir_era5        = Era5().fetch_wavedirection,
        waveheight_era5     = Era5().fetch_windwaveswellheight,
        waveperiod_era5     = Era5().fetch_waveperiod,
        windspeedU_era5     = Era5().fetch_wind_u,
        windspeedV_era5     = Era5().fetch_wind_v,
        wavedir_wwiii       = Wwiii().fetch_wavedirection,
        waveheight_wwiii    = Wwiii().fetch_windwaveheight,
        waveperiod_wwiii    = Wwiii().fetch_waveperiod,
        windspeedU_wwiii    = Wwiii().fetch_wind_u,
        windspeedV_wwiii    = Wwiii().fetch_wind_v
        )
load_map = dict(
        bathy_chs           = Chs().load_bathymetry,
        temp_hycom          = Hycom().load_temp,
        salinity_hycom      = Hycom().load_salinity,
        water_u_hycom       = Hycom().load_water_u,
        water_v_hycom       = Hycom().load_water_v,
        wavedir_era5        = Era5().load_wavedirection,
        wavedir_wwiii       = Wwiii().load_wavedirection,
        waveheight_era5     = Era5().load_windwaveswellheight,
        waveheight_wwiii    = Wwiii().load_windwaveheight,
        waveperiod_era5     = Era5().load_waveperiod,
        waveperiod_wwiii    = Wwiii().load_waveperiod,
        windspeed_era5      = Era5().load_wind_uv,
        windspeed_wwiii     = Wwiii().load_wind_uv
        )

# some reasonable default kwargs
default_val = dict(
        south=45.0, west=-68.5,
        north=51.5, east=-56.5,
        top=0, bottom=-5000,
        start=datetime(2015, 3, 1), end=datetime(2015, 3, 1, 3)
    )


def fetch_process(job, key):
    """ complete fetch requests in parallel for fetch_handler 
        job:
            job queue containing (callable, kwargs)
        key:
            multiprocessing lock (database access key)
    """
    while not job.empty():
        req = job.get()
        if not req[0](lock=key, **req[1]):
            #print('FETCH_PROCESS DEBUG MSG: fetch function returned false, '
            #        f'skipping fetch request\ndebug: {req[1]}')
            pass
    return


def fetch_handler(var, source, step=timedelta(days=1), parallel=3, **kwargs):
    """ check fetch query hash history and generate fetch requests

        requests are batched into 24h segments and paralellized.
        coordinates are rounded to nearest outer-boundary degree integer,
        a query hash is stored if a fetch request is successful

        args:
            var:
                variable type (string)
            source:
                data source (string)
            step:
                timestep size for batching fetch requests. by default, make 1 
                request per day of data
            parallel:
                number of processes to run fetch jobs

        example arguments:
            var='temp'
            source='hycom'
            step=timedelta(days=1)
            parallel=4
            kwargs=dict(
                start=datetime(2013, 3, 1), end=datetime(2013, 3, 31),
                south=45, west=-65.5, north=50.5, east=-56.5,
                top=0, bottom=100
            )

        return: nothing
    """

    assert f'{var}_{source}' in fetch_map.keys() \
            or f'{var}U_{source}' in fetch_map.keys(), 'invalid query, '\
        f'could not find source for variable. options are: '\
        f'{list(f.split("_")[::-1] for f in fetch_map.keys())}'

    if 'time' in kwargs.keys() and not 'start' in kwargs.keys():
        kwargs['start'] = kwargs['time']
        del kwargs['time']
    if not 'end' in kwargs.keys(): 
        kwargs['end'] = kwargs['start'] + timedelta(hours=3)

    np.array(list(x for x in range(100)))
    np.array(np.append([1], [x]) for x in range(10))

    key = Lock()
    job = Queue()

    # break request into gridded 24h chunks for querying
    num = 0
    qry = kwargs.copy()
    qry['south'], qry['west'] = np.floor([kwargs['south'], kwargs['west']])
    qry['north'], qry['east'] = np.ceil ([kwargs['north'], kwargs['east']])
    cur = datetime(qry['start'].year, qry['start'].month, qry['start'].day)

    # add chunks to job queue and assign processes
    while cur < kwargs['end']:
        qry['start'] = cur
        qry['end'] = cur + step 
        if var == 'bathy':  # no parallelization for non-temporal data 
            cur = kwargs['end']
            for k in ('start', 'end', 'top', 'bottom', 'lock'):
                if k in qry.keys(): del qry[k]  # trim hash indexing entropy
        if serialized(qry, f'fetch_{source}_{var}') is not False:
            #print(f'FETCH_HANDLER DEBUG MSG: already fetched '
            #      f'{source}_{var} {cur.date().isoformat()}! continuing...')
            pass
        else:
            if var == 'windspeed':
                job.put((fetch_map[f'{var}U_{source}'], qry.copy()))
                job.put((fetch_map[f'{var}V_{source}'], qry.copy()))
            else: job.put((fetch_map[f'{var}_{source}'], qry.copy()))
            num += 1
        cur += step

    pxs = [Process(target=fetch_process, args=(job,key)) 
            for n in range(min(num, parallel))]
    #print(f'FETCH_HANDLER DEBUG MSG: beginning downloads in {len(pxs)} processes')
    for p in pxs: p.start()
    for p in pxs: p.join()
    job.close()

    return 


source_map = (
    """
    CHS   (Canadian Hydrography Service)
          load_bathymetry:          bathymetric data in Canada's waterways. variable resolution \n
    ERA5  (Global environmental dataset from Copernicus Climate Data Store)
          load_windwaveswellheight: combined height of wind, waves, and swell. metres
          load_wavedirection:       mean wave direction, degrees
          load_waveperiod:          mean wave period, seconds
          load_wind_uv:             wind speed computed as sqrt(u^2 + v^2) / 2, where u, v are direction vectors
          load_wind_u:              wind speed coordinate U-vector, m/s
          load_wind_v:              wind speed coordinate V-vector, m/s \n
    HYCOM (Hybrid Coordinate Ocean Model)
          load_salinity:            g/kg salt in water
          load_temp:                degrees celsius
          load_water_uv:            ocean current computed as sqrt(u^2 + v^2) / 2, where u, v are direction vectors
          load_water_u:             ocean current coordinate U-vector, m/s
          load_water_v:             ocean current coordinate V-vector, m/s \n
    WWIII (WaveWatch Ocean Model Gen 3)
          load_wavedirection:       primary wave direction, degrees
          load_waveperiod:          primary mean wave period, seconds
          load_windwaveheight:      combined height of wind and waves, metres
          load_wind_uv:             wind speed computed as sqrt(u^2 + v^2) / 2, where u, v are direction vectors
          load_wind_u:              wind speed coordinate U-vector, m/s
          load_wind_v:              wind speed coordinate V-vector, m/s
    """)

