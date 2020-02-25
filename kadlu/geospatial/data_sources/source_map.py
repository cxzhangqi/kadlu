import numpy as np
from datetime import datetime, timedelta
from multiprocessing import Queue, Lock, Process
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
        top=0, bottom=-5000,
        start=datetime(2015, 3, 1), end=datetime(2015, 3, 1, 3)
        #water_density=1, seafloor_density=1
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
            print('FETCH_PROCESS DEBUG MSG: fetch function returned false, '
                 f'skipping fetch request\ndebug: {req[1]}')
    return


def fetch_handler(var, source, step=timedelta(days=1), parallel=8, **kwargs):
    """ check fetch query hash history and generate fetch requests

        requests are batched into 24h segments and paralellized.
        coordinates are rounded to nearest outer-boundary degree integer,
        a query hash is stored if a fetch request is successful

        kwargs=dict(
            start=datetime(2013, 3, 1), end=datetime(2013, 3, 31),
            south=45, west=-65.5, north=50.5, east=-56.5,
            top=0, bottom=100
        )

        var='temp'
        source='hycom'
        parallel=8

    """

    assert f'{var}_{source}' in fetch_map.keys(), 'invalid query, '\
           f'could not find source for variable. options are: {list(f.split("_") for f in fetch_map.keys())}'

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
        else: cur += step
        if serialized(qry, f'fetch_{source}_{var}') is not False:
            print(f'FETCH_HANDLER DEBUG MSG: already fetched '
                  f'{source}_{var} {qry["start"].date().isoformat()}! '
                  f'continuing...')
            continue
        job.put((fetch_map[f'{var}_{source}'], qry.copy()))
        num += 1

    pxs = [Process(target=fetch_process, args=(job,key)) 
           for n in range(min(num, parallel))]
    print(f'FETCH_HANDLER DEBUG MSG: beginning downloads in {len(pxs)} processes')
    for p in pxs: p.start()
    for p in pxs: p.join()
    job.close()

    return 

