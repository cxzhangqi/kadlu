""" enables automatic fetching of data """

from datetime import datetime, timedelta
from multiprocessing import Queue, Lock, Process

import numpy as np

from kadlu.geospatial.data_sources import source_map
from kadlu.geospatial.data_sources.data_util import serialized


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

    assert f'{var}_{source}' in source_map.fetch_map.keys() \
            or f'{var}U_{source}' in source_map.fetch_map.keys(), 'invalid query, '\
        f'could not find source for variable. options are: '\
        f'{list(f.split("_")[::-1] for f in source_map.fetch_map.keys())}'

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
                job.put((source_map.fetch_map[f'{var}U_{source}'], qry.copy()))
                job.put((source_map.fetch_map[f'{var}V_{source}'], qry.copy()))
            else: job.put((source_map.fetch_map[f'{var}_{source}'], qry.copy()))
            num += 1
        cur += step

    pxs = [Process(target=fetch_process, args=(job,key)) 
            for n in range(min(num, parallel))]
    #print(f'FETCH_HANDLER DEBUG MSG: beginning downloads in {len(pxs)} processes')
    for p in pxs: p.start()
    for p in pxs: p.join()
    job.close()

    return 

