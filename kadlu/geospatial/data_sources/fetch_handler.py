""" enables automatic fetching of data """

import time
import logging
from os import getpid
from datetime import datetime, timedelta

import numpy as np

from kadlu.geospatial.data_sources import source_map
from kadlu.geospatial.data_sources.data_util import serialized


def fetch_handler(var, src, dx=2, dy=2, dt=timedelta(days=1), **kwargs):
    """ check fetch query hash history and generate fetch requests

        requests are batched into dx° * dy° * dt request bins,
        with the entire range of depths included in each bin.
        coordinates are rounded to nearest outer-boundary degree integer,
        a query hash is stored if a fetch request is successful

        args:
            var:
                variable type (string)
                must be one of the variables listed in source_map
            src:
                data source (string)
                must be one of the sources listed in source_map
            dx:
                delta longitude bin size (int)
            dy: 
                delta latitude bin size (int)
            dt:
                delta time bin size (timedelta)

        return: nothing
    """

    assert f'{var}_{src}' in source_map.fetch_map.keys() \
            or f'{var}U_{src}' in source_map.fetch_map.keys(), 'invalid query, '\
        f'could not find source for variable. options are: '\
        f'{list(f.split("_")[::-1] for f in source_map.fetch_map.keys())}'


    # no request chunking for non-temporal data 
    if src == 'chs':  
        qry = kwargs.copy()
        for k in ('start', 'end', 'top', 'bottom', 'lock'):
            if k in qry.keys(): del qry[k]  # trim hash indexing entropy
        # TODO: split into 1-degree queries for better indexing
        source_map.fetch_map[f'{var}_{src}'](**qry.copy())
        return

    # break request into gridded dx*dy*dt chunks for querying
    xstep = lambda x, bound: int(x - (x % (-bound * dx)))
    ystep = lambda x, bound: int(x - (x % (-bound * dy)))
    kwargs['west'] = max(-180, xstep(kwargs['west'], -1))
    kwargs['east'] = min(+180, xstep(kwargs['east'], +1))
    kwargs['south'] = max(-90, ystep(kwargs['south'], -1))
    kwargs['north'] = min(+90, ystep(kwargs['north'], +1))

    # fetch data chunks
    t = datetime(kwargs['start'].year, kwargs['start'].month, kwargs['start'].day)
    while t < kwargs['end']:
        for x in range(kwargs['west'], kwargs['east'], dx):
            for y in range(kwargs['south'], kwargs['north'], dy):

                qry = {k : v for k,v in zip(
                        ('west', 'east', 'south', 'north', 'start', 'end'),
                        (x, x+dx, y, y+dy, t, t+dt))}

                if 'top' in kwargs.keys():  # get entire depth column
                    qry['top'] = 0  # kwargs['top']
                    qry['bottom'] = 5000  # kwargs['bottom']

                if not serialized(qry, f'fetch_{src}_{var}'):
                    source_map.fetch_map[f'{var}_{src}'](**qry.copy())
                else:
                    logging.debug(f'FETCH_HANDLER DEBUG MSG: already fetched '
                            f'{src}_{var} {t.date().isoformat()}! continuing...')
        t += dt
    
    return 

