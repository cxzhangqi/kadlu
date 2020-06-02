from datetime import datetime
from kadlu.geospatial.data_sources.fetch_handler import fetch_handler

kwargs = dict(
        start=datetime(2015, 2, 1), end=datetime(2015, 2, 1, 12),
        south=47,                   west=-61, 
        north=50,                   east=-58, 
        top=0,                      bottom=5000,
    )


def test_batch_wwiii():
    fetch_handler('wind_uv', 'wwiii', **kwargs)

def test_batch_hycom():
    fetch_handler('salinity', 'hycom', **kwargs)

def test_batch_era5():
    fetch_handler('wavedir', 'era5', **kwargs)

def test_batch_chs():
    fetch_handler('bathy', 'chs', south=45, west=-67, north=46, east=-66)

""" interactive testing


    kwargs = dict(
        start=datetime(2015, 3, 1), end=datetime(2015, 3, 3),
        south=45,                   west=-68.4, 
        north=51.5,                 east=-56.5, 
        top=0,                      bottom=5000
        )

"""

