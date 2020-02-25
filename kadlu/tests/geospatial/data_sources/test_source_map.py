from datetime import datetime
from kadlu.geospatial.data_sources.source_map import fetch_handler
kwargs = dict(
        start=datetime(2015, 3, 1), end=datetime(2015, 3, 3),
        south=45,                   west=-68.4, 
        north=51.5,                 east=-56.5, 
        top=0,                      bottom=100,
    )


def test_batch_wwiii():
    fetch_handler('waveheight', 'wwiii', parallel=4, **kwargs)

def test_batch_hycom():
    fetch_handler('temp', 'hycom', parallel=2, **kwargs)

def test_batch_era5():
    fetch_handler('wavedir', 'era5', parallel=4, **kwargs)

def test_batch_chs():
    # this is just a wrapper for Chs().fetch_bathymetry
    # with an additional hash index check.
    # non-temporal data not parallelized
    fetch_handler('bathy', 'chs', parallel=9999, 
            south=45, west=-67, north=46, east=-66,
            start=datetime.now(), end=datetime.now())

""" interactive testing

    var='significant_height_of_combined_wind_waves_and_swell'
    kwargs = dict(
        start=datetime(2015, 3, 1), end=datetime(2015, 3, 3),
        south=45,                   west=-68.4, 
        north=51.5,                 east=-56.5, 
        top=0,                      bottom=5000
        )

"""

