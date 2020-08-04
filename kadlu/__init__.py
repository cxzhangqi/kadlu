import os
import logging

LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO')
logging.basicConfig(format='%(asctime)s  %(message)s', level=LOGLEVEL, datefmt='%Y-%m-%d %I:%M:%S')

# data utils
from .geospatial.data_sources.data_util import (
        Capturing,
        database_cfg,
        dt_2_epoch,
        epoch_2_dt,
        ext,
        index,
        reshape_2D,
        reshape_3D,
        storage_cfg,
    )

# automatic fetching without loading
from .geospatial.data_sources.ifremer import Ifremer as ifremer 

# loading with automatic fetching
from .geospatial.data_sources.source_map import source_map 
from .geospatial.data_sources.chs import Chs as chs
from .geospatial.data_sources.era5 import Era5 as era5, era5_cfg
from .geospatial.data_sources.gebco import Gebco as gebco
from .geospatial.data_sources.hycom import Hycom as hycom
from .geospatial.data_sources.wwiii import Wwiii as wwiii

# load data from local files
from .geospatial.data_sources.load_from_file import load_netcdf
from .geospatial.data_sources.load_from_file import load_raster

# user-facing data loading API
from .geospatial.data_sources.source_map import load_map

# systematic file testing for all files in kadlu_data/testfiles/
from .tests.geospatial.data_sources.test_files import test_files

# plotting tools
from .plot_util import *


def load(source, var, **kwargs):
    """ automated fetching and loading from web sources 

        args
            source, var: strings
                to view the complete list of sources and variables:
                print(kadlu.source_map)

            kwargs: dictionary
                dict containing boundary coordinates. example:
                kwargs=dict(
                    south=44.25, west=-64.5,
                    north=44.70, east=-63.33,
                    top=0, bottom=5000,
                    start=datetime(2015, 3, 1), end=datetime(2015, 3, 1, 12)
                )

        returns an ND numpy array. arrays ordered by:
            val, lat, lon, [time, depth]
            times are in epoch format
    """

    source, var = source.lower(), var.lower()
    if var == 'bathymetry' or var == 'depth' or var == 'elevation': var = 'bathy'

    loadkey = f'{var}_{source}'
    assert loadkey in load_map.keys(), f'error: invalid source or variable. valid options include: \n\n'\
            f'{list(f.rsplit("_", 1)[::-1] for f in load_map.keys())}\n\n'\
            f'for more info, print(kadlu.source_map)'

    return load_map[loadkey](**kwargs)


def load_file(filepath, **kwargs):
    """ loading from local files 

        args
            filepath: string
                full path directory and filename of data to load

            kwargs: dictionary
                dict containing coordinate bounds, e.g. north, south, 
                east, west, top, bottom, start, end
                note that times are in datetime format

        returns an ND numpy array. arrays ordered by:
            val, lat, lon, [time, depth]
            times are in epoch format
    """
    assert os.path.isfile(filepath), f'error: could not find {filepath}'

    if ext(filepath, ('.nc',)):
        return load_netcdf(filepath, **kwargs)

    elif ext(filepath, ('.tif', '.tiff', '.gtiff',)):       
        return load_raster(filepath, **kwargs)

    else:
        assert False, f'error {filepath}: unknown file format - currently only .nc and .tif formats are accepted'

