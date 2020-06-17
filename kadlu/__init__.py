# data utils
from .geospatial.data_sources.data_util import era5_cfg
from .geospatial.data_sources.data_util import storage_cfg 
from .geospatial.data_sources.data_util import database_cfg 
from .geospatial.data_sources.data_util import epoch_2_dt
from .geospatial.data_sources.data_util import dt_2_epoch
from .geospatial.data_sources.data_util import index
from .geospatial.data_sources.data_util import reshape_2D
from .geospatial.data_sources.data_util import reshape_3D

# load data from file
from .geospatial.data_sources.load_from_file import load_netcdf_2D

# loading with automatic fetching
from .geospatial.data_sources.source_map import source_map 
from .geospatial.data_sources.chs import Chs as chs
from .geospatial.data_sources.era5 import Era5 as era5
from .geospatial.data_sources.gebco import Gebco as gebco
from .geospatial.data_sources.hycom import Hycom as hycom
from .geospatial.data_sources.wwiii import Wwiii as wwiii

# high-level data loading API
from .geospatial.data_sources.source_map import load_map
def load(source, var, **kwargs):

    if var == 'bathymetry' or var == 'depth' or var == 'elevation': var = 'bathy'

    loadkey = f'{var}_{source}'

    assert loadkey in load_map.keys(), 'invalid variable or source. valid options include: \n\n'\
            f'{list(f.rsplit("_", 1)[::-1] for f in load_map.keys())}\n\n'\
            f'for more info, print(kadlu.source_map)'

    return load_map[loadkey](**kwargs)

