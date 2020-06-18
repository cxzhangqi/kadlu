import logging
from PIL import Image

import netCDF4
import numpy as np

#from kadlu.geospatial.data_sources.fetch_handler import bin_request
from kadlu.geospatial.data_sources.data_util        import          \
        database_cfg,                                               \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized,                                                 \
        index


conn, db = database_cfg()


def load_raster(filepath, kwargs=dict(south=-90, west=-180, north=90, east=180, top=0, bottom=50000, step=0.1)):
    """ load 2D data from raster file """
    """
    #var = 'bathymetry'
    filepath = storage_cfg() + 'BlueMarbleNG_2004-12-01_rgb_3600x1800.TIFF'
    kwargs=dict(south=-90, west=-180, north=90, east=180, top=0, bottom=50000, step=0.1)

    """
    
    # open image and validate metadata
    im = Image.open(filepath)
    assert im.size == (np.arange(kwargs['west'], kwargs['east'], kwargs['step']).size, np.arange(kwargs['south'], kwargs['north'], kwargs['step']).size), 'metadata does not match data'

    # interpret pixels as elevation
    nan = float(im.tag[42113][0])
    val = np.ndarray((im.size[0], im.size[1]))
    for yi in range(im.size[1]): val[yi] = np.array(list(map(im.getpixel, zip([yi for xi in range(im.size[0])], range(im.size[1])))))
    mask = np.flip(val == nan, axis=0)

    # generate latlon arrays
    lon, lat = np.arange(kwargs['west'], kwargs['east'], kwargs['step']), np.arange(kwargs['south'], kwargs['north'], kwargs['step'])

    # select non-masked entries, remove missing, build grid
    z1 = np.flip(val, axis=0)
    x1, y1 = np.meshgrid(lon, lat)
    x2, y2, z2 = x1[~mask], y1[~mask], np.abs(z1[~mask])
    #source = ['chs' for z in z2]
    #grid = list(map(tuple, np.vstack((z2, y2, x2, source)).T))

    """
    # insert into db
    raster_table = lambda var: f'raster_{var}'
    n1 = db.execute(f"SELECT COUNT(*) FROM {raster_table(var)}").fetchall()[0][0]
    db.executemany(f"INSERT OR IGNORE INTO {raster_table(var)} VALUES (?,?,?,?)", grid)
    n2 = db.execute(f"SELECT COUNT(*) FROM {raster_table(var)}").fetchall()[0][0]
    db.execute("COMMIT")
    conn.commit()
    logging.info(f"RASTER {filepath.split('/')[-1]} {var} in region "
          f"{fmt_coords(dict(south=south,west=west,north=north,east=east))}. "
          f"processed and inserted {n2-n1} rows. "
          f"{len(z1[~mask]) - len(grid)} null values removed, "
          f"{len(grid) - (n2-n1)} duplicate rows ignored")
    """
    
    return z2, y2, x2


def load_netcdf(filename, var=None, **kwargs):
    """ read environmental data from netcdf and output to gridded 2D numpy array

        args:
            filename: string
                complete filepath descriptor of netcdf file to be read
            var: string (optional)
                the netcdf attribute to be read as the values.
                by default, a guess will be made based on the file metadata

        returns:
            values: numpy 2D array
            lats:   numpy 1D array
            lons:   numpy 1D array
    """

    ncfile = netCDF4.Dataset(filename)

    if var is None:
        assert 'lat' in ncfile.variables.keys()
        assert 'lon' in ncfile.variables.keys()
        assert len(ncfile.variables.keys()) == 3
        var = [_ for _ in ncfile.variables.keys() if _ != "lat" and _ != "lon"][0]

    assert var in ncfile.variables, f'variable {var} not in file. file contains {ncfile.variables.keys()}'

    logging.info(f'loading {var} from {ncfile.getncattr("title")}')

    rng_lat = index(kwargs['west'],  ncfile['lat'][:].data), index(kwargs['east'],  ncfile['lat'][:].data)
    rng_lon = index(kwargs['south'], ncfile['lon'][:].data), index(kwargs['north'], ncfile['lon'][:].data)

    val = ncfile[ var ][:].data[rng_lat[0]:rng_lat[1], rng_lon[0]:rng_lon[1]]
    lat = ncfile['lat'][:].data[rng_lat[0]:rng_lat[1]]
    lon = ncfile['lon'][:].data[rng_lon[0]:rng_lon[1]]

    return val, lat, lon

