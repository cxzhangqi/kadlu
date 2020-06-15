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

raster_table = lambda var: f'raster_{var}'


#def load_files(var, filenames, **kwargs):
""" this function will accept a list of files as string, and determine if 
        they have been extracted and processed already.
        if not, it will determine if they should be processed as netcdf
        or geotiff.

        the content within query bounds will then be returned
"""
"""

    if isinstance(filenames, str): filenames = [filenames]

    for fname in filenames:
        if not serialized(kwargs, fname): 

            # TODO:
            # process 3D data files

            if fname[-3:] == '.nc': process_netcdf_2D(var, fname, kwargs)

            elif fname[-3:] == '.tiff': process_rasters_2D(var, fname, kwargs)

            else: raise NotImplementedError('filetype not supported. valid types include .nc and .tiff')

            insert_hash(kwargs, fname)

            
    sql = ' AND '.join([f'SELECT * FROM {var} WHERE lat >= ?', 
        'lat <= ?',
        'lon >= ?',
        'lon <= ?',
        'time >= ?',
        'time <= ?',
        'depth >= ?',
        'depth <= ?'])

    db.execute(sql, tuple(map(str, [
            kwargs['south'],                kwargs['north'], 
            kwargs['west'],                 kwargs['east'],
            dt_2_epoch(kwargs['start']), dt_2_epoch(kwargs['end']),
            kwargs['top'],                  kwargs['bottom']
        ])))
    rowdata = np.array(db.fetchall(), dtype=object).T

    if len(rowdata) == 0:
        logging.warning(f'{fname} {var}: no data found in region {fmt_coords(kwargs)}, returning empty arrays')
        return np.array([[],[],[],[],[]])
        
    return rowdata[0:5].astype(float)
"""


def process_rasters_2D(var, filepath, kwargs=dict(south=-90, west=-180, north=90, east=180, top=0, bottom=50000, step=0.1)):
    """ process and store arbitrary 2D data from raster format """
    """
    var = 'bathymetry'
    filepath = storage_cfg() + 'bathy_2002.tiff'

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
    source = ['chs' for z in z2]
    grid = list(map(tuple, np.vstack((z2, y2, x2, source)).T))

    # insert into db
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


def load_netcdf_2D(filename, var=None, **kwargs):
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

