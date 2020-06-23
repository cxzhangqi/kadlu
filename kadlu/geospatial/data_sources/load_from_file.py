import logging
from PIL import Image
from functools import reduce
from xml.etree import ElementTree as ET
import json

import matplotlib
matplotlib.use('qt5agg')
import mpl_scatter_density
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

from kadlu.geospatial.data_sources.data_util        import          \
        database_cfg,                                               \
        storage_cfg,                                                \
        insert_hash,                                                \
        serialized,                                                 \
        index


def load_raster(filepath, plot=True, cmap=None, **kwargs):
    """ load data from raster file """
    """
    #var = 'bathymetry'
    filepath = storage_cfg() + 'gebco_2020_n0.0_s-90.0_w-180.0_e-90.0.tif'
    filepath = storage_cfg() + 'test.tif'
    filepath = storage_cfg() + 'GEBCO_BATHY_2002-01-01_rgb_3600x1800.TIFF'
    filepath = storage_cfg() + 'BlueMarbleNG_2004-12-01_rgb_3600x1800.TIFF'
    kwargs=dict(south=-90, west=-180, north=90, east=180, top=0, bottom=50000)
    """
    
    # suppress decompression bomb error and open raster
    Image.MAX_IMAGE_PIXELS = 500000000
    im = Image.open(filepath)
    nan = float(im.tag_v2[42113])

    # GDAL raster format
    # http://duff.ess.washington.edu/data/raster/drg/docs/geotiff.txt
    if 33922 in im.tag.tagdata.keys():
        i,j,k,x,y,z = im.tag_v2[33922]  # ModelTiepointTag
        dx, dy, dz  = im.tag_v2[33550]  # ModelPixelScaleTag
        meta        = im.tag_v2[42112]  # GdalMetadata
        tree        = ET.fromstring(meta)
        params      = {entry.attrib['name'] : entry.text for entry in tree}
        logging.info(f'{tree.tag}\nraster coordinate system: {im.tag_v2[34737]}\n{json.dumps(params, indent=2, sort_keys=True)}')
        lat = np.arange(y, y + (dy * im.size[1]), dy)[ :: -1]

    # NASA / jet propulsion labs raster format
    # https://landsat.usgs.gov/sites/default/files/documents/geotiff_spec.pdf
    elif 34264 in im.tag.tagdata.keys():
        a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p = im.tag_v2[34264]  # ModelTransformationTag
        dx,x,dy,y,dz,z = a,d,f,h,k,l  # refer to page 27 for transformation matrix
        lat = np.arange(y, y + (dy * im.size[1]), dy)

    else: assert False, 'unknown metadata tag encoding'
    assert not (z or dz), 'TODO: implement 3D raster support'

    # construct grid and decode pixel values
    if reduce(np.multiply, im.size) > 10000000: logging.info('this could take a few moments...')
    lon = np.arange(x, x + (dx * im.size[0]), dx)
    grid = np.ndarray(list(map(int, reduce(np.append, (im.size,np.array(im.getpixel((0,0))).shape)))))
    for yi in range(im.size[0]): grid[yi] = np.array(list(map(im.getpixel, zip([yi for xi in range(im.size[0])], range(im.size[1])))))
    mask = grid == nan
    val = np.ma.MaskedArray(grid, mask=mask)
    x1, y1 = np.meshgrid(lon, lat, indexing='ij')

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection='scatter_density')
        ax.set_title(filepath)
        density = ax.scatter_density(x1, y1, c=val, cmap=cmap)
        fig.colorbar(density, label='pixel value')
        plt.tight_layout()
        plt.show()
    
    return val, y1, x1


def load_netcdf(filename, var=None, **kwargs):
    """ read environmental data from netcdf and output to gridded numpy array

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

"""

filepath = storage_cfg() + 'gebco_2020_n0.0_s-90.0_w-180.0_e-90.0.tif'
filepath = storage_cfg() + 'test.tif'
filepath = storage_cfg() + 'GEBCO_BATHY_2002-01-01_rgb_3600x1800.TIFF'
filepath = storage_cfg() + 'BlueMarbleNG_2004-12-01_rgb_3600x1800.TIFF'

load_raster(filepath)

"""
