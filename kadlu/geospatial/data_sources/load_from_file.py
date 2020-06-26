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


def load_raster(filepath, plot=False, cmap=None, **kwargs):
    """ load data from raster file 

        args:
            filepath: string
                complete filepath descriptor of netcdf file to be read
            plot: boolean
                if True, a plot will be displayed using the qt5agg backend
            cmap: matplotlib colormap
                the full list of available colormaps can be viewed with:
                print(matplotlib.pyplot.colormaps())
                if None is supplied, pyplot will default to 
                matplotlib.pyplot.cm.cividis

        returns:
            values: numpy 2D array
            lats:   numpy 1D array
            lons:   numpy 1D array
    """

    # load raster
    Image.MAX_IMAGE_PIXELS = 500000000
    im = Image.open(filepath)

    # GDAL raster format
    # http://duff.ess.washington.edu/data/raster/drg/docs/geotiff.txt
    if 33922 in im.tag.tagdata.keys():
        i,j,k,x,y,z = im.tag_v2[33922]  # ModelTiepointTag
        dx, dy, dz  = im.tag_v2[33550]  # ModelPixelScaleTag
        meta        = im.tag_v2[42112]  # GdalMetadata
        xml         = ET.fromstring(meta)
        params      = {tag.attrib['name'] : tag.text for tag in xml}
        lat = np.arange(y, y + (dy * im.size[1]), dy)[ :: -1] -90
        rng_lat = (abs(index(kwargs['north'], lat[::-1])-len(lat)), 
                   abs(index(kwargs['south'], lat[::-1])-len(lat)))
        lon = np.arange(x, x + (dx * im.size[0]), dx)
        rng_lon = index(kwargs['east'], lon), index(kwargs['west'], lon)
        logging.info(f'{xml.tag}\nraster coordinate system: {im.tag_v2[34737]}'
                     f'\n{json.dumps(params, indent=2, sort_keys=True)}')

    # NASA / jet propulsion labs raster format (page 27)
    # https://landsat.usgs.gov/sites/default/files/documents/geotiff_spec.pdf
    elif 34264 in im.tag.tagdata.keys():
        dx,_,_,x,_,dy,_,y,_,_,dz,z,_,_,_,_ = im.tag_v2[34264]  # ModelTransformationTag
        lat = np.arange(y, y + (dy * im.size[1]), dy)
        rng_lat = index(kwargs['south'], -lat),  index(kwargs['north'], -lat)
        lon = np.arange(x, x + (dx * im.size[0]), dx)
        rng_lon = index(kwargs['west'], lon), index(kwargs['east'], lon)

    else: assert False, 'unknown metadata tag encoding'
    assert not (z or dz), '3D rasters not supported yet'

    # construct grid and decode pixel values
    if reduce(np.multiply, (rng_lon[1] - rng_lon[0], rng_lat[1] - rng_lat[0])) > 10000000: 
        logging.info('this could take a few moments...')
    grid = np.ndarray((len(np.arange(rng_lon[0], rng_lon[1])), len(np.arange(rng_lat[0], rng_lat[1]))))
    for yi in np.arange(rng_lon[0], rng_lon[1]) -rng_lon[0]: 
        grid[yi] = np.array(list(map(im.getpixel, zip(
            map(int, [yi for xi in np.arange(rng_lon[0], rng_lon[1]) -rng_lon[0]]), 
            map(int,               np.arange(rng_lat[0], rng_lat[1]) -rng_lat[0] )
        ))))
    mask = grid == float(im.tag_v2[42113])
    val = np.ma.MaskedArray(grid, mask=mask)

    # plot the data
    if plot:
        x1, y1 = np.meshgrid(lon[rng_lon[0]:rng_lon[1]], lat[rng_lat[0]:rng_lat[1]], indexing='ij')
        fig = plt.figure()
        if (rng_lon[1]-rng_lon[0]) * (rng_lat[1]-rng_lat[0]) >= 100000:
            ax = fig.add_subplot(1,1,1, projection='scatter_density')
            plt.axis('scaled')
            raster = ax.scatter_density(x1, y1, c=val, cmap=cmap)
            plt.tight_layout()
        else:
            ax = fig.add_subplot(1,1,1)
            ax.scatter(x1, y1, c=val, cmap=cmap)
        plt.show()

    #return val, y1, x1
    return val, lat[rng_lat[0]:rng_lat[1]], lon[rng_lon[0]:rng_lon[1]]


def load_netcdf(filename, var=None, plot=False, cmap=None, **kwargs):
    """ read environmental data from netcdf and output to gridded numpy array

        args:
            filename: string
                complete filepath descriptor of netcdf file to be read
            var: string (optional)
                the netcdf attribute to be read as the values.
                by default, a guess will be made based on the file metadata
            plot: boolean
                if True, a plot will be displayed using the qt5agg backend
            cmap: matplotlib colormap
                the full list of available colormaps can be viewed with:
                print(matplotlib.pyplot.colormaps())
                if None is supplied, pyplot will default to 
                matplotlib.pyplot.cm.cividis

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

    # plot the data
    if plot:
        x1, y1 = np.meshgrid(lon[rng_lon[0]:rng_lon[1]], lat[rng_lat[0]:rng_lat[1]], indexing='ij')
        fig = plt.figure()
        if (rng_lon[1]-rng_lon[0]) * (rng_lat[1]-rng_lat[0]) >= 100000:
            ax = fig.add_subplot(1,1,1, projection='scatter_density')
            plt.axis('scaled')
            raster = ax.scatter_density(x1, y1, c=val, cmap=cmap)
            plt.tight_layout()
        else:
            ax = fig.add_subplot(1,1,1)
            ax.scatter(x1, y1, c=val, cmap=cmap)
        plt.show()

    return val, lat, lon

