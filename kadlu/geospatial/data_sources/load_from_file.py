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
    """ load data from raster file """
    """
    #var = 'bathymetry'
    filepath = storage_cfg() + 'gebco_2020_n90.0_s0.0_w-90.0_e0.0.tif'
    filepath = storage_cfg() + 'test.tif'
    filepath = storage_cfg() + 'GEBCO_BATHY_2002-01-01_rgb_3600x1800.TIFF'
    filepath = storage_cfg() + 'BlueMarbleNG_2004-12-01_rgb_3600x1800.TIFF'

    plot=True
    cmap=None
    #kwargs=dict(south=-90, west=-180, north=90, east=180, top=0, bottom=50000)
    kwargs=dict(south=60, west=-, north=61, east=180, top=0, bottom=50000)

    kwargs['north'], kwargs['east'] = 60.4709, 61.26033
    kwargs['south'], kwargs['west'] = 60.46333,61.15905

    kwargs['north'], kwargs['east'] = 62.4709,  63.26033
    kwargs['south'], kwargs['west'] = 60.46333, 61.15905

We have two moorings 2019-2020  August 2019 to end September 2020  (yes, some futures in there as well).
  Mooring HiBioA 2019-2020 … located at:   60 degrees, 28.254 minutes N,  61 degrees, 15.620 W
  Mooring HiBioC 2019-2020 … located at:   60 degrees, 27.7998 minutes N,  61 degrees, 09.543 W

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
        logging.info(f'{tree.tag}\nraster coordinate system: {im.tag_v2[34737]}\n{json.dumps(params, indent=2, sort_keys=True)}')
        lat = np.arange(y, y + (dy * im.size[1]), dy)[ :: -1]
        #rng_lat = index(kwargs['west'], lat),  index(kwargs['east'], lat)
        rng_lat = index(kwargs['south'], lat),  index(kwargs['north'], lat)

    # NASA / jet propulsion labs raster format (page 27)
    # https://landsat.usgs.gov/sites/default/files/documents/geotiff_spec.pdf
    elif 34264 in im.tag.tagdata.keys():
        dx,_,_,x,_,dy,_,y,_,_,dz,z,_,_,_,_ = im.tag_v2[34264]  # ModelTransformationTag
        lat = np.arange(y, y + (dy * im.size[1]), dy)
        rng_lat = index(kwargs['south'], -lat),  index(kwargs['north'], -lat)

    else: assert False, 'unknown metadata tag encoding'

    lon = np.arange(x, x + (dx * im.size[0]), dx)
    rng_lon = index(kwargs['west'], lon), index(kwargs['east'], lon)

    assert not (z or dz), '3D rasters not supported yet'

    # construct grid and decode pixel values
    if reduce(np.multiply, (rng_lon[1] - rng_lon[0], rng_lat[1] - rng_lat[0])) > 10000000: logging.info('this could take a few moments...')
    #grid = np.ndarray(list(map(int, reduce(np.append, (im.size[0],np.array(im.getpixel((0,0))).shape)))))
    grid = np.ndarray((len(np.arange(rng_lon[0], rng_lon[1])), len(np.arange(rng_lat[0], rng_lat[1]))))
    for yi in np.arange(rng_lon[0], rng_lon[1]) -rng_lon[0]: 
        grid[yi] = np.array(list(map(im.getpixel, zip(
            map(int, [yi for xi in np.arange(rng_lon[0], rng_lon[1]) -rng_lon[0]]), 
            map(int,               np.arange(rng_lat[0], rng_lat[1]) -rng_lat[0] )
        ))))
    mask = grid == float(im.tag_v2[42113])
    val = np.ma.MaskedArray(grid, mask=mask)
    x1, y1 = np.meshgrid(lon[rng_lon[0]:rng_lon[1]], -1*lat[rng_lat[0]:rng_lat[1]], indexing='ij')

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection='scatter_density')
        if borderless:
            ax.set_xticks([]), ax.set_yticks([])
            plt.axis('scaled')
            raster = ax.scatter_density(x1, y1, c=val, cmap=cmap, vmin=0, vmax=256)
            #raster = ax.scatter_density(lon, lat, c=val, cmap=cmap, vmin=0, vmax=256)
            plt.tight_layout()
            plt.show()
            #ax.scatter(x1, y1, c=val)
            #fig.savefig(storage_cfg()+'figures/test_raster.png', figsize=(200,200), dpi=300)
        else:
            ax.set_title(filepath)
            raster = ax.scatter_density(x1, y1, c=val, cmap=cmap)
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

