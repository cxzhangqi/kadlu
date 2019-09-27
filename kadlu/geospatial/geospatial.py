""" Geospatial utilities module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import numpy as np
from osgeo import gdal
from netCDF4 import Dataset
import scipy.io as sio
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")

from matplotlib import pyplot as plt


def read_geotiff(path, band_id=1):
    """ Read data from a GeoTIFF file.

        Args: 
            path: str
                File path
            band_id: int
                Number of the GeoTIFF raster band containing the data values   

        Returns:
            values: masked numpy array
                Data values. The array has been masked where invalid values occur (NaNs or infs).
    """
    # load data
    data_set = gdal.Open(path)

    # access values
    band = data_set.GetRasterBand(band_id)
    values = data_set.ReadAsArray()

    # replace no-data value with nan
    nodata = band.GetNoDataValue()
    values[values == nodata] = np.nan
    values = np.ma.masked_invalid(values)

    return values


def crop(lat, lon, south, north, west, east, grid=False):
    """ Select rectangular region bounded by the geographical coordinates 
        latlon_SW to the south-west and latlon_NE to the north-east.

        If grid is False, lat and lon must have the same length.

        Args: 
            lat: 1d or 2d numpy array
                Latitude values
            lon: 1d or 2d numpy array
                Longitude values
            south: float
                Southern boundary of the region of interest.
            north: float
                Northern boundary of the region of interest.
            west: float
                Western boundary of the region of interest.
            east: float
                Eastern boundary of the region of interest.
            grid: bool
                Specify how to combine elements of lat and lon.

        Returns:
            ind: numpy array
                Selected indices. 1d if grid is False, 2d if grid is True.
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
    """
    ind_lat = np.argwhere((lat >= south) & (lat <= north))
    ind_lat = np.squeeze(ind_lat)

    ind_lon = np.argwhere((lon >= west) & (lon <= east))
    ind_lon = np.squeeze(ind_lon)

    if np.ndim(ind_lat) == 2:
        latc = ind_lat[:,0] + 1j * ind_lat[:,1]
        lonc = ind_lon[:,0] + 1j * ind_lon[:,1]
        x = np.intersect1d(latc, lonc)
        xr = np.real(x)
        xi = np.imag(x)
        xr = xr.astype(int)
        xi = xi.astype(int)
        ind_lat = np.arange(np.min(xr), np.max(xr)+1)
        ind_lon = np.arange(np.min(xi), np.max(xi)+1)
        lat = lat[:,0]
        lon = lon[0,:]
        ind = np.ix_(ind_lat, ind_lon)       

    if grid:
        ind = np.ix_(ind_lat, ind_lon)       
    else:
        ind = np.intersect1d(ind_lat, ind_lon)
        ind_lat = ind
        ind_lon = ind

    lat = lat[ind_lat]
    lon = lon[ind_lon]

    return ind, lat, lon


def load_data_from_file(path, val_name='bathy', lat_name='lat', lon_name='lon', lon_axis=1,\
    south=-90, north=90, west=-180, east=180):
    """ Load geospatial data from a single file. 

        Currently supported formats are NetCDF (*.nc) and MatLab (*.mat).

        The data can be cropped by speciyfing south/north/west/east 
        boundaries.

        Args: 
            path: str
                File path
            val_name: str
                Name of variable/field containing the data values
            lat_name: str
                Name of variable/field containing the latitude values
            lon_name: str
                Name of variable/field containing the longitude values
            lon_axis: int
                Specify if the longitude dimension is the second (1, default) 
                or first (0) axis.
            south: float
                Southern boundary of the region of interest.
            north: float
                Northern boundary of the region of interest.
            west: float
                Western boundary of the region of interest.
            east: float
                Eastern boundary of the region of interest.

        Returns:
            val: 1d or 2d numpy array
                Data values
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
    """
    # detect format
    ext = path[path.rfind('.'):]

    # load data
    if ext == '.nc': # NetCDF
        d = Dataset(path)
        val = np.array(d[val_name])
        lat = np.array(d[lat_name])
        lon = np.array(d[lon_name])

    elif ext == '.mat': # MatLab
        d = sio.loadmat(path)
        val = np.array(d[val_name])
        lat = np.squeeze(np.array(d[lat_name]))
        lon = np.squeeze(np.array(d[lon_name]))

    else:
        print('Unknown file format *{0}'.format(ext))
        exit(1)

    # crop the region of interest
    grid = (np.ndim(val) == 2)
    indices, lat, lon = crop(lat, lon, south, north, west, east, grid=grid)
    val = val[indices]

    # ensure that lat and lon are strictly increasing
    if np.all(np.diff(lat) < 0):
        lat = np.flip(lat, axis=0)
        val = np.flip(val, axis=0)
    if np.all(np.diff(lon) < 0):
        lon = np.flip(lon, axis=0)
        val = np.flip(val, axis=1)

    # flip axes, if necessary
    if lon_axis == 0:
        val = np.swapaxes(val, 0, 1)

    # if axis size are inconsistent, try swapping
    if val.shape[0] != lat.shape[0]:
        val = np.swapaxes(val, 0, 1)

    return val, lat, lon


def write_data_to_file(lat, lon, val, destination, compression=False):
    """ Write latitude, longitude, and geospatial data to a file.

        The current implementation only supports MatLab format (*.mat)

        Returns AssertionError if the destination path does not have *.mat extension.

        Args:
            lat: 1d numpy array
                Latitude values
            lon: 1d numpy array
                Longitude values
            val: 1d or 2d numpy array
                Geospatial data values
            destination: str
                Name of output file. Must have extension *.mat
            compression: bool
                Compress matrices on write. Default is False
    """
    # parse file format
    p = destination.rfind('.')
    ext = destination[p:]
    
    assert ext == '.mat', 'Destination file must have extension *.mat (MATLAB file)'

    output = {'lat': lat, 'lon': lon, 'data': val}
    sio.savemat(file_name=destination, mdict=output, do_compression=compression)

    print('Data saved to ' + destination)


def plot(x, y, z, x_label='x (m)', y_label='y (m)', z_label='Elevation (m)'):
    """ Plot 2d geospatial data by drawing a color heat map.

        Args:
            x: 1d numpy array
                x-coordinates or longitudes
            y: 1d numpy array
                y-coordinates or latitudes
            z: 2d numpy array
                data values
            x_label: str
                x-axis label
            y_label: str
                y-axis label
            z_label: str
                z-axis label

        Returns:
            fig: matplotlib.figure.Figure
                A figure object.
    """
    # axes ranges
    x_min = np.min(x)
    x_max = np.max(x)
    y_min = np.min(y)
    y_max = np.max(y)

    # meshgrid
    x,y = np.meshgrid(x,y)

    # plot
    fig, ax = plt.subplots(figsize=(8,6))
    img = ax.imshow(z.T, aspect='auto', origin='lower', extent=(x_min, x_max, y_min, y_max))

    # axes titles
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    # Add a color bar which maps values to colors
    fig.colorbar(img, format='%.02f', label='Elevation (m)')

    return fig


