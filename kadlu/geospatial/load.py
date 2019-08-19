""" Module within the kadlu package for loading geospatial data 

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import numpy as np
from kadlu.geospatial.read import read_netcdf
import kadlu.geospatial.bathy_chs as chs
from kadlu.geospatial.utils import crop
from kadlu.geospatial.bathy_reader import LatLon


def load_bathy(storage_location, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180), source="CHS"):
    """ Load bathymetry data within specified geographical region.

        Possible data sources are: CHS, GEBCO

        TODO: Get rid of the storage_location argument and instead use the config.ini file

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.
            source: str
                Bathymetry data source(s).               

        Returns:
            bathy: numpy array
                Bathymetry values
            lats: numpy array
                Latitude values
            lons: numpy array
                Longitude values
    """
    if source == "CHS":
        bathy, lats, lons = load_bathy_chs(storage_location, latlon_SW, latlon_NE)

    elif source == "GEBCO":
        bathy, lats, lons = load_bathy_gebco(storage_location, latlon_SW, latlon_NE)

    else:
        print('Unknown bathymetry data source')
        exit(1)

    # crop the region of interest
    indices, lats, lons = crop(lats, lons, latlon_SW, latlon_NE)
    bathy = bathy[indices]

    return bathy, lats, lons


def load_bathy_chs(storage_location, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
    """ Load Non-Navigational NONNA-100 bathymetric data from the Canadian Hydrographic 
        Service (CHS) within specified geographical region.

        TODO: Get rid of the storage_location argument and instead use the config.ini file

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.

        Returns:
            bathy: 1d numpy array
                Bathymetry values
            lats: 1d numpy array
                Latitude values
            lons: 1d numpy array
                Longitude values
    """
    # fetch relevant data files
    files = chs.fetch(storage_location, latlon_SW, latlon_NE)

    bathy, lats, lons = list(), list(), list()        

    # loop over geotiff files
    for f in files:

        # read data from geotiff file
        z = chs.read(path=f)

        # create lat-lon arrays
        lat, lon = chs.latlon(path=f)

        # make a grid
        x, y = np.meshgrid(lon, lat)

        # select non-masked entries
        x = x[~z.mask]
        y = y[~z.mask]
        z = z[~z.mask]

        lats.append(y)
        lons.append(x)
        bathy.append(z)

    # concatenate
    bathy = np.ma.concatenate(bathy)
    lats = np.concatenate(lats)
    lons = np.concatenate(lons)

    return bathy, lats, lons


def load_bathy_gebco(path, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
    """ Load GEBCO bathymetry data within specified geographical region.

        TODO: Get rid of the path argument and instead use the config.ini file

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.

        Returns:
            bathy: 2d numpy array
                Bathymetry values
            lats: 1d numpy array
                Latitude values
            lons: 1d numpy array
                Longitude values
    """
    # read data
    bathy, lats, lons = read_netcdf_2d(path=path, val_name="bathy", lat_name="lat", lon_name="lon")        

    return bathy, lats, lons


def crop_bathy(bathy, lat, lon, latlon_SW, latlon_NE):
    """ Select rectangular region bounded by the geographical coordinates 
        latlon_SW to the south-west and latlon_NE to the north-east.

        Latitudes and longitudes can be given as 2d numpy arrays, but
        must be in regularly spaced grid.

        Args: 
            bathy: 2d numpy array
                Bathymetry values
            lat: 1d or 2d numpy array
                Latitude values
            lon: 1d or 2d numpy array
                Longitude values
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (NE) boundary of the region of interest.

        Returns:
            bathy: numpy array
                Bathymetry values
            lat: numpy array
                Latitude values
            lon: numpy array
                Longitude values
    """
    ind_lat = np.argwhere((lat >= latlon_SW.latitude) & (lat <= latlon_NE.latitude))
    ind_lat = np.squeeze(ind_lat)

    ind_lon = np.argwhere((lon >= latlon_SW.longitude) & (lon <= latlon_NE.longitude))
    ind_lon = np.squeeze(ind_lon)

    if np.ndim(ind_lat) == 1:
        ind = np.ix_(ind_lat, ind_lon)       

    else:
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

    bathy = bathy[ind]
    lat = lat[ind_lat]
    lon = lon[ind_lon]

    return bathy, lat, lon