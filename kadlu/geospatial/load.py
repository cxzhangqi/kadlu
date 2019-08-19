""" Module within the kadlu package for loading geospatial data 

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
from kadlu.geospatial.read import read_geotiff_2d
from kadlu.geospatial.fetch import fetch_bathy_chs
from kadlu.geospatial.bathy_reader import LatLon


def load_bathy(latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180), source="CHS"):
    """ Load bathymetry data within specified geographical region.

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.
            source: str
                Bathymetry data source(s).               

        Returns:
            lat: 1d numpy array
                Latitude values
            lon: 1d numpy array
                Longitude values
            bathy: 1d numpy array
                Bathymetry values
    """
    assert source == "CHS", "The only available bathymetry data source is CHS"

    # fetch relevant data files
    files = fetch_bathy_chs(latlon_SW, latlon_NE)

    bathy, lat, lon = list(), list(), list()        

    # loop over geotiff files
    for f in files:

        # read data from geotiff file
        v, _, _ = read_geotiff_2d(path=f)

        # create lat,lon arrays
        x = np.ones(1)
        y = np.ones(1)

        lat.append(x)
        lon.append(y)
        bathy.append(v)

    # concatenate
    bathy = np.concatenate(bathy)
    lat = np.concatenate(lat)
    lon = np.concatenate(lon)

    return bathy, lat, lon
