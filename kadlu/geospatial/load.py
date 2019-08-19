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
import kadlu.geospatial.bathy_chs as chs
from kadlu.geospatial.bathy_reader import LatLon


def load_bathy(storage_location, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180), source="CHS"):
    """ Load bathymetry data within specified geographical region.

        TODO: Get rid of the storage_location argument and instead use the config.ini file

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.
            source: str
                Bathymetry data source(s).               

        Returns:
            bathy: 1d numpy array
                Bathymetry values
            lats: 1d numpy array
                Latitude values
            lons: 1d numpy array
                Longitude values
    """
    assert source == "CHS", "The only available bathymetry data source is CHS"

    # fetch relevant data files
    files = chs.fetch(latlon_SW, latlon_NE)

    bathy, lats, lons = list(), list(), list()        

    # loop over geotiff files
    for f in files:

        # read data from geotiff file
        v = chs.read(path=f)

        # create lat-lon arrays
        lat, lon = chs.latlon(path=f)

        lats.append(lat)
        lons.append(lon)
        bathy.append(v)

    # concatenate
    bathy = np.concatenate(bathy)
    lats = np.concatenate(lats)
    lons = np.concatenate(lons)

    return bathy, lats, lons
