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
import kadlu.geospatial.data_sources.chs as chs
import kadlu.geospatial.data_sources.gebco as gebco 
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
        bathy, lats, lons = chs.load(storage_location, latlon_SW, latlon_NE)

    elif source == "GEBCO":
        bathy, lats, lons = gebco.load(storage_location, latlon_SW, latlon_NE)

    else:
        print('Unknown bathymetry data source')
        exit(1)

    return bathy, lats, lons