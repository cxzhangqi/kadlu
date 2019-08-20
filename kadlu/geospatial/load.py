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


def load_bathy(storage_location, south=-90, north=90, west=-180, east=180, source="CHS"):
    """ Load bathymetry data within specified geographical region.

        Possible data sources are: CHS, GEBCO

        TODO: Get rid of the storage_location argument and instead use the config.ini file

        Args: 
            south: float
                Southern boundary of the region of interest.
            north: float
                Northern boundary of the region of interest.
            west: float
                Western boundary of the region of interest.
            east: float
                Eastern boundary of the region of interest.
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
        bathy, lats, lons = chs.load(storage_location, south, north, west, east)

    elif source == "GEBCO":
        bathy, lats, lons = gebco.load(storage_location, south, north, west, east)

    else:
        print('Unknown bathymetry data source')
        exit(1)

    return bathy, lats, lons