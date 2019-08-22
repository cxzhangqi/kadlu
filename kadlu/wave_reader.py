""" Wave data reader module within the kadlu package

    Authors: Casey Hilliard
    contact: r.casey.hilliard@gmail.com
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from enum import Enum
import numpy as np
import pygrib

from kadlu.utils import LatLon

class Source(Enum):
    """ Enum class for wave data sources
    """
    RDWPS_St_Lawrence = 1
    GDWPS = 2
    NOAA_Wavewatch_Global = 3
    ECMWF_ERA5 = 4

class WaveReader():

    """ Class for reading Wave data
    
        Given a target file of GRIB data, containing wave parameters (fore- or
        hind- cast, or observed ), a timestamp, and grid defining parameters, 
        constructs an interpolated, regularly gridded representation.
    
    
        Attributes: 
            input: str
                The path to an input (GRIB2) datafile, containing
                data to span the required spatio-temporal range required.
            wave_parameter: str
                Name of the variable that contains the wave parameter values.
            lat_name: str
                Name of the variable that contains the latitude values.
            lon_name: str
                Name of the variable that contains the longitude values.

    """
    def __init__(self, input_files, source, lat_name='nj', lon_name='ni', param_name='HTSGW', lon_axis=1):

        # get files
        self.files = self._get_files(input_files)
        assert len(self.files) > 0, "You must provide at least one data file" 

        # file format
        self.format = self._detect_format(self.files)

        self.source = source
        self.param_name = param_name
        self.lon_axis = lon_axis    

    def _read_grib(self, date, date_interval, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
        """ Read longitude, latitude, and wave parameter matrices from file, given a specific timestamp and
            region.

            Args: 
                date: Timestamp
                    Timestamp of data to be extracted from GRIB2 input file(s).
                date_interval: Timedelta
                    Interval of time over which GRIB file is stepped.
                latlon_SW: LatLon
                    South-western (SW) boundary of the region of interest.
                latlon_NE: LatLon
                    North-eastern (NE) boundary of the region of interest.

            Returns:
                lat: 1d numpy array
                    Latitude values
                lon: 1d numpy array
                    Longitude values
                param: 2d numpy array
                    Parameter values
        """
        # load data
        grib = Dataset(self.files[0])

        grbs=pygrib.open(grib)
        param_values = grbs.select(self.param_name)[0]
        bathy = np.array(d.variables[self.bathy_name])
        
        lat = np.array(d.variables[self.lat_name])
        lon = np.array(d.variables[self.lon_name])

        # select region of interest
        lat, lon, bathy = self._select_region(lat, lon, bathy, latlon_SW, latlon_NE)

        return lat, lon, bathy
