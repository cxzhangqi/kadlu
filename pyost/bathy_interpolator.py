""" Bathymetry interpolation module within the pyost package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/pyost
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import numpy as np
from collections import namedtuple
from scipy.interpolate import RectBivariateSpline
from pyost.bathy_reader import BathyReader, LatLon


class BathyInterpolator():
    """ Class for interpolating bathymetry data.

        Args: 
            bathy_reader: BathyReader
                Bathymetry data file reader
            latlon_SW: LatLon
                South-western (SW) boundary of the interpolation region.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the interpolation region.
    """
    def __init__(self, bathy_reader, latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180), origin=None):
        self.origin = origin
    
        lat, lon, bathy = bathy_reader.read(latlon_SW, latlon_NE)

        self.interp_ll = RectBivariateSpline(x=lat, y=lon, z=bathy)
    
        self.interp_xy = None

    def eval_xy(self, x, y, grid=False):
        """ Evaluate interpolated bathymetry in position coordinates (XY).

            Returns:
                zi: Interpolated bathymetry values
        """
        assert self.origin is not None, 'Evaluation by position coordinates requires that the origin has been specified.'
        assert self.interp_xy is not None, 'xy interpolation not available'

        zi = self.interp_xy.__call__(x=x, y=y, grid=grid)
        return zi

    def eval_ll(self, lat, lon, grid=False):
        """ Interpolate bathymetry grid in latitude and longitude coordinates (LL).

                Returns:
                    zi: Interpolated bathymetry values
        """
        zi = self.interp_ll.__call__(x=lat, y=lon, grid=grid)
        return zi