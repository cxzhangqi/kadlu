""" Module within the kadlu package for fetching geospatial data 

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
from kadlu.geospatial.bathy_reader import LatLon

def fetch_bathy_chs(latlon_SW=LatLon(-90,-180), latlon_NE=LatLon(90,180)):
    """Fetch the Canadian Hydrographic Service Non-Navigational (NONNA-100) Bathymetric Data products.

        Args: 
            latlon_SW: LatLon
                South-western (SW) boundary of the region of interest.
            latlon_NE: LatLon
                North-eastern (SE) boundary of the region of interest.

        Returns:
            paths: list
                Paths to the data files that were retrieved.
    """
    fnames = ["CA2_4300N06000W.tif",  "CA2_4400N06000W.tif"]
    paths = "/home/oliskir/src/meridian/kadlu/kadlu/tests/assets/tif" + fnames
    return paths