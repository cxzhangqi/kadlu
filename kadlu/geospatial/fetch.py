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
    fnames = ["CA2_4300N06000W.tif",  "CA2_4400N06000W.tif"]
    paths = "/home/oliskir/src/meridian/kadlu/kadlu/tests/assets/tif" + fnames
    return paths