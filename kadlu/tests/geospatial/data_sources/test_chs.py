""" Unit tests for the the 'geospatial.bathy_chs' module in the 'kadlu' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import pytest
import numpy as np
import os
from kadlu.geospatial.data_sources import chs
from kadlu.geospatial.data_sources.chs import Chs
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
from datetime import datetime, timedelta


path_to_assets = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),"assets")

# the script will now force-fetch regardless of whether local files exist, so this is no longer needed
# local file checking is now performed inside the load function instead of the fetch function
"""
def unfetch():
    for f in os.listdir(storage_cfg()):
        if "CA2_" in f:
            print(f"Removing {f}")
            os.remove(f"{storage_cfg()}{f}")
    return
"""


def test_fetch_returns_expected_number_of_files():
    paths = Chs().fetch_bathymetry(south=42.1, north=44.6, west=-60.9, east=-59.1)
    assert len(paths) == 6
    assert os.path.basename(paths[0]) == "CA2_4200N06000W.tif"

def test_filename():
    fname = chs.filename(south=50, west=-62)
    assert fname == "CA2_5000N06200W.tif"

def test_parse_sw_corner():
    s,w = chs.parse_sw_corner(f"{storage_cfg()}CA2_5000N06200W.tif")
    assert s == 50
    assert w == -62

# matt_s 2019-11-12
# this function has been moved to inside the load_chs_file function 
# to explicitly force loose coupling. to test latlon, write tests for 
# load_chs_file function 
"""
def test_latlon():
    Chs().fetch_bathymetry(south=49, north=51, west=-63, east=-61)
    lats, lons = chs.latlon(f"{storage_cfg()}CA2_5000N06200W.tif")
    assert lats.shape[0] == 1001
    assert lons.shape[0] == 1001
    assert lats[0] == 50
    assert lats[1] == 50.001
    assert lons[0] == -62
    assert lons[1] == -61.999
    lats, lons = chs.latlon(f"{storage_cfg()}CA2_5000N06200W.tif")

    # removing these tests since latlon function now dynamically generates array size

    #assert lats.shape[0] == 501
    #assert lats[0] == 50
    #assert lats[1] == 50.001
"""

def test_load_single_chs_file():
    bathy, lats, lons = Chs().load_bathymetry(south=43, west=-60, north=44, east=-59)
    #assert np.ma.min(bathy) == pytest.approx(-3257.100, abs=0.001)
    #assert np.ma.max(bathy) == pytest.approx(1.645, abs=0.001)
    assert bathy.shape[0] == lats.shape[0]
    assert bathy.shape[0] == lons.shape[0]

def test_fetch_correct_number_of_files():
    south = 44.006
    north = 45.994
    west= -59.994
    east= -58.006
    filenames = Chs().fetch_bathymetry(south, north, west, east)
    assert(len(filenames) == 4)
    filenames_str = '\n'.join(filenames)
    """
    assert("4400N05800W.tif" in filenames_str)
    assert("4400N05900W.tif" in filenames_str)
    assert("4500N05800W.tif" in filenames_str)
    assert("4500N05900W.tif" in filenames_str)
    """

