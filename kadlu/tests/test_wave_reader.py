"""
    Unit tests for 'wave reader' module in the 'kadlu' package
    
    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
"""

import pytest
import os
import numpy as np
from datetime import datetime, timedelta
from enum import Enum
# TODO:
# this has changed with refactoring to wave_fetch
# test should now use waveSources and region dictionaries
#from kadlu.wave_fetch import WaveFetch, WaveSources, RDWPSRegion
from kadlu.wave_reader import WaveReader
from mpl_toolkits.basemap import Basemap
import urllib.request

import pygrib
import matplotlib.pyplot as plt
import tempfile


path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_read_wave_from_file():
    #storage_location = path_to_assets + "/somedata.file"
    #timestamp = "00000000"  # unknown timestamp format
    #wave_source = WaveSources()
    #reader = WaveFetch(storage_location, fetch_date, wave_source)
    pass


def test_can_read_ERA5():
    pass

def test_can_read_WWIII():
    pass

def test_can_read_RDWPS():
    pass

