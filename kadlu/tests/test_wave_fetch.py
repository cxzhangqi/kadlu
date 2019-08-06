"""
    Unit tests for 'wave fetch' module in the 'kadlu' package
    
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
from kadlu.wave_fetch import WaveFetch, WaveSources, RDWPSRegion
from kadlu.wave_reader import WaveReader
from mpl_toolkits.basemap import Basemap
import urllib.request

import pygrib



def test_can_download_gribfile():
    fetch_year = "2019"
    fetch_month = "01"

    regionAndWaveName = "ak_4m.dp"
    fetchURLString = 'https://data.nodc.noaa.gov/ncep/nww3/' + fetch_year + '/' + fetch_month + '/gribs/multi_1.' + regionAndWaveName + '.' + fetch_year + fetch_month + '.grb2'

    response = urllib.request.urlopen(fetchURLString)
    data = response.read()
    assert(data[0:4] == b'GRIB')  # make sure file is binary grib data
