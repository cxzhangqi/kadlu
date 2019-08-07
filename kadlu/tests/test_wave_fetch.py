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
from kadlu.wave_fetch import WaveFetch, WaveSources

# currently the default storage location is ~/kadlu/storage
# maybe in the future this can be configured by the user
storage = (os.getenv("HOME") + "/kadlu/storage")


def test_ERA5_fetch():
    wfTestERA5 = WaveFetch(storage, datetime(2018, 1, 1, 0, 0, 0, 0), WaveSources.ECMWF_ERA5)
    wfTestERA5.fetchERA5()
    (grbsample, samp_title_text) = wfTestERA5.loadERA5()
    wfTestERA5.plotSampleGrib(grbsample, samp_title_text)

def test_WWIII_fetch():
    wfTestWWIII = WaveFetch(storage, datetime(2017, 2, 3, 0, 0, 0, 0), WaveSources.NOAA_WWIII)
    wfTestWWIII.fetchWWIII()
    (grbsampleWWIII, samp_title_text_WWIII) = wfTestWWIII.loadWWIII()
    wfTestWWIII.plotSampleGrib(grbsampleWWIII, samp_title_text_WWIII)

def test_RDWPS_fetch():
    wfTestRDWPS = WaveFetch(storage, datetime.now()  - timedelta(hours=3), WaveSources.EC_RDWPS)
    wfTestRDWPS.fetchRDWPS()
    (grbsampleRDWPS, samp_title_textRDWPS) = wfTestRDWPS.loadRDWPS()
    wfTestRDWPS.plotSampleGrib(grbsampleRDWPS, samp_title_textRDWPS)


