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
from os import path
from os.path import dirname
import numpy as np
from datetime import datetime, timedelta
from enum import Enum
from kadlu.wave_fetch import WaveFetch, WaveSources


def test_ERA5_fetch():
    wfTestERA5 = WaveFetch(fetch_datetimestamp=datetime(2018, 1, 1, 0, 0, 0, 0), wave_source=WaveSources.ECMWF_ERA5)
    gribdata = wfTestERA5.fetchERA5()
    assert(isinstance(gribdata[0][0], np.ndarray))
    (grbsample, samp_title_text) = wfTestERA5.loadERA5()
    #wfTestERA5.plotSampleGrib(grbsample, samp_title_text)

def test_WWIII_fetch():
    wfTestWWIII = WaveFetch(fetch_datetimestamp=datetime(2017, 2, 3, 0, 0, 0, 0), wave_source=WaveSources.NOAA_WWIII)
    gribdata = wfTestWWIII.fetchWWIII()
    assert(isinstance(gribdata[0][0], np.ndarray))
    (grbsampleWWIII, samp_title_text_WWIII) = wfTestWWIII.loadWWIII()
    #wfTestWWIII.plotSampleGrib(grbsampleWWIII, samp_title_text_WWIII)

def test_RDWPS_fetch():
    wfTestRDWPS = WaveFetch(fetch_datetimestamp=datetime.now()-timedelta(hours=3), wave_source=WaveSources.EC_RDWPS)
    gribdata = wfTestRDWPS.fetchRDWPS()
    assert(isinstance(gribdata[0][0], np.ndarray))
    (grbsampleRDWPS, samp_title_textRDWPS) = wfTestRDWPS.loadRDWPS()
    #wfTestRDWPS.plotSampleGrib(grbsampleRDWPS, samp_title_textRDWPS)

def test_DalCoast_fetch():
    wfTestDalCoast = WaveFetch(fetch_datetimestamp=datetime.now(), wave_source=None)
    wfTestDalCoast.fetchDalCoast()
    (grbsampleDalCoast, samp_title_text_DalCoast) = wfTestDalCoast.loadDalCoast()
    #wfTest.plotSampleGrib(grbsampleDalCoast, samp_title_text_DalCoast)



