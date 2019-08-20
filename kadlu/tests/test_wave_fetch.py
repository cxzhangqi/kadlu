import pytest
import os
from os import path
from os.path import dirname
from datetime import datetime, timedelta
from kadlu.wave_fetch import WaveFetch

# matt_s 2019-08
# note: still need to find user friendly way to expose the 
# waveSources and regions interfaces to the user

from kadlu.wave_fetch import waveSources
from kadlu.wave_fetch import regions
print(f"Valid wavesouces are:\n{waveSources}")
print(f"Valid regions are:\n{regions}")

def test_ERA5_fetch():
    wf = WaveFetch()
    wf.fetchERA5(wavevar=waveSources['ERA5']['swh'], time=datetime(2018, 1, 1, 0, 0, 0, 0))

def test_WWIII_fetch():
    wf = WaveFetch()
    wf.fetchWWIII(wavevar=waveSources['WWIII']['swh'], time=datetime(2017, 2, 3, 0, 0, 0, 0), region=regions['WWIII']['global'])

def test_RDWPS_fetch():
    wf = WaveFetch()
    wf.fetchRDWPS(wavevar=waveSources['RDWPS']['swh'], time=datetime.now()-timedelta(hours=3), region=regions['RDWPS'][0])

