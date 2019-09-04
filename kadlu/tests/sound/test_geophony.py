""" Unit tests for the the 'sound.geophony' module in the 'kadlu' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import pytest
import math
import numpy as np
from kadlu.sound.geophony import Geophony
from kadlu.sound.sound_propagation import TLCalculator, Seafloor
from kadlu.geospatial.ocean import Ocean


def test_initialize_geophony():
    s = Seafloor()
    o = Ocean()
    tl = TLCalculator(ocean=o, seafloor=s)
    geo = Geophony(tl_calculator=tl, depth=[-100, -200, -300])

def test_model_geophony():
    s = Seafloor(thickness=2000)
    o = Ocean(bathy=-10000, wave=1.0)
    tl = TLCalculator(ocean=o, seafloor=s, sound_speed=1480, radial_bin=100, radial_range=50e3, angular_bin=45, vertical_bin=100)
    geo = Geophony(tl_calculator=tl, depth=[-100, -2000, -9000])
    geo.model(frequency=1000, south=44, north=45, west=60, east=61)