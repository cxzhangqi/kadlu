""" Unit tests for the the 'sound.sound_propagation' module in the 'kadlu' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""
import pytest
import os
import numpy as np
from kadlu.sound.sound_propagation import TLCalculator, Seafloor

def test_initialize_seafloor_with_default_args():
    s = Seafloor()
    assert s.c == 1700
    assert s.density == 1.5
    assert s.thickness == 2000
    assert s.loss == 0.5

def test_initialize_seafloor_with_user_specified_args():
    s = Seafloor(c=1555, density=0.8, thickness=1000, loss=0.8)
    assert s.c == 1555
    assert s.density == 0.8
    assert s.thickness == 1000
    assert s.loss == 0.8

def test_attempt_to_compute_nsq_gives_error():
    s = Seafloor()
    with pytest.raises(AssertionError):
        _ = s.nsq()

def test_compute_nsq():
    s = Seafloor()
    s.c0 = 1500
    s.frequency = 20
    n = s.nsq()
    assert np.real(n) == pytest.approx(0.77835066, abs=1E-6)
    assert np.imag(n) == pytest.approx(0.01426442, abs=1E-6)