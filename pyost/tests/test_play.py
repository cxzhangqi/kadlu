""" Unit tests for the the 'play' module in the 'pyost' package


    Authors: Mark Thomas and Oliver Kirsebom
    contact: markthomas@dal.ca and oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/pyost
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import pytest
import os

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")


def test_dummy(one):
    assert one == 1
