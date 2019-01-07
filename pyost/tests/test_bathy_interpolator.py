""" Unit tests for the the 'bathy_interpolator' module in the 'pyost' package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/pyost
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import pytest
import os
import numpy as np
from bathy_reader import BathyReader, LatLon
from bathy_interpolator import BathyInterpolator

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_interpolate_bornholm():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(path=path, bathy_name='bathy')
    # interpolate without lat-lon constraints
    interp = BathyInterpolator(bathy_reader=reader)
    lat, lon, bathy = reader.read()
    z = interp.eval_ll(lat=lat[0], lon=lon[0])
    print(z.shape)
    assert z == bathy[0,0]
    # interpolate with lat-lon constraints
    interp = BathyInterpolator(bathy_reader=reader, latlon_SW=LatLon(45.5,4.5), latlon_NE=LatLon(47.5,6.5))
