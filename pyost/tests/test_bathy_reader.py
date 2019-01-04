""" Unit tests for the the 'bathy_reader' module in the 'pyost' package

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
from bathy_reader import BathyNetCDFReader, LatLon

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_read_bathymetry_from_netcdf_file(one):
    path = path_to_assets + '/bornholm.nc'
    reader = BathyNetCDFReader(path=path, bathy_name='bathy')
    _, _, bathy = reader.read()
    assert np.min(bathy) == -100
    assert np.max(bathy) == 159