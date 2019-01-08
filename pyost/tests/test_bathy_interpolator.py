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
    # interpolate at single grid point
    z = interp.eval_ll(lat=lat[0], lon=lon[0]) 
    z = int(z)
    assert z == bathy[0,0]
    # interpolate at single point between two grid points
    x = (lat[1]+lat[2])/2
    z = interp.eval_ll(lat=x, lon=lon[0]) 
    z = float(z)
    zmin = min(bathy[1,0], bathy[2,0])
    zmax = max(bathy[1,0], bathy[2,0])
    assert z >= zmin
    assert z <= zmax
    # interpolate at two points
    x1 = (lat[1]+lat[2])/2
    x2 = (lat[2]+lat[3])/2
    z = interp.eval_ll(lat=[x1,x2], lon=lon[0]) 
    zmin = min(bathy[1,0], bathy[2,0])
    zmax = max(bathy[1,0], bathy[2,0])
    assert z[0] >= zmin
    assert z[0] <= zmax
    zmin = min(bathy[2,0], bathy[3,0])
    zmax = max(bathy[2,0], bathy[3,0])
    assert z[1] >= zmin
    assert z[1] <= zmax
    # interpolate with grid = True/False
    x1 = (lat[1]+lat[2])/2
    x2 = (lat[2]+lat[3])/2
    y1 = (lon[1]+lon[2])/2
    y2 = (lon[2]+lon[3])/2
    z = interp.eval_ll(lat=[x1,x2], lon=[y1,y2], grid=False)
    assert np.ndim(z) == 1
    assert z.shape[0] == 2 
    z = interp.eval_ll(lat=[x1,x2], lon=[y1,y2], grid=True) 
    assert np.ndim(z) == 2
    assert z.shape[0] == 2 
    assert z.shape[1] == 2 
    # interpolate with lat-lon constraints
    interp = BathyInterpolator(bathy_reader=reader, latlon_SW=LatLon(55.10,14.80), latlon_NE=LatLon(55.30,15.10))
    m = len(lat)
    lat = lat[lat >= 55.10]
    m1 = m - len(lat)
    n = len(lon)
    lon = lon[lon >= 14.80]
    n1 = n - len(lon)
    bathy = bathy[m1:,n1:]
    z = interp.eval_ll(lat=lat[0], lon=lon[0]) 
    z = int(z)
    assert z == bathy[0,0]
