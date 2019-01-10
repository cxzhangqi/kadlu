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
from pyost.bathy_reader import BathyReader, LatLon
from pyost.bathy_interpolator import BathyInterpolator
from pyost.util import R1_IUGG, deg2rad, XYtoLL, LLtoXY

# Degree to radian conversion factor
deg2rad = np.pi / 180.
path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_interpolate_latlon():
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

def test_interpolation_grids_are_what_they_should_be():
    path = path_to_assets + '/bornholm.mat'
    rebin = 3
    reader = BathyReader(path=path, bathy_name='bathy')
    interp = BathyInterpolator(bathy_reader=reader, rebin_xy=rebin)
    lat, lon, _ = reader.read()
    lat_c = 0.5 * (lat[0] + lat[-1])
    lon_c = 0.5 * (lon[0] + lon[-1])
    # origin is at center
    assert interp.origin.latitude == lat_c
    assert interp.origin.longitude == lon_c
    # x,y have same number of nodes as lon,lat
    assert len(interp.x_nodes) == rebin * len(interp.lon_nodes)
    assert len(interp.y_nodes) == rebin * len(interp.lat_nodes)
    # x and y are symmetric around 0
    assert interp.x_nodes[0] == -interp.x_nodes[-1]
    assert interp.y_nodes[0] == -interp.y_nodes[-1]
    # x and y are regularly spaced
    xdiff = np.diff(interp.x_nodes)
    ydiff = np.diff(interp.y_nodes)
    assert np.all(pytest.approx(xdiff[0] == xdiff, rel=1E-9))
    assert np.all(pytest.approx(ydiff[0] == ydiff, rel=1E-9))

def test_interpolation_tables_agree_on_xy_grid():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(path=path, bathy_name='bathy')
    interp = BathyInterpolator(bathy_reader=reader)
    # x fixed
    ix = int(len(interp.x_nodes)/2)
    x = interp.x_nodes[ix]
    for y in interp.y_nodes: 
        bxy = interp.eval_xy(x=x, y=y)
        la, lo = XYtoLL(x=x, y=y, lat_ref=interp.origin.latitude, lon_ref=interp.origin.longitude)
        bll = interp.eval_ll(lat=la, lon=lo)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)
    # y fixed
    iy = int(len(interp.y_nodes)/2)
    y = interp.y_nodes[iy]
    for x in interp.x_nodes: 
        bxy = interp.eval_xy(x=x, y=y)
        la, lo = XYtoLL(x=x, y=y, lat_ref=interp.origin.latitude, lon_ref=interp.origin.longitude)
        bll = interp.eval_ll(lat=la, lon=lo)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)

def test_interpolation_tables_agree_on_ll_grid():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(path=path, bathy_name='bathy')
    interp = BathyInterpolator(bathy_reader=reader)
    # lat fixed
    ilat = int(len(interp.lat_nodes)/2)
    lat = interp.lat_nodes[ilat]
    for lon in interp.lon_nodes: 
        bll = interp.eval_ll(lat=lat, lon=lon)
        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=interp.origin.latitude, lon_ref=interp.origin.longitude)
        bxy = interp.eval_xy(x=x, y=y)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)
    # lon fixed
    ilon = int(len(interp.lon_nodes)/2)
    lon = interp.lon_nodes[ilon]
    for lat in interp.lat_nodes: 
        bll = interp.eval_ll(lat=lat, lon=lon)
        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=interp.origin.latitude, lon_ref=interp.origin.longitude)
        bxy = interp.eval_xy(x=x, y=y)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)

def test_interpolation_tables_agree_anywhere():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(path=path, bathy_name='bathy')
    interp = BathyInterpolator(bathy_reader=reader)
    # --- at origo ---
    lat_c = interp.origin.latitude
    lon_c = interp.origin.longitude
    z_ll = interp.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = interp.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1e-3) or z_xy == pytest.approx(z_ll, abs=0.1)
    # --- 0.1 degrees north of origo ---
    lat = lat_c + 0.1
    lon = lon_c
    x,y = LLtoXY(lat=lat, lon=lon, lat_ref=lat_c, lon_ref=lon_c)
    z_ll = interp.eval_ll(lat=lat, lon=lon)
    z_ll = float(z_ll)
    z_xy = interp.eval_xy(x=x, y=y) 
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1e-3) or z_xy == pytest.approx(z_ll, abs=0.1)    
    # --- 0.08 degrees south of origo ---
    lat = lat_c - 0.08
    lon = lon_c
    x,y = LLtoXY(lat=lat, lon=lon, lat_ref=lat_c, lon_ref=lon_c)
    z_ll = interp.eval_ll(lat=lat, lon=lon)
    z_ll = float(z_ll)
    z_xy = interp.eval_xy(x=x, y=y) 
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1e-3) or z_xy == pytest.approx(z_ll, abs=0.1)    
    # --- at shifted origo ---
    interp = BathyInterpolator(bathy_reader=reader, origin=LatLon(55.30,15.10))
    lat_c = interp.origin.latitude
    lon_c = interp.origin.longitude
    z_ll = interp.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = interp.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1e-3) or z_xy == pytest.approx(z_ll, abs=0.1)

def test_interpolation_tables_agree_on_xy_grid_for_dbarclays_data():
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    reader = BathyReader(path=path, lat_name='latgrat', lon_name='longrat', bathy_name='mat')
    interp = BathyInterpolator(bathy_reader=reader, rebin_xy=4)
    # x fixed
    ix = int(len(interp.x_nodes)/2)
    x = interp.x_nodes[ix]
    for y in interp.y_nodes: 
        bxy = interp.eval_xy(x=x, y=y)
        la, lo = XYtoLL(x=x, y=y, lat_ref=interp.origin.latitude, lon_ref=interp.origin.longitude)
        bll = interp.eval_ll(lat=la, lon=lo)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)
    # y fixed
    iy = int(len(interp.y_nodes)/2)
    y = interp.y_nodes[iy]
    for x in interp.x_nodes: 
        bxy = interp.eval_xy(x=x, y=y)
        la, lo = XYtoLL(x=x, y=y, lat_ref=interp.origin.latitude, lon_ref=interp.origin.longitude)
        bll = interp.eval_ll(lat=la, lon=lo)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)

def test_interpolation_tables_agree_on_ll_grid_for_dbarclays_data():
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    reader = BathyReader(path=path, lat_name='latgrat', lon_name='longrat', bathy_name='mat')
    interp = BathyInterpolator(bathy_reader=reader, rebin_xy=8)
    # lat fixed
    ilat = int(len(interp.lat_nodes)/2)
    lat = interp.lat_nodes[ilat]
    for lon in interp.lon_nodes: 
        bll = interp.eval_ll(lat=lat, lon=lon)
        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=interp.origin.latitude, lon_ref=interp.origin.longitude)
        bxy = interp.eval_xy(x=x, y=y)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)
    # lon fixed
    ilon = int(len(interp.lon_nodes)/2)
    lon = interp.lon_nodes[ilon]
    for lat in interp.lat_nodes: 
        bll = interp.eval_ll(lat=lat, lon=lon)
        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=interp.origin.latitude, lon_ref=interp.origin.longitude)
        bxy = interp.eval_xy(x=x, y=y)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)

def test_interpolation_tables_agree_anywhere_for_dbarclays_data():
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    reader = BathyReader(path=path, lat_name='latgrat', lon_name='longrat', bathy_name='mat')
    interp = BathyInterpolator(bathy_reader=reader, rebin_xy=4)
    # --- at origo ---
    lat_c = interp.origin.latitude
    lon_c = interp.origin.longitude
    z_ll = interp.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = interp.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1E-3) or z_ll == pytest.approx(z_xy, abs=0.1)
    # --- at shifted origo ---
    interp = BathyInterpolator(bathy_reader=reader, rebin_xy=4, origin=LatLon(9.,140.))
    lat_c = interp.origin.latitude
    lon_c = interp.origin.longitude
    z_ll = interp.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = interp.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1E-3) or z_ll == pytest.approx(z_xy, abs=0.1)