""" Unit tests for the the 'bathy_interpolator' module in the 'kadlu' package

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
from kadlu.bathy_reader import BathyReader, LatLon
from kadlu.bathy_interpolator import BathyInterpolator
from kadlu.utils import R1_IUGG, deg2rad, XYtoLL, LLtoXY

# Degree to radian conversion factor
deg2rad = np.pi / 180.
path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

def test_can_interpolate_latlon():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(input=path, bathy_name='bathy')
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
    reader = BathyReader(input=path, bathy_name='bathy')
    interp = BathyInterpolator(bathy_reader=reader)
    lat, lon, _ = reader.read()
    lat_c = 0.5 * (lat[0] + lat[-1])
    lon_c = 0.5 * (lon[0] + lon[-1])
    # origin is at center
    assert interp.origin.latitude == lat_c
    assert interp.origin.longitude == lon_c

def test_interpolation_tables_agree_on_ll_grid():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(input=path, bathy_name='bathy')
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
    reader = BathyReader(input=path, bathy_name='bathy')
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

def test_interpolation_tables_agree_on_ll_grid_for_dbarclays_data():
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    reader = BathyReader(input=path, lat_name='latgrat', lon_name='longrat', bathy_name='mat', lon_axis=0)
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

def test_interpolation_tables_agree_anywhere_for_dbarclays_data():
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    reader = BathyReader(input=path, lat_name='latgrat', lon_name='longrat', bathy_name='mat', lon_axis=0)
    interp = BathyInterpolator(bathy_reader=reader)
    # --- at origo ---
    lat_c = interp.origin.latitude
    lon_c = interp.origin.longitude
    z_ll = interp.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = interp.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1E-3) or z_ll == pytest.approx(z_xy, abs=0.1)
    # --- at shifted origo ---
    interp = BathyInterpolator(bathy_reader=reader, origin=LatLon(9.,140.))
    lat_c = interp.origin.latitude
    lon_c = interp.origin.longitude
    z_ll = interp.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = interp.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1E-3) or z_ll == pytest.approx(z_xy, abs=0.1)

def test_mariana_trench_is_in_correct_location():
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    reader = BathyReader(input=path, lat_name='latgrat', lon_name='longrat', bathy_name='mat', lon_axis=0)
    interp = BathyInterpolator(bathy_reader=reader)
    d = interp.eval_ll(lat=11.3733, lon=142.5917)
    assert d < -10770
    d = interp.eval_ll(lat=12.0, lon=142.4)
    assert d > -3000
    d = interp.eval_ll(lat=11.4, lon=143.1)
    assert d < -9000

def test_can_interpolate_multiple_points_in_ll():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(input=path, bathy_name='bathy')
    interp = BathyInterpolator(bathy_reader=reader)
    lat_c = interp.origin.latitude
    lon_c = interp.origin.longitude
    # --- 4 latitudes ---
    lats = [lat_c, lat_c+0.1, lat_c-0.2, lat_c+0.03]
    # --- 4 longitudes --- 
    lons = [lon_c, lon_c+0.15, lon_c-0.08, lon_c-0.12]
    # interpolate
    depths = interp.eval_ll(lat=lats, lon=lons)
    zi = list()
    for lat, lon in zip(lats, lons):
        zi.append(interp.eval_ll(lat=lat, lon=lon))
    for z,d in zip(zi, depths):
        assert z == pytest.approx(d, rel=1e-3)

def test_can_interpolate_multiple_points_in_xx():
    path = path_to_assets + '/bornholm.mat'
    reader = BathyReader(input=path, bathy_name='bathy')
    interp = BathyInterpolator(bathy_reader=reader)
    # --- 4 x coordinates ---
    xs = [0, 1000, -2000, 300]
    # --- 4 y coordinates --- 
    ys = [0, 1500, 800, -120]
    # interpolate
    depths = interp.eval_xy(x=xs, y=ys)
    zi = list()
    for x, y in zip(xs, ys):
        zi.append(interp.eval_xy(x=x, y=y))
    for z,d in zip(zi, depths):
        assert z == pytest.approx(d, rel=1e-3)

def test_can_interpolate_unstructured_grid():
    class Reader():
        def __init__(self, latlon_SW=None, latlon_NE=None):
            _=None
        def read(self, latlon_SW=None, latlon_NE=None):
            lats = np.array([0.0, 1.0, 1.5, 2.1, 3.0])
            lons = np.array([0.0, 2.0, 0.2, 0.7, 1.2])
            depths = np.array([-90.0, -200.0, -140.0, -44.0, -301.0])
            return lats, lons, depths

    reader = Reader()
    interp = BathyInterpolator(bathy_reader=reader)
    # --- 4 latitudes ---
    lats = [0.01, 0.1, 0.2, 2.1]
    # --- 4 longitudes --- 
    lons = [0.01, 0.15, 0.08, 0.71]
    # interpolate
    depths = interp.eval_ll(lat=lats, lon=lons)
    zi = list()
    for lat, lon in zip(lats, lons):
        zi.append(interp.eval_ll(lat=lat, lon=lon))
    for z,d in zip(zi, depths):
        assert z == pytest.approx(d, rel=1e-3)

def test_can_interpolate_geotiff_data():
    path = path_to_assets + '/tif/CA2_4300N06000W.tif'
    reader = BathyReader(path)
    interp = BathyInterpolator(bathy_reader=reader)
    # --- 4 latitudes ---
    lats = [43.3, 43.2, 43.7, 43.5]
    # --- 4 longitudes --- 
    lons = [-59.6, -59.8, -59.2, -59.3]
    # interpolate
    depths = interp.eval_ll(lat=lats, lon=lons)
    zi = list()
    for lat, lon in zip(lats, lons):
        zi.append(interp.eval_ll(lat=lat, lon=lon))
    for z,d in zip(zi, depths):
        assert z == pytest.approx(d, rel=1e-3)
    # interpolate on grid
    depths_grid = interp.eval_ll(lat=lats, lon=lons, grid=True)
    assert depths_grid.shape[0] == 4
    assert depths_grid.shape[1] == 4
    for i in range(4):
        assert depths_grid[i,i] == depths[i]

def test_can_interpolate_geotiff_data_in_xy():
    path = path_to_assets + '/tif/CA2_4300N06000W.tif'
    reader = BathyReader(path)
    interp = BathyInterpolator(bathy_reader=reader, method='nearest')
    depth = interp.eval_xy(x=-20e3, y=20e3)
    assert depth == -82.
    depth = interp.eval_ll(lat=43.69, lon=-59.75)
    assert depth == -82.
    depth = interp.eval_xy(x=-30e3, y=40e3)
    assert depth == -43.
    depth = interp.eval_xy(x=-20e3, y=40e3)
    assert depth == -38.
    depth = interp.eval_xy(x=-30e3, y=20e3)
    assert depth == -80.
    a = [-20e3, -30e3]
    b = [20e3, 40e3]
    depth = interp.eval_xy(x=a, y=b)
    assert depth[0] == -82.
    assert depth[1] == -43.
    depth = interp.eval_xy(x=a, y=b, grid=True)
    assert depth[0][0] == -82.
    assert depth[1][1] == -43.
    assert depth[0][1] == -38.
    assert depth[1][0] == -80.

def test_can_interpolate_geotiff_data_in_ll():
    path = path_to_assets + '/tif/CA2_4300N06000W.tif'
    reader = BathyReader(path)
    interp = BathyInterpolator(bathy_reader=reader, method='nearest')
    depth = interp.eval_ll(lat=43.6, lon=-59.8)
    assert depth == -188.
    depth = interp.eval_ll(lat=43.2, lon=-59.8)
    assert depth == pytest.approx(-2011.7, abs=0.1)
    a = [43.6, 43.2]
    b = [-59.8, -59.8]
    depth = interp.eval_ll(lat=a, lon=b)
    assert depth[0] == -188.
    assert depth[1] == pytest.approx(-2011.7, abs=0.1)
    depth = interp.eval_ll(lat=a, lon=b, grid=True)
    assert depth[0][0] == -188.
    assert depth[1][1] == pytest.approx(-2011.7, abs=0.1)
    assert depth[0][1] == -188.
    assert depth[1][0] == pytest.approx(-2011.7, abs=0.1)

def test_can_interpolate_gradient_ll():
    class Reader():
        def __init__(self, latlon_SW=None, latlon_NE=None):
            _=None
        def read(self, latlon_SW=None, latlon_NE=None):
            lats = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
            lons = np.array([0.0, 2.0, 4.0, 6.0, 8.0])
            depths = np.array([[-100, -100, -100, -100, -100],\
                            [-200, -200, -200, -200, -200],\
                            [-300, -300, -300, -300, -300],\
                            [-400, -400, -400, -400, -400],\
                            [-500, -500, -500, -500, -500]])   
            return lats, lons, depths

    reader = Reader()
    interp = BathyInterpolator(bathy_reader=reader)
    # interpolate as usual
    depth = interp.eval_ll(lat=1.0, lon=1.0)
    assert depth == pytest.approx(-200, abs=0.001)
    # interpolate lat-gradient
    dzdy = interp.eval_ll(lat=1.0, lon=1.0, lat_deriv_order=1)
    assert dzdy == pytest.approx(-100./(1.0*np.pi/180), abs=0.001)
    # interpolate lat+lon gradient
    dzdydx = interp.eval_ll(lat=1.0, lon=1.0, lat_deriv_order=1, lon_deriv_order=1)
    assert dzdydx == pytest.approx(0, abs=0.001)

def test_can_interpolate_gradient_xy():
    class Reader():
        def __init__(self, latlon_SW=None, latlon_NE=None):
            _=None
        def read(self, latlon_SW=None, latlon_NE=None):
            lats = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
            lons = np.array([0.0, 2.0, 4.0, 6.0, 8.0])
            depths = np.array([[-100, -100, -100, -100, -100],\
                            [-200, -200, -200, -200, -200],\
                            [-300, -300, -300, -300, -300],\
                            [-400, -400, -400, -400, -400],\
                            [-500, -500, -500, -500, -500]])   
            return lats, lons, depths

    reader = Reader()
    interp = BathyInterpolator(bathy_reader=reader, origin=LatLon(0,0))
    # convert lat-lon to x-y
    x,y = LLtoXY(lat=1.0, lon=1.0, lat_ref=0, lon_ref=0)
    _,y1 = LLtoXY(lat=2.0, lon=1.0, lat_ref=0, lon_ref=0)
    dy = y1 - y
    # interpolate as usual
    depth = interp.eval_xy(x=x, y=y)
    assert depth == pytest.approx(-200, abs=0.001)
    # interpolate y-gradient
    dzdy = interp.eval_xy(x=x, y=y, y_deriv_order=1)
    assert dzdy == pytest.approx(-100./dy, abs=0.001)
    # interpolate y+x gradient
    dzdydx = interp.eval_xy(x=x, y=y, y_deriv_order=1, x_deriv_order=1)
    assert dzdydx == pytest.approx(0, abs=0.001)
