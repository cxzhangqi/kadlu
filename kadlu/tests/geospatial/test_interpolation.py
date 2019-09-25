""" Unit tests for the the 'geospatial.interpolation' module in the 'kadlu' package

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
from kadlu.geospatial.geospatial import load_data_from_file
from kadlu.geospatial.interpolation import Interpolator2D, Interpolator3D, Uniform2D, Uniform3D, DepthInterpolator3D
import kadlu.geospatial.data_sources.chs as chs
from kadlu.utils import deg2rad, LLtoXY, XYtoLL, LatLon


path_to_assets = os.path.join(os.path.dirname(os.path.dirname(__file__)),"assets")


def test_interpolate_bathymetry_using_latlon_coordinates():
    
    # load bathy data
    path = path_to_assets + '/bornholm.mat'
    bathy, lat, lon = load_data_from_file(path)

    # initialize interpolator
    ip = Interpolator2D(bathy, lat, lon)

    # interpolate at a single point on the grid
    z = ip.eval_ll(lat=lat[0], lon=lon[0]) 
    z = int(z)
    assert z == bathy[0,0]

    # interpolate at a single point between two grid points
    x = (lat[1] + lat[2]) / 2
    z = ip.eval_ll(lat=x, lon=lon[0]) 
    z = float(z)
    zmin = min(bathy[1,0], bathy[2,0])
    zmax = max(bathy[1,0], bathy[2,0])
    assert z >= zmin
    assert z <= zmax

    # interpolate at two points
    x1 = (lat[1] + lat[2])/2
    x2 = (lat[2] + lat[3])/2
    z = ip.eval_ll(lat=[x1,x2], lon=lon[0])
    zmin = min(bathy[1,0], bathy[2,0])
    zmax = max(bathy[1,0], bathy[2,0])
    assert z[0] >= zmin
    assert z[0] <= zmax
    zmin = min(bathy[2,0], bathy[3,0])
    zmax = max(bathy[2,0], bathy[3,0])
    assert z[1] >= zmin
    assert z[1] <= zmax

    # interpolate with grid = True/False
    x1 = (lat[1] + lat[2])/2
    x2 = (lat[2] + lat[3])/2
    y1 = (lon[1] + lon[2])/2
    y2 = (lon[2] + lon[3])/2
    z = ip.eval_ll(lat=[x1,x2], lon=[y1,y2], grid=False)
    assert np.ndim(z) == 1
    assert z.shape[0] == 2 
    z = ip.eval_ll(lat=[x1,x2], lon=[y1,y2], grid=True) 
    assert np.ndim(z) == 2
    assert z.shape[0] == 2 
    assert z.shape[1] == 2 


def test_interpolation_grids_are_what_they_should_be():
    # load data and initialize interpolator
    path = path_to_assets + '/bornholm.mat'
    bathy, lat, lon = load_data_from_file(path)
    ip = Interpolator2D(bathy, lat, lon)
    lat_c = 0.5 * (lat[0] + lat[-1])
    lon_c = 0.5 * (lon[0] + lon[-1])
    assert ip.origin.latitude == lat_c
    assert ip.origin.longitude == lon_c


def test_interpolation_tables_agree_on_latlon_grid():
    # load data and initialize interpolator
    path = path_to_assets + '/bornholm.mat'
    bathy, lat, lon = load_data_from_file(path)
    ip = Interpolator2D(bathy, lat, lon)

    # lat fixed
    ilat = int(len(ip.lat_nodes)/2)
    lat = ip.lat_nodes[ilat]
    for lon in ip.lon_nodes: 
        bll = ip.eval_ll(lat=lat, lon=lon)
        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=ip.origin.latitude, lon_ref=ip.origin.longitude)
        bxy = ip.eval_xy(x=x, y=y)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)

    # lon fixed
    ilon = int(len(ip.lon_nodes)/2)
    lon = ip.lon_nodes[ilon]
    for lat in ip.lat_nodes: 
        bll = ip.eval_ll(lat=lat, lon=lon)
        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=ip.origin.latitude, lon_ref=ip.origin.longitude)
        bxy = ip.eval_xy(x=x, y=y)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)


def test_interpolation_tables_agree_anywhere():
    # load data and initialize interpolator
    path = path_to_assets + '/bornholm.mat'
    bathy, lat, lon = load_data_from_file(path)
    ip = Interpolator2D(bathy, lat, lon)

    # --- at origo ---
    lat_c = ip.origin.latitude
    lon_c = ip.origin.longitude
    z_ll = ip.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = ip.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1e-3) or z_xy == pytest.approx(z_ll, abs=0.1)

    # --- 0.1 degrees north of origo ---
    lat = lat_c + 0.1
    lon = lon_c
    x,y = LLtoXY(lat=lat, lon=lon, lat_ref=lat_c, lon_ref=lon_c)
    z_ll = ip.eval_ll(lat=lat, lon=lon)
    z_ll = float(z_ll)
    z_xy = ip.eval_xy(x=x, y=y) 
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1e-3) or z_xy == pytest.approx(z_ll, abs=0.1)    

    # --- 0.08 degrees south of origo ---
    lat = lat_c - 0.08
    lon = lon_c
    x,y = LLtoXY(lat=lat, lon=lon, lat_ref=lat_c, lon_ref=lon_c)
    z_ll = ip.eval_ll(lat=lat, lon=lon)
    z_ll = float(z_ll)
    z_xy = ip.eval_xy(x=x, y=y) 
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1e-3) or z_xy == pytest.approx(z_ll, abs=0.1)   

    # --- at shifted origo ---
    bathy, lat, lon = load_data_from_file(path)
    ip = Interpolator2D(bathy, lat, lon, origin=LatLon(55.30,15.10))
    lat_c = ip.origin.latitude
    lon_c = ip.origin.longitude
    z_ll = ip.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = ip.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1e-3) or z_xy == pytest.approx(z_ll, abs=0.1)


def test_interpolation_tables_agree_on_ll_grid_for_dbarclays_data():
    # load data and initialize interpolator
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    bathy, lat, lon = load_data_from_file(path, lat_name='latgrat', lon_name='longrat', val_name='mat', lon_axis=0)
    ip = Interpolator2D(bathy, lat, lon)

    # lat fixed
    ilat = int(len(ip.lat_nodes)/2)
    lat = ip.lat_nodes[ilat]
    for lon in ip.lon_nodes: 
        bll = ip.eval_ll(lat=lat, lon=lon)
        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=ip.origin.latitude, lon_ref=ip.origin.longitude)
        bxy = ip.eval_xy(x=x, y=y)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)

    # lon fixed
    ilon = int(len(ip.lon_nodes)/2)
    lon = ip.lon_nodes[ilon]
    for lat in ip.lat_nodes: 
        bll = ip.eval_ll(lat=lat, lon=lon)
        x, y = LLtoXY(lat=lat, lon=lon, lat_ref=ip.origin.latitude, lon_ref=ip.origin.longitude)
        bxy = ip.eval_xy(x=x, y=y)
        assert bxy == pytest.approx(bll, rel=1e-3) or bxy == pytest.approx(bll, abs=0.1)


def test_interpolation_tables_agree_anywhere_for_dbarclays_data():
    # load data and initialize interpolator
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    bathy, lat, lon = load_data_from_file(path, lat_name='latgrat', lon_name='longrat', val_name='mat', lon_axis=0)
    ip = Interpolator2D(bathy, lat, lon)

    # --- at origo ---
    lat_c = ip.origin.latitude
    lon_c = ip.origin.longitude
    z_ll = ip.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = ip.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1E-3) or z_ll == pytest.approx(z_xy, abs=0.1)

    # --- at shifted origo ---
    bathy, lat, lon = load_data_from_file(path, lat_name='latgrat', lon_name='longrat', val_name='mat', lon_axis=0)
    ip = Interpolator2D(bathy, lat, lon, origin=LatLon(9.,140.))
    lat_c = ip.origin.latitude
    lon_c = ip.origin.longitude
    z_ll = ip.eval_ll(lat=lat_c, lon=lon_c) # interpolate using lat-lon
    z_ll = float(z_ll)
    z_xy = ip.eval_xy(x=0, y=0) # interpolate using x-y
    z_xy = float(z_xy)
    assert z_ll == pytest.approx(z_xy, rel=1E-3) or z_ll == pytest.approx(z_xy, abs=0.1)


def test_mariana_trench_is_in_correct_location():
    # load data and initialize interpolator
    path = path_to_assets + '/BathyData_Mariana_500kmx500km.mat'
    bathy, lat, lon = load_data_from_file(path, lat_name='latgrat', lon_name='longrat', val_name='mat', lon_axis=0)
    ip = Interpolator2D(bathy, lat, lon)
    d = ip.eval_ll(lat=11.3733, lon=142.5917)
    assert d < -10770
    d = ip.eval_ll(lat=12.0, lon=142.4)
    assert d > -3000
    d = ip.eval_ll(lat=11.4, lon=143.1)
    assert d < -9000


def test_can_interpolate_multiple_points_in_ll():
    # load data and initialize interpolator
    path = path_to_assets + '/bornholm.mat'
    bathy, lat, lon = load_data_from_file(path)
    ip = Interpolator2D(bathy, lat, lon)
    # coordinates of origin
    lat_c = ip.origin.latitude
    lon_c = ip.origin.longitude
    # --- 4 latitudes ---
    lats = [lat_c, lat_c+0.1, lat_c-0.2, lat_c+0.03]
    # --- 4 longitudes --- 
    lons = [lon_c, lon_c+0.15, lon_c-0.08, lon_c-0.12]
    # interpolate
    depths = ip.eval_ll(lat=lats, lon=lons)
    zi = list()
    for lat, lon in zip(lats, lons):
        zi.append(ip.eval_ll(lat=lat, lon=lon))
    for z,d in zip(zi, depths):
        assert z == pytest.approx(d, rel=1e-3)


def test_can_interpolate_multiple_points_in_xx():
    # load data and initialize interpolator
    path = path_to_assets + '/bornholm.mat'
    bathy, lat, lon = load_data_from_file(path)
    ip = Interpolator2D(bathy, lat, lon)
    # --- 4 x coordinates ---
    xs = [0, 1000, -2000, 300]
    # --- 4 y coordinates --- 
    ys = [0, 1500, 800, -120]
    # interpolate
    depths = ip.eval_xy(x=xs, y=ys)
    zi = list()
    for x, y in zip(xs, ys):
        zi.append(ip.eval_xy(x=x, y=y))
    for z,d in zip(zi, depths):
        assert z == pytest.approx(d, rel=1e-3)


def test_can_interpolate_regular_grid():
    # create fake data
    lat = np.array([44, 45, 46, 47, 48])
    lon = np.array([60, 61, 62, 63])
    bathy = np.random.rand(len(lat), len(lon))
    # initialize interpolator
    ip = Interpolator2D(bathy, lat, lon)
    # check value at grid point
    b = ip.eval_ll(lat=45, lon=62)
    assert b == pytest.approx(bathy[1,2], abs=1E-9)


def test_can_interpolate_irregular_grid():
    # create fake data
    lat = np.array([0.0, 1.0, 1.5, 2.1, 3.0])
    lon = np.array([0.0, 2.0, 0.2, 0.7, 1.2])
    bathy = np.array([-90.0, -200.0, -140.0, -44.0, -301.0])
    # initialize interpolator
    ip = Interpolator2D(bathy, lat, lon)
    # --- 4 latitudes ---
    lats = [0.01, 1.0, 0.5, 2.1]
    # --- 4 longitudes --- 
    lons = [0.01, 2.0, 1.0, 0.71]
    # interpolate all at once
    depths = ip.eval_ll(lat=lats, lon=lons)
    assert depths[1] == pytest.approx(-200, abs=1E-9)
    assert depths[2] < -90 and depths[2] > -200
    # interpolate one at a time
    zi = list()
    for lat, lon in zip(lats, lons):
        zi.append(ip.eval_ll(lat=lat, lon=lon))
    # check that the all-at-once and one-at-a-time 
    # approaches give the same result
    for z,d in zip(zi, depths):
        assert z == pytest.approx(d, rel=1e-3)


def test_can_interpolate_irregular_grid_by_mapping_to_regular_grid():
    # create fake data
    lat = np.array([0.0, 1.0, 1.5, 2.1, 3.0])
    lon = np.array([0.0, 2.0, 0.2, 0.7, 1.2])
    bathy = np.array([-90.0, -200.0, -140.0, -44.0, -501.0])
    # regular grid that data will be mapped to
    lat_reg = np.array([-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
    lon_reg = np.array([-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    # initialize interpolator
    ip = Interpolator2D(bathy, lat, lon, method_irreg='regular', lats_reg=lat_reg, lons_reg=lon_reg)
    # --- 4 latitudes ---
    lats = [0.01, 1.0, 0.5, 2.1]
    # --- 4 longitudes --- 
    lons = [0.01, 2.0, 1.0, 0.71]
    # interpolate all at once
    depths = ip.eval_ll(lat=lats, lon=lons)
    assert depths[1] == pytest.approx(-200, abs=1E-9)
    assert depths[2] < -90 and depths[2] > -200


def test_can_interpolate_geotiff_data():
    # load data and initialize interpolator
    path = path_to_assets + '/tif/CA2_4300N06000W.tif'
    bathy, lat, lon = chs.load_chs_file(path)
    ip = Interpolator2D(bathy, lat, lon, method_irreg='nearest')
    # --- 4 latitudes ---
    lats = [43.3, 43.2, 43.7, 43.5]
    # --- 4 longitudes --- 
    lons = [-59.6, -59.8, -59.2, -59.3]
    # interpolate
    depths = ip.eval_ll(lat=lats, lon=lons)
    zi = list()
    for lat, lon in zip(lats, lons):
        zi.append(ip.eval_ll(lat=lat, lon=lon))
    for z,d in zip(zi, depths):
        assert z == pytest.approx(d, rel=1e-3)
    # interpolate on grid
    depths_grid = ip.eval_ll(lat=lats, lon=lons, grid=True)
    assert depths_grid.shape[0] == 4
    assert depths_grid.shape[1] == 4
    for i in range(4):
        assert depths_grid[i,i] == depths[i]


def test_interpolate_uniform_3d_data():
    N = 10
    np.random.seed(1)
    # create fake data
    val = np.ones(shape=(N,N,N))
    lat = np.arange(N)
    lon = np.arange(N)
    depth = np.arange(N)
    # initialize interpolator
    ip = Interpolator3D(val, lat, lon, depth)
    # check interpolation at a few random points
    lats = np.random.rand(3) * (N - 1)
    lons = np.random.rand(3) * (N - 1)
    depths = np.random.rand(3) * (N - 1)
    vi = ip.eval_ll(lat=lats, lon=lons, z=depths)
    for v in vi:
        assert v == pytest.approx(1, abs=1E-9)
    # check interpolation on a grid
    lats = np.random.rand(3) * (N - 1)
    lons = np.random.rand(4) * (N - 1)
    depths = np.random.rand(5) * (N - 1)
    vi = ip.eval_ll(lat=lats, lon=lons, z=depths, grid=True)
    assert vi.shape[0] == 3
    assert vi.shape[1] == 4
    assert vi.shape[2] == 5
    assert np.all(np.abs(vi - 1.0) < 1E-9)


def test_interpolate_3d_data_with_constant_slope():
    N = 10
    np.random.seed(1)
    # create fake data
    val = np.ones(shape=(N,N,N))
    for k in range(N):
        val[:,:,k] = k * val[:,:,k]        
    lat = np.arange(N)
    lon = np.arange(N)
    depth = np.arange(N)
    # initialize interpolator
    ip = Interpolator3D(val, lat, lon, depth)
    # check interpolation
    lats = np.array([4, 4, 4])
    lons = np.array([4, 4, 4])
    depths = np.array([4, 4.5, 5])
    vi = ip.eval_ll(lat=lats, lon=lons, z=depths)
    assert vi[0] == pytest.approx(4, abs=1E-9)
    assert vi[1] == pytest.approx(4.5, abs=1E-9)
    assert vi[2] == pytest.approx(5, abs=1E-9)
    # check interpolation on grid
    lats = np.array([2, 3, 4])
    lons = np.array([1, 2, 3])
    depths = np.array([4, 4.5, 5])
    vi = ip.eval_ll(lat=lats, lon=lons, z=depths, grid=True)
    assert np.all(np.abs(vi[:,:,0] - 4) < 1E-9)
    assert np.all(np.abs(vi[:,:,1] - 4.5) < 1E-9)
    assert np.all(np.abs(vi[:,:,2] - 5) < 1E-9)

def test_interpolate_3d_data_using_xy_coordinates():
    N = 11
    np.random.seed(1)
    # create fake data
    val = np.ones(shape=(N,N,N))
    for k in range(N):
        val[:,:,k] = k * val[:,:,k]        
    lat = np.arange(N)
    lon = np.arange(N)
    depth = np.arange(N)
    # initialize interpolator
    ip = Interpolator3D(val, lat, lon, depth)
    # check interpolation
    x = np.array([0, 100, 200])
    y = np.array([0, 100, 200])
    depths = np.array([4, 4.5, 5])
    vi = ip.eval_xy(x=x, y=y, z=depths)
    assert vi[0] == pytest.approx(4, abs=1E-9)
    assert vi[1] == pytest.approx(4.5, abs=1E-9)
    assert vi[2] == pytest.approx(5, abs=1E-9)


def test_interpolate_3d_outside_grid():
    N = 10
    np.random.seed(1)
    # create fake data
    val = np.ones(shape=(N,N,N))
    lat = np.arange(N)
    lon = np.arange(N)
    depth = np.arange(N)
    # initialize interpolator
    ip = Interpolator3D(val, lat, lon, depth)
    # check interpolation outside grid
    lats = 20
    lons = 5
    depths = 5
    vi = ip.eval_ll(lat=lats, lon=lons, z=depths)
    assert vi == 1


def test_interpolate_uniform_2d():
    ip = Uniform2D(17)
    v = ip.eval_ll(lat=5, lon=2000)
    assert v == 17
    v = ip.eval_ll(lat=[5, 12, 13], lon=[2000, 0, 1])
    assert np.all(v == 17)
    assert v.shape[0] == 3
    v = ip.eval_ll(lat=[5, 12, 13], lon=[2000, 0], grid=True)
    assert np.all(v == 17)
    assert v.shape[0] == 3
    assert v.shape[1] == 2


def test_interpolate_uniform_3d():
    ip = Uniform3D(17)
    v = ip.eval_ll(lat=5, lon=2000, z=-10)
    assert v == 17
    v = ip.eval_ll(lat=[5, 12, 13], lon=[2000, 0, 1], z=[0, 2, -3])
    assert np.all(v == 17)
    assert v.shape[0] == 3
    v = ip.eval_ll(lat=[5, 12, 13], lon=[2000, 0], z=[0, 2, -3], grid=True)
    assert np.all(v == 17)
    assert v.shape[0] == 3
    assert v.shape[1] == 2
    assert v.shape[2] == 3


def test_interpolate_depth_3d():
    ip = DepthInterpolator3D(values=[0,1,4,9], depths=[0,1,2,3], method='quadratic')
    # inside range
    v = ip.eval_ll(lat=5, lon=2000, z=1.5)
    assert v == 1.5*1.5
    # outside range
    v = ip.eval_ll(lat=5, lon=2000, z=3.5)
    assert v == 3.5*3.5