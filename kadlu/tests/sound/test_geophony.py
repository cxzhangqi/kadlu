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
import os
import numpy as np
from kadlu.sound.geophony import geophony, transmission_loss, kewley_sl_func, source_level
from kadlu.geospatial.ocean import Ocean
from kadlu.utils import R1_IUGG, deg2rad

current_dir = os.path.dirname(os.path.realpath(__file__))
path_to_assets = os.path.join(os.path.dirname(current_dir),"assets")

def test_kewley_sl_func():
    sl1 = kewley_sl_func(freq=10, wind_uv=0)
    sl2 = kewley_sl_func(freq=40, wind_uv=2.57)
    assert sl1 == sl2
    assert sl2 == 40.0
    sl3 = kewley_sl_func(freq=40, wind_uv=5.14)
    assert sl3 == 44.0
    sl4 = kewley_sl_func(freq=100, wind_uv=5.14)
    assert sl4 == 42.5

def test_source_level():
    ok = {'load_bathymetry': 10000, 'load_wind_uv': 5.14}
    o = Ocean(**ok)
    sl = source_level(freq=10, x=0, y=0, area=1, ocean=o, sl_func=kewley_sl_func)
    assert sl == 44.0
    sl = source_level(freq=100, x=[0,100], y=[0,100], area=[1,2], ocean=o, sl_func=kewley_sl_func)
    assert sl[0] == 42.5
    assert sl[1] == sl[0] + 10*np.log10(2)

def test_geophony_flat_seafloor():
    """ Check that we can execute the geophony method for a 
        flat seafloor and uniform sound speed profile"""
    kwargs = {'load_bathymetry':10000, 'load_wind_uv':1.0, 'ssp':1480, 'angular_bin':90, 'dr':1000, 'dz':1000}
    geo = geophony(freq=100, south=44, north=46, west=-60, east=-58, depth=[100, 2000], xy_res=71, **kwargs)
    spl = geo['spl']
    x = geo['x']
    y = geo['y']
    assert x.shape[0] == 3
    assert y.shape[0] == 5
    assert spl.shape[0] == 3
    assert spl.shape[1] == 5
    assert spl.shape[2] == 2
    assert np.all(np.diff(x) == 71e3)
    assert np.all(np.diff(y) == 71e3)
    # try again, but this time for specific location
    kwargs = {'load_bathymetry':10000, 'load_wind_uv':1.0, 'ssp':1480, 'angular_bin':90, 'dr':1000, 'dz':1000, 'propagation_range':50}
    geo = geophony(freq=100, lat=45, lon=-59, depth=[100, 2000], **kwargs)

def test_geophony_in_canyon(bathy_canyon):
    """ Check that we can execute the geophony method for a 
        canyon-shaped bathymetry and uniform sound speed profile"""
    kwargs = {'load_bathymetry':bathy_canyon, 'load_wind_uv':1.0, 'ssp':1480, 'angular_bin':90, 'dr':1000, 'dz':1000}
    z = [100, 1500, 3000]
    geo = geophony(freq=10, south=43, north=46, west=60, east=62, depth=z, xy_res=71, **kwargs)
    spl = geo['spl']
    x = geo['x']
    y = geo['y']
    assert spl.shape[0] == x.shape[0]
    assert spl.shape[1] == y.shape[0]
    assert spl.shape[2] == len(z)
    assert np.all(np.diff(x) == 71e3)
    assert np.all(np.diff(y) == 71e3)    
    # check that noise is NaN below seafloor and non Nan above
    bathy = np.swapaxes(np.reshape(geo['bathy'], newshape=(y.shape[0], x.shape[0])), 0, 1)
    bathy = bathy[:,:,np.newaxis]
    xyz = np.ones(shape=bathy.shape) * z
    idx = np.nonzero(xyz >= bathy)
    assert np.all(np.isnan(spl[idx]))
    idx = np.nonzero(xyz < bathy)
    assert np.all(~np.isnan(spl[idx]))

def test_transmission_loss_real_world_env():
    """ Check that we can initialize a transmission loss object 
        for a real-world environment and obtain the expected result """
    from datetime import datetime
    bounds = dict(
               start=datetime(2015,1,1), end=datetime(2015,1,2), 
               top=0, bottom=10000
             )    
    src = dict(load_bathymetry='chs', load_temp='hycom', load_salinity='hycom')
    sound_source = {'freq': 200, 'lat': 43.8, 'lon': -59.04, 'source_depth': 12}
    seafloor = {'sound_speed':1700,'density':1.5,'attenuation':0.5}
    transm_loss = transmission_loss(seafloor=seafloor, propagation_range=20, **src, **bounds, **sound_source, dr=100, angular_bin=45, dz=50)
    tl_h, ax_h, tl_v, ax_v = transm_loss.calc(vertical=True)
    assert tl_h.shape == (1,1,8,200)
    assert tl_v.shape == (1,73,201,8)
    answ_h = np.array([[-97.5552,-117.2225,-114.9747,-120.1722,-119.6334,-120.807,-121.7397,-122.591,-123.3596,-124.0323],
                        [-97.5733,-116.2934,-117.8427,-117.4504,-120.2387,-122.1493,-121.2045,-123.6816,-120.3431,-122.9175],
                        [-97.5836,-114.2256,-117.1567,-118.3938,-120.652,-120.261,-122.1713,-120.3823,-125.3299,-121.1278],
                        [-97.5726,-117.5522,-118.0031,-117.7187,-117.7675,-121.0919,-119.8807,-120.6923,-122.5087,-123.5594],
                        [-97.5595,-114.8125,-117.4613,-117.5968,-116.9806,-120.1556,-121.8201,-120.0543,-119.8424,-121.3673],
                        [-97.5614,-110.4824,-115.4014,-119.6987,-119.0094,-117.893,-117.7077,-123.3449,-120.5522,-122.1771],
                        [-97.5671,-110.2498,-114.8029,-117.4865,-123.2917,-120.1903,-121.0737,-121.3798,-123.2269,-120.2656],
                        [-97.5643,-112.0558,-115.5498,-120.1163,-119.0064,-120.1453,-121.1145,-121.961,-122.7209,-123.4159]])
    answ_v = np.array([[-31.9244,-65.4078,-68.0748,-70.133,-71.7806,-72.9215,-73.8679,-74.7315,-75.4971,-76.1683,-76.7786],
                        [-53.4189,-147.1575,-145.934,-151.2488,-152.2335,-151.0661,-151.3919,-150.7127,-151.2096,-152.0068,-152.6981],
                        [-59.4951,-164.7136,-162.7216,-168.1186,-169.0975,-168.0223,-168.3488,-167.6671,-168.1616,-168.9581,-169.6493],
                        [-63.699,-175.0372,-172.9661,-178.374,-179.3523,-178.2953,-178.622,-177.9398,-178.4339,-179.2303,-179.9214],
                        [-67.4673,-182.6341,-180.5805,-185.9919,-186.9699,-185.9195,-186.2463,-185.5639,-186.0578,-186.8542,-187.5453],
                        [-71.5852,-189.1836,-187.1677,-192.5806,-193.5562,-192.5106,-192.8377,-192.156,-192.6501,-193.4464,-194.1375],
                        [-77.339,-195.4623,-193.3593,-200.5519,-201.3583,-200.5171,-201.08,-201.0214,-201.6483,-202.3987,-203.0568],
                        [-93.0928,-208.2705,-209.4284,-217.3889,-218.8682,-218.0186,-218.911,-219.691,-220.4787,-221.1746,-221.7932]])
    np.testing.assert_array_almost_equal(-tl_h[0,0,:,::20], answ_h, decimal=3) 
    np.testing.assert_array_almost_equal(-tl_v[0,1::10,::20,0], answ_v, decimal=3) 

def test_transmission_loss_flat_seafloor():
    """ Check that we can initialize a transmission loss object 
        for a flat seafloor and uniform sound speed profile """
    transm_loss = transmission_loss(freq=100, source_depth=75, propagation_range=0.5, load_bathymetry=2000, ssp=1480, angular_bin=10)
    tl_h, ax_h, tl_v, ax_v = transm_loss.calc(vertical=True)
    answ = np.genfromtxt(os.path.join(path_to_assets, 'lloyd_mirror_f100Hz_SD75m.csv'), delimiter=",")
    assert answ.shape == tl_v[0,:,:,0].shape
    np.testing.assert_array_almost_equal(-tl_v[0,1:,:,0], answ[1:,:], decimal=3) 

#def test_test():
#    from datetime import datetime
#    bounds = dict(
#               south=43.53, north=44.29, west=-59.84, east=-58.48,
#               start=datetime(2015,1,1), end=datetime(2015,1,2), 
#               top=0, bottom=10000
#             )    
#    src = dict(load_bathymetry='chs', load_temp='hycom', load_salinity='hycom')
#    sound_source = {'freq': 200, 'lat': 43.8, 'lon': -59.04, 'source_depth': 12}
#    o = Ocean(**src, **bounds)
#    seafloor = {'sound_speed':1700,'density':1.5,'attenuation':0.5}
#    transm_loss = transmission_loss(seafloor=seafloor, propagation_range=20, **src, **bounds, **sound_source, ssp=1480, dz=50)
#    transm_loss.calc(vertical=False)
