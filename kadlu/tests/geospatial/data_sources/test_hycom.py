import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import hycom
from kadlu.geospatial.data_sources.hycom import Hycom
from os.path import isfile

# run pytest with '--capture=no' arg to print class function info 
print(Hycom())

# mahone bay test area:
south =  44.4
north =  44.7
west  = -64.4
east  = -63.8
time  = datetime(2015, 1, 1)


def test_fetch_salinity():
    fnames = Hycom().fetch_salinity(south=south, north=north, west=west, east=east, time=time)
    for f in fnames: assert(isfile(f))

def test_load_salinity():
    salinity, lat, lon = Hycom().load_salinity(south=south, north=north, west=west, east=east, time=time)

def test_fetch_temp():
    fnames = Hycom().fetch_temp(south=south, north=north, west=west, east=east, time=time)
    for f in fnames: assert(isfile(f))

def test_load_temp():
    salinity, lat, lon = Hycom().load_temp(south=south, north=north, west=west, east=east, time=time)

def test_fetch_water_u():
    fnames = Hycom().fetch_water_u(south=south, north=north, west=west, east=east, time=time)
    for f in fnames: assert(isfile(f))

def test_load_water_u():
    salinity, lat, lon = Hycom().load_water_u(south=south, north=north, west=west, east=east, time=time)

def test_fetch_water_v():
    fnames = Hycom().fetch_water_v(south=south, north=north, west=west, east=east, time=time)
    for f in fnames: assert(isfile(f))

def test_load_water_v():
    salinity, lat, lon = Hycom().load_water_u(south=south, north=north, west=west, east=east, time=time)

def test_hycom_dt_2_t():
    time = datetime(2015, 1, 1)
    year_str, daterange = hycom.dt_2_t(time)
    assert(year_str == "2015")
    assert(daterange == (0, 0))
    # more test cases for time conversion here


