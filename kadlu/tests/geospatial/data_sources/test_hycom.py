import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import hycom
from kadlu.geospatial.data_sources.hycom import Hycom
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
import os
from os.path import isfile

# run pytest with '--capture=no' arg to print class function info 
print(Hycom())

# mahone bay test area:
south =  44.4
north =  44.7
west  = -64.4
east  = -63.8
start = datetime(2015, 1, 1)
end   = datetime(2015, 1, 1, 2)

test_fetch = True

# remove fetched files to test fetching
def unfetch():
    return 
    for f in os.listdir(storage_cfg()):
        if "hycom" in f:
            print(f"Removing {f}")
            os.remove(f"{storage_cfg()}{f}")


def test_fetch_salinity():
    if not test_fetch: return
    unfetch()
    fnames = Hycom().fetch_salinity(south=south, north=north, west=west, east=east, start=start, end=end)
    for f in fnames: assert(isfile(f))

def test_load_salinity():
    salinity, lat, lon = Hycom().load_salinity(south=south, north=north, west=west, east=east, start=start, end=end)

def test_fetch_temp():
    if not test_fetch: return
    unfetch()
    fnames = Hycom().fetch_temp(south=south, north=north, west=west, east=east, start=start, end=end)
    for f in fnames: assert(isfile(f))

def test_load_temp():
    salinity, lat, lon = Hycom().load_temp(south=south, north=north, west=west, east=east, start=start, end=end)

def test_fetch_water_u():
    if not test_fetch: return
    unfetch()
    fnames = Hycom().fetch_water_u(south=south, north=north, west=west, east=east, start=start, end=end)
    for f in fnames: assert(isfile(f))

def test_load_water_u():
    salinity, lat, lon = Hycom().load_water_u(south=south, north=north, west=west, east=east, start=start, end=end)

def test_fetch_water_v():
    if not test_fetch: return
    unfetch()
    fnames = Hycom().fetch_water_v(south=south, north=north, west=west, east=east, start=start, end=end)
    for f in fnames: assert(isfile(f))

def test_load_water_v():
    salinity, lat, lon = Hycom().load_water_u(south=south, north=north, west=west, east=east, start=start, end=end)

def test_hycom_dt_2_tslice():
    start = datetime(2015, 1, 1)
    end = datetime(2015, 1, 7)
    dateslice = hycom.dt_2_tslice(start, end)
    assert(dateslice[0] == 0)
    end = datetime(2015, 12, 31, 23, 59)
    dateslice = hycom.dt_2_tslice(start, end)
    assert(2859 <= dateslice[1] <= 2860)

    # known issue: hycom time conversion is only accurate within 3 hours

    # more test cases for time conversion here


