import pytest
from datetime import datetime, timedelta
from kadlu.geospatial.data_sources import hycom
from kadlu.geospatial.data_sources.hycom import Hycom
from kadlu.geospatial.data_sources.fetch_util import storage_cfg
import os
from os.path import isfile

# run pytest with '--capture=no' arg to print class function info 
print(Hycom())

# gulf st lawrence test area:
south =  46
north =  52
west  = -70
east  = -56
start = datetime(2015, 1, 1)
end   = datetime(2015, 1, 1, 6)

# disable automated testing of fetching to avoid slamming the API with
# requests in the automated development pipeline
test_fetch = False

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
    val, lat, lon, time, depth = Hycom().load_salinity(south=south, north=north, west=west, east=east, start=start, end=end)
    # leaving values as 4d array for now
    pass
    return
    assert(len(val) == len(lat) == len(lon) == len(time))
    assert(sum(lat <= 90) == sum(lat >= -90) == len(lat))
    assert(sum(lon <= 180) == sum(lon >= -180) == len(lon))

def test_fetch_temp():
    if not test_fetch: return
    unfetch()
    fnames = Hycom().fetch_temp(south=south, north=north, west=west, east=east, start=start, end=end)
    for f in fnames: assert(isfile(f))

def test_load_temp():
    val, lat, lon, time, depth = Hycom().load_temp(south=south, north=north, west=west, east=east, start=start, end=end)
    # leaving values as 4d array for now
    pass
    return 
    assert(len(val) == len(lat) == len(lon) == len(time))
    assert(sum(lat <= 90) == sum(lat >= -90) == len(lat))
    assert(sum(lon <= 180) == sum(lon >= -180) == len(lon))

def test_fetch_water_u():
    if not test_fetch: return
    unfetch()
    fnames = Hycom().fetch_water_u(south=south, north=north, west=west, east=east, start=start, end=end)
    for f in fnames: assert(isfile(f))

def test_load_water_u():
    val, lat, lon, time, depth = Hycom().load_water_u(south=south, north=north, west=west, east=east, start=start, end=end)
    # leaving values as 4d array for now
    pass
    return 
    assert(len(val) == len(lat) == len(lon) == len(time))
    assert(sum(lat <= 90) == sum(lat >= -90) == len(lat))
    assert(sum(lon <= 180) == sum(lon >= -180) == len(lon))

def test_fetch_water_v():
    if not test_fetch: return
    unfetch()
    fnames = Hycom().fetch_water_v(south=south, north=north, west=west, east=east, start=start, end=end)
    for f in fnames: assert(isfile(f))

def test_load_water_v():
    val, lat, lon, time, depth = Hycom().load_water_u(south=south, north=north, west=west, east=east, start=start, end=end)
    # leaving values as 4d array for now
    pass
    return 
    assert(len(val) == len(lat) == len(lon) == len(time))
    assert(sum(lat <= 90) == sum(lat >= -90) == len(lat))
    assert(sum(lon <= 180) == sum(lon >= -180) == len(lon))

def test_hycom_dt_2_tslice():
    start = datetime(2015, 1, 1)
    end = datetime(2015, 1, 7)
    dateslice = hycom.dt_2_tslice(start, end, Hycom().time)
    assert(dateslice[0] == 0)
    end = datetime(2015, 12, 31, 23, 59)
    dateslice = hycom.dt_2_tslice(start, end, Hycom().time)
    assert(2859 <= dateslice[1] <= 2860)

