import pytest
from datetime import datetime, timedelta
from kadlu import hycom
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
start = datetime(2000, 1, 2)
end   = datetime(2000, 1, 2, 4)

# disable automated testing of fetching to avoid slamming the API with
# requests in the automated development pipeline
test_fetch = True

# smaller test area for faster tests
south =  46 
north =  48
west  = -60
east  = -58
top   = 0
bottom = 50

def test_fetch_salinity():
    if not test_fetch: return
    Hycom().fetch_salinity(south=south, north=north, west=west, east=east, start=start, end=end, top=top, bottom=bottom)
    """

    fetch_water_u(Hycom(), south=south, north=north, west=west, east=east, top=top, bottom=bottom, start=start, end=end)

    def fetch_salinity(self, **kwargs):
        return kwargs
    qry = fetch_salinity(Hycom(), south=south, north=north, west=west, east=east, start=start, end=end, top=top, bottom=bottom)


    fetch_hycom(year=year, slices=slices, fetchvar=fetchvar, lat=lat, lon=lon, epoch=epoch, depth=depth)



    qry = {'south': 46.0, 'north': 48.0, 'west': 175.0, 'east': 179.9200439453125, 'top': 0, 'bottom': 1, 'start': datetime(2000, 1, 2, 0, 0), 'end': datetime(2000, 1, 2, 4, 0)}
    """

def test_load_salinity():
    val, lat, lon, time, depth = Hycom().load_salinity(south=south, north=north, west=west, east=east, start=start, end=end)
    return
    assert(len(val) == len(lat) == len(lon) == len(time))
    assert(sum(lat <= 90) == sum(lat >= -90) == len(lat))
    assert(sum(lon <= 180) == sum(lon >= -180) == len(lon))

def test_fetch_temp():
    if not test_fetch: return
    Hycom().fetch_temp(south=south, north=north, west=west, east=east, start=start, end=end, top=top, bottom=bottom)

def test_load_temp():
    val, lat, lon, time, depth = Hycom().load_temp(south=south, north=north, west=west, east=east, start=start, end=end, top=top, bottom=bottom)
    assert(len(val) == len(lat) == len(lon) == len(time))
    assert(sum(lat <= 90) == sum(lat >= -90) == len(lat))
    assert(sum(lon <= 180) == sum(lon >= -180) == len(lon))

def test_fetch_water_u():
    if not test_fetch: return
    Hycom().fetch_water_u(south=south, north=north, west=west, east=east, start=start, end=end, top=top, bottom=bottom)

def test_load_water_u():
    val, lat, lon, time, depth = Hycom().load_water_u(south=south, north=north, west=west, east=east, start=start, end=end, top=top, bottom=bottom)
    assert(len(val) == len(lat) == len(lon) == len(time))
    assert(sum(lat <= 90) == sum(lat >= -90) == len(lat))
    assert(sum(lon <= 180) == sum(lon >= -180) == len(lon))

def test_fetch_water_v():
    if not test_fetch: return
    Hycom().fetch_water_v(south=south, north=north, west=west, east=east, start=start, end=end, top=top, bottom=bottom)

def test_load_water_v():
    val, lat, lon, time, depth = Hycom().load_water_u(south=south, north=north, west=west, east=east, start=start, end=end, top=top, bottom=bottom)
    assert(len(val) == len(lat) == len(lon) == len(time))
    assert(sum(lat <= 90) == sum(lat >= -90) == len(lat))
    assert(sum(lon <= 180) == sum(lon >= -180) == len(lon))


""" interactive mode debugging: assert db ordering is correct

    step through fetch_hycom() and put output and grid arrays into memory. 
    example test input:
>>>     
        year = '2000'
        fetchvar = 'salinity'
        slices = [
            (0, 2),         # time: start, end 
            (0, 3),         # depth: top, bottom
            (800, 840),     # x grid index: lon min, lon max
            (900, 1000)     # y grid index: lat min, lat max
        ]
        lat, lon = load_grid()
        epoch = load_times()
        depth = load_depth()

    run through the output builder loop again. this time, add an assertion to check 
    that the 4D array was flattened correctly
>>>     
        ix = 0  # debug index: assert order is correct
        for arr in payload.split("\n"):
            ix_str, row_csv = arr.split(", ", 1)
            a, b, c = [int(x) for x in ix_str[1:-1].split("][")]
            # output[a][b][c] = np.array(row_csv.split(", "), dtype=np.int)
            assert((output[a][b][c] == grid[ix:ix+len(output[a][b][c]), 0]).all())
            ix += len(output[a][b][c])

"""
"""
    import timeit 

    lat, sorted_arr = load_grid()
    def fcn():
        return index(np.random.rand()*360 - 180, sorted_arr)
    
    timeit.timeit(fcn, number=1000000)
    
"""

