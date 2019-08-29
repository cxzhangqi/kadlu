import numpy as np
import requests
from kadlu.geospatial.data_sources import fetch_util


def fetch(south=-90, north=90, west=-180, east=180):
    """
        Fetches salinity for the given area.

        TODO:
            allow time range input
            allow depth input ???
            fetch temperature
    """
    year = "2015"
    source = f"https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/{year}.ascii?"

    lat, lon = load_grid()

    t = (0, 2)      # time
    d = (0, 3)      # depth
    x = (ll_2_xy(west, lon),    ll_2_xy(east, lon))
    y = (ll_2_xy(south, lat),   ll_2_xy(north, lat))

    slx         = lambda tup, step=1: f"[{tup[0]}:{step}:{tup[1]}]"
    varslices   = lambda var, slxs  : f"{var}{''.join([slx(v) for v in slxs])}"

    salinity_txt = requests.get(f"{source}{varslices('salinity', [t, d, x, y])}")
    assert(salinity_txt.status_code == 200)

    meta, data = salinity_txt.text.split("---------------------------------------------\n")
    arrs = data.split("\n\n")[:-1]

    shape_str, payload = arrs[0].split("\n", 1)
    shape = tuple([int(x) for x in shape_str.split("[", 1)[1][:-1].split("][")])
    salinity = np.ndarray(shape, dtype=np.int)
    for arr in payload.split("\n"):
        ix_str, row_csv = arr.split(", ", 1)
        a, b, c = [int(x) for x in ix_str[1:-1].split("][")]
        salinity[a][b][c] = np.array(row_csv.split(", "), dtype=np.int)

    # not sure if we need the other values also
    #for arr in arrs[1:]:
    #    header, vals = arr.split("\n")
    #    print(f"{header}:\t{vals}")

    return salinity, lat[y[0]:y[1]], lon[x[0]:x[1]]


def fetch_grid():
    print("Fetching latitude and longitude grid from HYCOM...")
    url = "https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015.ascii?lat%5B0:1:3250%5D,lon%5B0:1:4499%5D"
    grid_txt = requests.get(url)
    assert(grid_txt.status_code == 200)

    meta, data = grid_txt.text.split("---------------------------------------------\n")
    lat_csv, lon_csv = data.split("\n\n")[:-1]
    lat = np.array(lat_csv.split("\n")[1].split(", "), dtype=np.float)
    lon = np.array(lon_csv.split("\n")[1].split(", "), dtype=np.float)

    stor_loc = fetch_util.instantiate_storage_config()
    np.save(f"{stor_loc}hycom_lats.npy", lat, allow_pickle=False)
    np.save(f"{stor_loc}hycom_lons.npy", lon, allow_pickle=False)


def load_grid():
    storage = fetch_util.instantiate_storage_config()
    def loadLatLon(storage): return np.load(f"{storage}hycom_lats.npy"), np.load(f"{storage}hycom_lons.npy")
    try:
        return loadLatLon(storage)
    except FileNotFoundError:
        fetch_grid()
        return loadLatLon(storage)


def ll_2_xy(key, array): 
    """ converts lat/lon values to grid index """
    if key > array[-1]: return len(array) - 1
    return max(np.nonzero(array >= key)[0][0] - 1, 0)

