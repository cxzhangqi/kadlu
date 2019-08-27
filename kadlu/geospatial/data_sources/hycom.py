import numpy as np
import requests

def fetch():
    # in the future there will be a function to abstract time, depth, lats, 
    # and lons based on inputs to this function 
    # matt_s 2019-08
    year = "2015"
    source = f"https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/{year}.ascii?"

    t = (0, 2)      # time
    d = (0, 3)      # depth
    x = (900, 940)  # lat
    y = (800, 830)  # lon

    slx = lambda tup, step=1: f"[{tup[0]}:{step}:{tup[1]}]"
    varslices = lambda var, slxs : f"{var}{''.join([slx(v) for v in slxs])}"
    salinity_txt = requests.get(f"{source}{varslices('salinity', [t, d, x, y])}")
    assert(salinity_txt.status_code == 200)

    meta, data = salinity_txt.text.split("---------------------------------------------\n")
    arrs = data.split("\n\n")[:-1]

    shape_str, payload = arrs[0].split("\n", 1)
    shape = tuple([int(x) for x in shape_str.split("[", 1)[1][:-1].split("][")])
    salinity = np.ndarray(shape, dtype=np.int)
    for arr in payload.split("\n"):
        ix_str, row_csv = arr.split(", ", 1)
        ix = [int(x) for x in ix_str[1:-1].split("][")]
        salinity[ix[0]][ix[1]][ix[2]] = np.array(row_csv.split(", "), dtype=np.int)
        print(f"{ix}\t{row}")
        break

    #for arr in arrs[1:]:
    #    header, vals = arr.split("\n")
    #    print(f"{header}:\t{vals}")

    return salinity
