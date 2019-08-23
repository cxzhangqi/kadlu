import numpy as np
from pyproj import Proj, transform
from kadlu.geospatial.data_sources import chs
from osgeo import gdal

# converting from spherical mercator to WGS 84 (global)
inproj = Proj(init="epsg:3857")
outproj = Proj(init="epsg:4326")

south=45.1549 
north=46.6209 
west=-84.164 
east=-82.3553

# fetch CHS data, filter out .tfw files, grab the first .tif
fpaths = chs.fetch(south=45.1549, north=46.6209, west=-84.164, east=-82.3553)
mask = ["tfw" not in x for x in fpaths]
tifs = np.array(fpaths)[mask]

# read values and convert
# no idea what z1 and z2 are, looks like all zeroes ??
for tif in tifs:
    bathy = gdal.Open(tif)
    x1, x1binsize, z1, y1, z2, y1binsize = bathy.GetGeoTransform()
    if "CA2" in tif.split('/')[-1]: x2, y2 = x1, y1
    else: x2, y2 = transform(inproj, outproj, x1, y1)
    print(f"lat: {y2:.5f}\tlon: {x2:.5f}\t\tfile: {tif.split('/')[-1]}")
