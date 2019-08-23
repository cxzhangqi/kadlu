import numpy as np
from pyproj import Proj, transform
from kadlu.geospatial.data_sources import chs
from osgeo import gdal

# convert from spherical mercator to WGS 84 (global)
inproj = Proj(init="epsg:3857")
outproj = Proj(init="epsg:4326")

# fetch CHS data, filter out .tfw files, grab the first .tif
fpaths = chs.fetch(south=45, west=-84, north=47, east=-81)
mask = ["tfw" not in x for x in fpaths]
tifs = np.array(fpaths)[mask]
bathy = gdal.Open(tifs[0])

# read values and convert
# no idea what z1 and z2 are, looks like all zeroes ??
x1, x1binsize, z1, y1, z2, y1binsize = bathy.GetGeoTransform()
x2, y2 = transform(inproj, outproj, x1, y1)
print(x2, y2)
