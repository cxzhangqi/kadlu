from osgeo import gdal
import numpy as np

path = '/home/oliskir/Meridian/CHS/CHS_NONNA_100/CA2_4300N06000W.tif'

# open dataset
ds = gdal.Open(path)

print(ds.RasterCount)
print(ds.GetMetadata())

band = ds.GetRasterBand(1)

arr = ds.ReadAsArray()
print(arr.shape)


# Get nodata value from the GDAL band object
nodata = band.GetNoDataValue()

#Create a masked array for making calculations without nodata values
arr = np.ma.masked_equal(arr, nodata)
type(arr)

print(arr.min())

stats = band.GetStatistics( True, True )

print("[ STATS ] =  Minimum={0:.3f}, Maximum={1:.3f}, Mean={2:.3f}, StdDev={3:.3f}".format( \
                stats[0], stats[1], stats[2], stats[3] ))



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm

M = arr.shape[0]
N = arr.shape[1]

X = np.arange(M)
Y = np.arange(N)
X, Y = np.meshgrid(X, Y)

arr = np.flip(arr, axis=0)


fig, ax = plt.subplots()

p = ax.pcolor(X, Y, arr[:M,:N])
cb = fig.colorbar(p)

plt.show()


# close dataset
ds = None
band = None
