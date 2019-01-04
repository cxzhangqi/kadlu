
# According to 
# https://stackoverflow.com/questions/31524188/speeding-up-scipys-rectbivariatespline-evaluation-function
# the ev() function is much slower than direct call to __call__ with grid=True

from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline
import numpy as np
import time


d = Dataset('../res/GEBCO_2014_2D.nc')

Nx=len(d.variables['lat'])
Ny=len(d.variables['lon'])
dN=10

x = d.variables['lat'][int(Nx/2-dN):int(Nx/2+dN)]
y = d.variables['lon'][int(Ny/2-dN):int(Ny/2+dN)]
z = d.variables['elevation'][int(Nx/2-dN):int(Nx/2+dN),int(Ny/2-dN):int(Ny/2+dN)]

bs = RectBivariateSpline(x=x, y=y, z=z)

xi = np.arange(x[0], x[-1], 0.00005)
yi = np.arange(y[0], y[-1], 0.00005)


start = time.time()

xim, yim = np.meshgrid(xi, yi)
zim = bs.ev(xi=xim, yi=yim)

end = time.time()

print(zim.shape)
print('ev time: ', end - start)


start = time.time()

zi = bs.__call__(x=xi, y=yi, grid=True)  # this method is faster by a factor of 10 !

end = time.time()

print(zi.shape)
print('call time: ', end - start)


from bathy_reader import BathyNetCDFReader, LatLon

reader = BathyNetCDFReader(path='../res/GEBCO_2014_2D.nc')
lat,lon,bathy = reader.read(latlon_SW=LatLon(54.95,14.55), latlon_NE=LatLon(55.35,15.25))
print(np.min(bathy), np.max(bathy))
print(np.min(lat), np.max(lat))

rootgrp = Dataset("bornholm.nc", "w", format="NETCDF4")
lat_dim = rootgrp.createDimension("lat", lat.shape[0])
lon_dim = rootgrp.createDimension("lon", lon.shape[0])

latitudes = rootgrp.createVariable("lat","f4",("lat",))
longitudes = rootgrp.createVariable("lon","f4",("lon",))
bathymetries = rootgrp.createVariable("bathy","f4",("lat","lon",))

latitudes[:] = lat
longitudes[:] = lon
bathymetries[:,:] = bathy

rootgrp.description = "Bathymetry for the island of Bornholm, Denmark."

rootgrp.close()

exit(1)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=10)
plt.show()
