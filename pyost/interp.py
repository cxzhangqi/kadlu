
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
