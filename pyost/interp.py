
# According to 
# https://stackoverflow.com/questions/31524188/speeding-up-scipys-rectbivariatespline-evaluation-function
# the ev() function is much slower than direct call to __call__ with grid=True


from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline

d = Dataset('../res/GEBCO_2014_2D.nc')

Nx=len(d.variables['lat'])
Ny=len(d.variables['lon'])
dN=10

x = d.variables['lat'][int(Nx/2-dN):int(Nx/2+dN)]
y = d.variables['lon'][int(Ny/2-dN):int(Ny/2+dN)]
z = d.variables['elevation'][int(Nx/2-dN):int(Nx/2+dN),int(Ny/2-dN):int(Ny/2+dN)]

bs = RectBivariateSpline(x=x, y=y, z=z)




import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm
import numpy as np


fig = plt.figure()

ax = fig.add_subplot(1, 1, 1, projection='3d')

xi = np.arange(x[0], x[-1], 0.002)
yi = np.arange(y[0], y[-1], 0.002)
xi, yi = np.meshgrid(xi, yi)

zi = bs.ev(xi=xi, yi=yi)

surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
#ax.set_zlim(-, 1.01)
fig.colorbar(surf, shrink=0.5, aspect=10)

plt.show()
