from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import Rbf
import numpy as np


def f(x,y,z):
    return 2 * x**3 + 3 * y**2 - z

N = 10
r = 10.
pts = np.array([[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]])
x = y = z = np.linspace(0, r, N)


# RegularGridInterpolator

data = f(*np.meshgrid(x, y, z, indexing='ij', sparse=True))
rgi = RegularGridInterpolator((x, y, z), data)
di = rgi(pts)
print(di)


# Rbf

X,Y,Z = np.meshgrid(x,y,z,indexing='ij')

d = f(X,Y,Z)

X =  X.flatten()
Y =  Y.flatten()
Z =  Z.flatten()
d = d.flatten()

rbfi = Rbf(X, Y, Z, d, function='cubic', smooth=0)  # radial basis function interpolator instance

di = rbfi(pts[:,0], pts[:,1], pts[:,2])   # interpolated values

print(di)