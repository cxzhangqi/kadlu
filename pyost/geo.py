from geopy.distance import lonlat, distance



a = (10, 45)  # lon, lat
b = (20, 45)
c = (10, 46)
aa = lonlat(*a)
bb = lonlat(*b)
d = distance(aa, bb)
print(d)



import numpy as np

x_ref = 0
y_ref = 0

lat_ref = a[1]
lon_ref = a[0]
lat = np.array([b[1], c[1]])
lon = np.array([b[0], c[0]])

rot = 90
deg2rad = np.pi / 180.
R = 6371009 # mean earth radius in meters
R2 = R * np.cos(lat_ref * deg2rad)

x = (lon - lon_ref) * deg2rad * R2
y = (lat - lat_ref) * deg2rad * R

SINROT = np.sin(rot * deg2rad)
COSROT = np.cos(rot * deg2rad)
rotmat = np.array([[COSROT, -SINROT], [SINROT, COSROT]]) 
xy = np.array([x, y])
xy = np.swapaxes(xy, 0,1)
#print('\nxy: ', xy)
XY = rotmat.dot(xy)
#print('\n ',XY)
x = XY[:,0]
y = XY[:,1]

x = x + x_ref
y = y + y_ref

print(x/1e3, y/1e3)
