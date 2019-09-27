from kadlu.bathy_reader import BathyReader

filename = 'assets/BathyData_Mariana_500kmx500km.mat'

reader = BathyReader(input=filename, lat_name='latgrat', lon_name='longrat', bathy_name='mat')

lat, lon, bathy = reader.read()

import numpy as np

# obtain min and max values
lat_min = np.min(lat)
lat_max = np.max(lat)
lon_min = np.min(lon)
lon_max = np.max(lon)
bathy_min = np.min(bathy)

# print in an easily readable format
print('Latitudes: {0:.2f}-{1:.2f} degrees north'.format(lat_min, lat_max))
print('Longitudes: {0:.2f}-{1:.2f} degress east'.format(lon_min, lon_max))
print('Deepest point: {0:.1f} meters below sea surface'.format(-bathy_min))

from kadlu.bathy_reader import LatLon

# specify rectangular bounding box
SW = LatLon(10.0, 141.0)
NE = LatLon(11.0, 142.0)

# read data for selected region
lat, lon, bathy = reader.read(SW, NE)

# print deepest point for selected region
print('Deepest point: {0:.1f} meters below sea surface'.format(-np.min(bathy)))


### Tutorial 3: Interpolate bathy data ### 

from kadlu.bathy_interpolator import BathyInterpolator
#interp = BathyInterpolator(bathy_reader=reader)
refloc = LatLon(9.0, 140)  # reference location at 9 deg N and 140 deg E
interp = BathyInterpolator(bathy_reader=reader, origin=refloc)

# nodes of the latitude-longitude grid
lats = interp.lat_nodes
lons = interp.lon_nodes

# print the min and max values and the number of nodes.
# (note that the nodes are always strictly ascending)
print('Lat range: {0:.2f} to {1:.2f} degrees north ({2} nodes)'.format(lats[0], lats[-1], len(lats)))
print('Lon range: {0:.2f} to {1:.2f} degrees east ({2} nodes)'.format(lons[0], lons[-1], len(lons)))

depth_ll = interp.eval_ll(lat=9, lon=140) # using polar coordinates
depth_xy = interp.eval_xy(x=0, y=0) # using planar coordinates

# check that the two interpolation schemes agree
print('Polar: {0:.1f} m'.format(depth_ll))
print('Planar: {0:.1f} m'.format(depth_xy))

x = [0, 1000, 2000]  # x coordinates in meters
y = [6000, 7000, 8000]  # y coordinates in meters
bathy_i = interp.eval_xy(x=x, y=y)  # evaluate bathymetry at positions (x[i],y[i]); i=0,1,2
bathy_ij = interp.eval_xy(x=x, y=y, grid=True)  # evaluate bathymetry at positions (x[i],y[j]); i=0,1,2; j=0,1,2

print(bathy_i,'\n')
print(bathy_ij)

import matplotlib.pyplot as plt
fig = interp.plot_ll()  # draw elevation map using polar coordinates
plt.show()

fig = interp.plot_xy()  # draw elevation map using planar coordinates
plt.show()

fig = interp.plot_xy(x_bins=10, x_min=100e3, x_max=300e3)  # plot only x-values between 100 and 300 km, divided into 10 bins
plt.show()
