"""
https://rabernat.github.io/research_computing_2018/maps-with-cartopy.html
https://www.naturalearthdata.com/features/
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime
from kadlu.geospatial.data_sources.chs import Chs
from scipy.interpolate import griddata
#from matplotlib.transforms import offset_copy
import cartopy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature


""" load the data """

kwargs = dict(
        start=datetime(2015, 1, 9), end=datetime(2015, 1, 10, 12),
        south=44,                   west=-64, 
        north=45,                   east=-63, 
        top=0,                      bottom=5000
    )

#Chs().fetch_bathymetry(**kwargs)

val1, lat1, lon1 =  np.array(Chs().load_bathymetry(**kwargs)).astype(float)


""" perform projections on data """

terrain = cimgt.Stamen('terrain')
#terrain = cimgt.OSM()

#proj = ccrs.Miller()
proj = ccrs.Mercator()

extent = proj.transform_points(
        ccrs.Geodetic(),
        np.array([kwargs['west'], kwargs['east']]), 
        np.array([kwargs['south'], kwargs['north']])
        )[:,:-1]
lonlat = proj.transform_points(
        ccrs.Geodetic(),
        lon1,
        lat1
    )
lon = lonlat[:,0]
lat = lonlat[:,1]

###
#num_lats = int(np.ceil((kwargs['north']- kwargs['south']) / 0.001)) + 1
#num_lons = abs(int(np.ceil((kwargs['west']- kwargs['east']) / 0.001)) + 1)
#lons = np.linspace(start=min(lon1), stop=max(lon1), num=num_lons)
#lats = np.linspace(start=min(lat1), stop=max(lat1), num=num_lats)
#data = griddata(points=(lon1, lat1), values=val1, xi=(lons[None,:],lats[:,None]), method='linear')
###
#num_lats = int(np.ceil((extent[1][1] - extent[0][1]) / 0.001)) + 1
#num_lons = abs(int(np.ceil((extent[0][0] -extent[1][0]) / 0.001)) + 1)
num_lats = int(np.ceil((extent[1][1] - extent[0][1]) / 100)) + 1
num_lons = abs(int(np.ceil((extent[0][0] -extent[1][0]) / 100)) + 1)
lons = np.linspace(start=min(lon), stop=max(lon), num=num_lons)
lats = np.linspace(start=min(lat), stop=max(lat), num=num_lats)
data = griddata(points=(lon, lat), values=val1, xi=(lons[None,:],lats[:,None]), method='linear')
###


""" prepare filigree """

coast = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
rivers_lakes = cfeature.NaturalEarthFeature('physical', 'lakes_north_america', '10m')
#cmap='Greys_r'
#cmap='ocean'
#cmap='gray'
#cmap='binary_r'
cmap='bone'

""" draw the plot """

"""
from shapely.ops import cascaded_union
from cartopy.io import shapereader
from shapely.geometry import Point

pts = [Point(xp, yp) for xp, yp in zip(lon, lat)]
mpts = cascaded_union(pts)
mask = mpts.convex_hull
"""

#fig = plt.figure(figsize=(10, 5))
fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax = fig.add_subplot(1, 1, 1, projection=proj)
#ax.add_image(terrain, 11, zorder=5)
#ax = plt.axes(projection=terrain.crs)
#ax.set_extent((extent[0][0], extent[1][0], extent[0][1], extent[1][0]))
#ax.add_image(terrain, 8)
ax.add_wms(wms='http://vmap0.tiles.osgeo.org/wms/vmap0', layers=['basic'], zorder=0)
ax.contourf(lons, lats, data,
            transform=proj,
            cmap=cmap
        )
ax.add_feature(cartopy.feature.LAKES, edgecolor='black', projection=proj)
ax.add_feature(cartopy.feature.RIVERS, edgecolor='black')
#ax.add_feature(cartopy.feature.LAND, edgecolor='black')
#ax.add_feature(coast, facecolor=(0.93,0.93,0.93,1), edgecolor='black')
ax.add_feature(rivers_lakes, facecolor='black', edgecolor='black', projection=proj)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER
#ax.set_xticks([-63.75, -63.5, -63.25, -63], crs=proj)
#ax.set_yticks([47.25, 47.5, 47.75, 48], crs=ccrs.Miller())
ax.set_title('bathymetry (metres)')
norm = matplotlib.colors.Normalize(vmin=min(val1), vmax=max(val1))
plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap))
#ax.set_global()
plt.show()

