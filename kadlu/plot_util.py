import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import griddata
import cartopy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from scipy.spatial import ConvexHull

"""

from kadlu.geospatial.data_sources.chs import Chs
kwargs = dict(
        start=datetime(2015, 1, 9), end=datetime(2015, 1, 10, 12),
        south=44,                   west=-64, 
        north=45,                   east=-63, 
        top=0,                      bottom=5000
    )
#Chs().fetch_bathymetry(**kwargs)
val, lat, lon =  np.array(Chs().load_bathymetry(**kwargs)).astype(float)

"""

proj = ccrs.Mercator()

#terrain = cimgt.Stamen('terrain')
#terrain = cimgt.StamenTerrain()
#terrain = cimgt.Stamen('toner-lite')
#terrain = cimgt.MapQuestOpenAerial()

#cmap='Greys_r'
#cmap='ocean'
#cmap='gray'
#cmap='binary_r'
cmap=plt.cm.bone
#cmap.set_over((0,0,0,0))


def plot2D(val, lat, lon, title='bathymetry (metres)', **kwargs):
    # project data onto coordinate space
    extent = proj.transform_points(
            ccrs.Geodetic(),
            np.array([kwargs['west'], kwargs['east']]), 
            np.array([kwargs['south'], kwargs['north']])
        )[:,:-1]
    projected_lonlat = proj.transform_points(
            ccrs.Geodetic(),
            lon,
            lat
        )

    plon = projected_lonlat[:,0]
    plat = projected_lonlat[:,1]
    

    # create interpolation grid, perform interpolation
    num_lats = int(np.ceil((extent[1][1] - extent[0][1]) / 100)) + 1
    num_lons = abs(int(np.ceil((extent[0][0] -extent[1][0]) / 100)) + 1)
    lons = np.linspace(start=min(plon), stop=max(plon), num=num_lons)
    lats = np.linspace(start=min(plat), stop=max(plat), num=num_lats)
    data = griddata(points=(plon, plat), values=val, xi=(lons[None,:],lats[:,None]), method='linear')

    # map filigree
    coast = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
    #rivers_lakes = cfeature.NaturalEarthFeature('physical', 'lakes_north_america', '10m')

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    #ax.add_image(terrain, 9)
    #ax.add_wms(wms='http://vmap0.tiles.osgeo.org/wms/vmap0', layers=['basic'])
    ax.add_feature(coast, facecolor=(0.89,0.92,0.93,0.8), edgecolor=(0,0,0,0.6))
    #ax.add_feature(rivers_lakes, facecolor=(0.89,0.92,0.94,1), edgecolor='black')
    ax.contourf(lons, lats, data,
                transform=proj,
                levels=range(int(min(val)), -1, 12),
                cmap=cmap
            )
    ax.add_feature(cartopy.feature.LAKES, edgecolor='black', projection=proj)
    ax.add_feature(cartopy.feature.RIVERS, edgecolor='black')
    ax.add_feature(rivers_lakes, facecolor='black', edgecolor='black', projection=proj)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
    gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER
    #ax.set_xticks([-63.75, -63.5, -63.25, -63], crs=proj)
    #ax.set_yticks([47.25, 47.5, 47.75, 48], crs=ccrs.Miller())
    ax.set_title(title)
    norm = matplotlib.colors.Normalize(vmin=min(val), vmax=max(val))
    plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap))
    #ax.set_global()
    plt.show()

    return 




