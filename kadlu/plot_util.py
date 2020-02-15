import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import griddata
import cartopy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature



#terrain = cimgt.Stamen('terrain')
#terrain = cimgt.StamenTerrain()
#terrain = cimgt.Stamen('toner-lite')
#terrain = cimgt.MapQuestOpenAerial()

#cmap='Greys_r'
#cmap='ocean'
#cmap='gray'
#cmap='binary_r'

#cmap=plt.cm.ocean
#cmap=plt.cm.bone
#cmap=plt.cm.binary_r
#cmap=plt.cm.viridis


proj = ccrs.Mercator()
config = dict(
        bgcontour   = lambda v: np.linspace(min(v)-.1, max(v)+.1, 3),
        bathy=dict(
            cm      = plt.cm.bone,
            alpha   = 0.9,
            levels  = lambda v, n=25: np.linspace(min(v), max(v)-2, n),
            title   = 'bathymetry (metres)'),
        temp=dict(
            #cm      = plt.cm.jet, 
            cm      = plt.cm.coolwarm, 
            alpha   = 0.8,
            levels  = lambda v, n=25: np.linspace(min(v)+.1, max(v), n),
            title   = 'temperature (celsius)'),
        salinity=dict(
            cm      = plt.cm.viridis,
            alpha   = 0.7,
            levels  = lambda v, n=25: np.linspace(min(v)+.1, max(v), n),
            title   = 'salinity (g/kg of salt in water)'),
        #wavedir=dict(
        #    cm=
        waveheight=dict(
            #cm      = plt.cm.Spectral_r,
            #cm      = plt.cm.BrBG,
            cm      = plt.cm.BuPu,
            #cm      = sns.cubehelix_palette(8, start=2, rot=0, dark=0, light=.95, reverse=True),
            alpha   = 0.85,
            levels  = lambda v, n=20: np.linspace(min(v)+.1, max(v), n),
            title   = 'wave height (metres)')
    )



def plot2D(val, lat, lon, var, title='bathymetry (metres)', **kwargs):
    """

    from kadlu.geospatial.data_sources.chs import Chs
    from kadlu.geospatial.data_sources.hycom import Hycom
    from kadlu.geospatial.data_sources.era5 import Era5
    kwargs = dict(
            start=datetime(2015, 1, 9), end=datetime(2015, 1, 9, 3),
            south=44,                   west=-64.5, 
            north=45,                   east=-62.5, 
            top=0,                      bottom=0
        )
    kwargs = dict(
            start=datetime(2015, 1, 9), end=datetime(2015, 1, 9, 3),
            south=45,                   west=-68.4, 
            north=51.5,                 east=-56.5, 
            top=0,                      bottom=5000
        )

    #Chs().fetch_bathymetry(**kwargs)
    #Hycom().fetch_temp(**kwargs)
    #Hycom().fetch_salinity(**kwargs)
    #Era5().fetch_windwaveswellheight(**kwargs)
    val, lat, lon, time = Era5().load_windwaveswellheight(**kwargs)
    var = 'waveheight'
    val, lat, lon, time, depth =  Hycom().load_temp(**kwargs)
    var = 'temp'
    val, lat, lon, time, depth =  Hycom().load_salinity(**kwargs)
    var = 'salinity'
    val, lat, lon =  Chs().load_bathymetry(**kwargs)
    val, lat, lon = [c[::1000] for c in Chs().load_bathymetry(**kwargs)]
    var = 'bathy'

    """
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
    #num_lats = int(np.ceil((extent[1][1] - extent[0][1]) / 100)) + 1
    #num_lons = abs(int(np.ceil((extent[0][0] -extent[1][0]) / 100)) + 1)
    num_lats = 1000
    num_lons = 1000
    lons = np.linspace(start=min(plon), stop=max(plon), num=num_lons)
    lats = np.linspace(start=min(plat), stop=max(plat), num=num_lats)
    data = griddata(points=(plon, plat), values=val, xi=(lons[None,:],lats[:,None]), method='linear')
    # map content settings
    coast = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
    #ocean = cfeature.NaturalEarthFeature('physical', 'ocean', '10m')
    norm = matplotlib.colors.Normalize(vmin=min(val), vmax=max(val))
    fg = (.92, .92, .92, 1)

    cm = config[var]['cm']
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, 
            title=config[var]['title'],
            projection=proj, 
            facecolor=cm(256), 
            frameon=True
        )
    ax.contourf(lons, lats, data,
                transform=proj,
                levels=config[var]['levels'](val),
                cmap=cm,
                alpha=config[var]['alpha'],
                zorder=8
            )
    ax.contour(lons, lats, data,
                transform=proj,
                levels=config[var]['levels'](val),
                cmap=cm,
                alpha=1,
                linewidths=2,
                zorder=9
            )
    #ax.add_feature(ocean, facecolor=cm(0), edgecolor=(0,0,0,1), zorder=10)
    ax.add_feature(coast, facecolor=fg, edgecolor=(0,0,0,1), zorder=10)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', zorder=11)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
    gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER
    plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cm))
    manager = plt.get_current_fig_manager()
    s = manager.window.maximumSize()
    manager.resize(s.width(), s.height())
    plt.tight_layout()
    plt.show()

    return 



def fetch_topo():
    url = 'ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/shp/Elevation/canvec_50K_NS_Elevation_shp.zip'
    with requests.get(url, stream=True) as payload:
        assert payload.status_code == 200, 'error bad request'
        fname = storage_cfg()+'topo_NS.zip'
        with open(fname, 'wb') as f: f.write(payload.content)
        


#ax.add_feature(cartopy.feature.LAKES, edgecolor='black', projection=proj)
#ax.add_feature(cartopy.feature.RIVERS, edgecolor='black')
#ax.add_feature(rivers_lakes, facecolor='black', edgecolor='black', projection=proj)

#ax.set_xticks([-63.75, -63.5, -63.25, -63], crs=proj)
#ax.set_yticks([47.25, 47.5, 47.75, 48], crs=ccrs.Miller())

#ax.add_image(terrain, 9)
#ax.add_wms(wms='http://vmap0.tiles.osgeo.org/wms/vmap0', layers=['basic'])
#ax.add_feature(rivers_lakes, facecolor=(0.89,0.92,0.94,1), edgecolor='black')
#ax.add_feature(coast, facecolor=(0.89,0.92,0.93,0.8), edgecolor=(0,0,0,0.6))

