import numpy as np
from datetime import datetime
from kadlu.geospatial.interpolation             import     \
        Interpolator2D,                                    \
        Interpolator3D,                                    \
        Uniform2D,                                         \
        Uniform3D
from kadlu.geospatial.data_sources.data_util    import index
from kadlu.geospatial.data_sources.chs          import Chs
from kadlu.geospatial.data_sources.hycom        import Hycom
from kadlu.geospatial.data_sources.era5         import Era5
from kadlu.geospatial.data_sources.wwiii        import Wwiii
#import kadlu.geospatial.data_sources.gebco as gebco 
#from kadlu.utils import LatLon


def flatten(cols, frame_ix):
    """ reduce 4D to 3D by averaging over time dimension """
    assert reduce(lambda a, b: (a==b)*a, frame_ix[1:] - frame_ix[:-1])

    ix = range(0, len(frame_ix) -1)
    frames = np.array([cols[0][frame_ix[f] : frame_ix[f +1]] for f in ix])
    vals = (reduce(np.add, frames) / len(frames))
    _, y, x, _, z = cols[:, frame_ix[0] : frame_ix[1]]

    warnings.warn("query data has been averaged across the time dimension "
                  "for 3D interpolation.\nto avoid this behaviour, "
                  "use keyword argument 'time' instead of start/end")

    return vals, y, x, frames, z


def reshape_3D(callback, **kwargs):
    """ load 3D data from database and prepare it for interpolation """
    cols = callback(**kwargs).astype(np.float)
    frame_ix = np.append(np.nonzero(cols[3][1:] > cols[3][:-1])[0] + 1, len(cols[3]))
    vals, y, x, _, z = flatten(cols, frame_ix) if len(frame_ix) > 1 else cols
    rows = np.array((vals, y, x, z)).T

    # reshape row data to 3D array
    xgrid, ygrid, zgrid = np.unique(x), np.unique(y), np.unique(z)
    gridspace = np.full((len(ygrid), len(xgrid), len(zgrid)), fill_value=-30000)
    # this could be optimized to avoid an index lookup cost maybe
    for row in rows:
        x_ix = index(row[2], xgrid)
        y_ix = index(row[1], ygrid)
        z_ix = index(row[3], zgrid)
        gridspace[y_ix, x_ix, z_ix] = row[0]

    # remove -30000 values for interpolation:
    # fill missing depth values with last value in each column
    # this section could be cleaned up
    for xi in range(0, gridspace.shape[0]):
        for yi in range(0, gridspace.shape[1]):
            col = gridspace[xi, yi]
            if sum(col == -30000) > 0 and sum(col == -30000) < len(col):
                col[col == -30000] = col[col != -30000][-1]
                gridspace[xi, yi] = col

    # TODO:
    # create default values for columns that are entirely null

    return dict(
            values=gridspace, lats=ygrid, lons=xgrid, depths=zgrid, 
            origin=kwargs['origin'], method=kwargs['method']
        )
                    

def reshape_2D(callback, **kwargs):
    """ load 2D data from database and prepare it for interpolation """
    # do something similar to reshape_3D but for 2D 
    pass
    

# TODO:
# add params to output kwargs before passing to interpolator, possibly by
# configuration in the config.ini file. 
# this function is a temporary and quick fix
def default(kwargs):
    kwargs['water_density']=1.0
    kwargs['origin'] = 0, 0 #lat_ref, lon_ref ???
    kwargs['method'] = 'linear'
    return kwargs


class Ocean():
    def __init__(self):
        # later we can define the preferred data sources on class initialization
        # or config.ini. for now, just use these ones
        self.load_bathymetry = Chs().load_bathymetry
        self.load_temp = Hycom().load_temp
        self.load_salinity = Hycom().load_salinity
        self.load_wave = Era5().load_windwaveswellheight

    def interp_bathy(self, **kwargs):
        kwargs = default(kwargs)  # temporary fix
        return Interpolator2D(**reshape_2D(self.load_bathymetry, **kwargs))

    def interp_temp(self, **kwargs):
        kwargs = default(kwargs)  # temporary fix
        return Interpolator3D(**reshape_3D(self.load_temp, **kwargs))

    def interp_salinity(self, **kwargs):
        kwargs = default(kwargs)  # temporary fix
        return Interpolator3D(**reshape_3D(self.load_salinity, **kwargs))

    def interp_wave(self, **kwargs):
        kwargs = default(kwargs)  # temporary fix
        return Interpolator2D(**reshape_2D(self.load_wave, **kwargs))


    # functions for evaluating the interpolated data go here
    def bathy():
        pass

    def bathy_gradient():
        pass

    def temp():
        pass

    def salinity():
        pass

    def wave():
        # the wave variable to be used here can be selected in the class 
        # initialization, or possibly in the config.ini file
        pass

    def wind_speed():
        pass

