import numpy as np

# Equatorial radius (6,378.1370 km)
# Polar radius (6,356.7523 km)
# The International Union of Geodesy and Geophysics 
# (IUGG) defines the mean radius (denoted R1) to be                        
# 6,371.009 km
R1_IUGG = 6371009

# Degree to radian conversion factor
deg2rad = np.pi / 180.


def LLtoXY(lat, lon, lat_ref=0, lon_ref=0, rot=0, grid=False):
    """ Transform lat-lon coordinates to xy position coordinates.

        By default, the origin of the xy coordinate system is 
        set to 0 deg latitude and 0 deg longitude. 
        
        By default, the x-axis is aligned with the longitude axis 
        (west to east) and the y-axis is aligned with the latitude
        axis (south to north).

        Args: 
            lat: float or numpy array
                latitude coordinate of a the location(s) of interest in degrees
            lon: float or numpy array
                longitude coordinate of a the location(s) of interest in degrees
            lat_ref: float
                latitude reference coordinate in degrees
            lon_ref: float
                longitude reference coordinate in degrees
            rot: float
                Rotation angle in degrees for the xy coordinate system 
                (clockwise rotation).

        Returns:
            x: float or numpy array
                x position coordinates in meters
            y: float or numpy array
                y positions coordinates in meters
    """
    global R1_IUGG, deg2rad

    if grid:
        lat, lon = np.meshgrid(lat, lon)

    else:
        lat = np.array([lat])
        lon = np.array([lon])
        if np.ndim(lat) == 2:
            lat = np.squeeze(lat)
            lon = np.squeeze(lon)

        assert lat.shape[0] == lon.shape[0], 'lat and lon must have same length'

    R = R1_IUGG
    R2 = R * np.cos(lat_ref * deg2rad)

    x = (lon - lon_ref) * deg2rad * R2
    y = (lat - lat_ref) * deg2rad * R

    if rot is not 0:
        s = np.sin(rot * deg2rad)
        c = np.cos(rot * deg2rad)
        rotmat = np.array([[c, -s], [s, c]]) 
        xy = np.array([x, y])
        xy = np.swapaxes(xy, 0,1)
        xy = rotmat.dot(xy)
        x = xy[:,0]
        y = xy[:,1]

    if len(x) == 1:
        x = float(x)
        y = float(y)

    return x, y

def XYtoLL(x, y, lat_ref=0, lon_ref=0, rot=0, grid=False):
    """ Transform xy position coordinates to lat-lon coordinates.

        By default, the origin of the xy coordinate system is 
        set to 0 deg latitude and 0 deg longitude. 
        
        By default, the x-axis is aligned with the longitude axis 
        (west to east) and the y-axis is aligned with the latitude
        axis (south to north).

        Args: 
            x: float or numpy array
                x coordinate of a the location(s) of interest in meters
            y: float or numpy array
                y coordinate of a the location(s) of interest in meters
            lat_ref: float
                latitude reference coordinate in degrees
            lon_ref: float
                longitude reference coordinate in degrees
            rot: float
                Rotation angle in degrees for the xy coordinate system 
                (clockwise rotation).

        Returns:
            lat: float or numpy array
                latitude coordinates in degrees
            lon: float or numpy array
                longitude coordinates in degrees
    """
    global R1_IUGG, deg2rad

    if grid:
        x, y = np.meshgrid(x, y)

    else:
        x = np.array([x])
        y = np.array([y])
        if np.ndim(x) == 2:
            x = np.squeeze(x)
            y = np.squeeze(y)

        assert x.shape[0] == y.shape[0], 'x and y must have same length'


    R = R1_IUGG
    R2 = R * np.cos(lat_ref * deg2rad)

    if rot is not 0:
        s = np.sin(rot * deg2rad)
        c = np.cos(rot * deg2rad)
        rotmat = np.array([[c, s], [-s, c]]) 
        xy = np.array([x, y])
        xy = np.swapaxes(xy, 0,1)
        xy = rotmat.dot(xy)
        x = xy[:,0]
        y = xy[:,1]

    lon = lon_ref + x / deg2rad / R2
    lat = lat_ref + y / deg2rad / R

    if len(lat) == 1:
        lat = float(lat)
        lon = float(lon)

    return lat, lon

def regXYgrid(lat, lon, lat_ref=0, lon_ref=0, rebin=1):
    """ Transform a regular polar grid (lat,lon) into a regular planar 
        grid (x,y) for the same area.

        Args: 
            lat: array
                latitudes in degrees from -90 to +90.
            lon: array
                longitudes in degrees from -180 to +180.
            lat_ref: float
                reference latitude defining the origin of the y axis.
            lon_ref: float
                reference longitude defining the origin of the x axis.
            rebin: int
                Re-binning factor for the x-y grid. For example, if 
                rebin=4, the binning used for the x and y axes will 
                be 4 times finer than that of the longitude and 
                latitude axes, respectively. 

        Returns:
            x: array
                x coordinates of the planar grid
            y: array
                y coordinates of the planar grid
    """
    Nx = rebin * len(lon)
    Ny = rebin * len(lat)

    # determine extent and bin size of regular x-y grid
    lats = [lat[0], lat[0], lat[-1], lat[-1]]
    lons = [lon[0], lon[-1], lon[0], lon[-1]]

    x, y = LLtoXY(lat=lats, lon=lons, lat_ref=lat_ref, lon_ref=lon_ref)

    x_length_S = x[1] - x[0]
    x_length_N = x[3] - x[2]

    if x_length_S < x_length_N:
        x0 = x[0]
        dx = (x[1] - x[0]) / (Nx - 1)
    else:
        x0 = x[2]
        dx = (x[3] - x[2]) / (Nx - 1)

    y0 = y[0]
    dy = (y[2] - y[0]) / (Ny - 1)

    # create regular x-y grid
    x = np.arange(Nx, dtype=np.float)
    x *= dx
    x += x0
    y = np.arange(Ny, dtype=np.float)
    y *= dy
    y += y0

    return x, y