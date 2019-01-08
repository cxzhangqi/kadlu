import numpy as np

R1_IUGG = 6371009
# Equatorial radius (6,378.1370 km)
# Polar radius (6,356.7523 km)
# The International Union of Geodesy and Geophysics 
# (IUGG) defines the mean radius (denoted R1) to be                        
# 6,371.009 km


def LLtoXY(lat, lon, lat_ref=0, lon_ref=0, rot=0):
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
    global R1_IUGG

    assert np.ndim(lat) == np.ndim(lon), 'lat and lon must have same dimension'

    if np.ndim(lat) == 0:
        lat = np.array([lat])
        lon = np.array([lon])
    else:
        lat = np.array(lat)
        lon = np.array(lon)

    deg2rad = np.pi / 180.

    R = R1_IUGG
    R2 = R * np.cos(lat_ref * deg2rad)

    x = (lon - lon_ref) * deg2rad * R2
    y = (lat - lat_ref) * deg2rad * R

    if rot is not 0:
        SINROT = np.sin(rot * deg2rad)
        COSROT = np.cos(rot * deg2rad)
        rotmat = np.array([[COSROT, -SINROT], [SINROT, COSROT]]) 
        xy = np.array([x, y])
        xy = np.swapaxes(xy, 0,1)
        XY = rotmat.dot(xy)
        x = XY[:,0]
        y = XY[:,1]

    if len(x) == 1:
        x = float(x)
        y = float(y)

    return x, y