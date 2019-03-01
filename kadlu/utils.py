import numpy as np
import os

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

def torad(lat, lon):
    """ Convert latitute and longitude values from degrees to radians.

        The method expects the latitude to be in the range (-90,90) and
        the longitude to be in the range (-180,180).

        The output latitude is in the range (0,pi) and the output 
        longitude is in the range (-pi,pi).

        Args: 
            lat: float or array
                latitude(s) in degrees from -90 to +90.
            lon: float or array
                longitude(s) in degrees from -180 to +180.

        Returns:
            lat_rad: float or array
                latitude(s) in radians from 0 to pi.
            lon_rad: float or array
                longitude(s) in radians from -pi to +pi.
    """
    lat_rad = (lat + 90) * deg2rad
    lon_rad = lon * deg2rad
    return lat_rad, lon_rad


def get_files(path, substr, fullpath=True, subdirs=False):
    """ Find all files in the specified directory containing the specified substring in their file name

        Args:
            path: str
                Directory path
            substr: str
                Substring contained in file name
            fullpath: bool
                Return full path to each file or just the file name 
            subdirs: bool
                Also search all subdirectories

        Returns:
            files: list (str)
                Alphabetically sorted list of file names
    """
    # find all files
    allfiles = list()
    if not subdirs:
        f = os.listdir(path)
        for fil in f:
            if fullpath:
                x = path
                if path[-1] is not '/':
                    x += '/'
                allfiles.append(os.path.join(x, fil))
            else:
                allfiles.append(fil)
    else:
        for r, _, f in os.walk(path):
            for fil in f:
                if fullpath:
                    allfiles.append(os.path.join(r, fil))
                else:
                    allfiles.append(fil)

    # select those that contain specified substring
    files = list()
    for f in allfiles:
        if substr in f:
            files.append(f)

    # sort alphabetically
    files.sort()

    return files


def get_member(cls, member_name):
    for name, member in cls.__members__.items():
        if member_name == name:
            return member

    s = ", ".join(name for name, _ in cls.__members__.items())
    raise ValueError("Unknown value \'{0}\'. Select between: {1}".format(member_name, s))