""" Ocean module within the kadlu package

    Authors: Oliver Kirsebom
    contact: oliver.kirsebom@dal.ca
    Organization: MERIDIAN-Intitute for Big Data Analytics
    Team: Acoustic data Analytics, Dalhousie University
    Project: packages/kadlu
             Project goal: Tools for underwater soundscape modeling
     
    License:

"""

import numpy as np
from kadlu.utils import LatLon
from kadlu.geospatial.data_sources.chs import Chs
from kadlu.geospatial.data_sources.hycom import Hycom
from kadlu.geospatial.data_sources.era5 import Era5
from kadlu.geospatial.data_sources.wwiii import Wwiii
import kadlu.geospatial.data_sources.gebco as gebco 
from kadlu.geospatial.interpolation import Interpolator2D, Interpolator3D, Uniform2D, Uniform3D
from kadlu.geospatial.data_sources.fetch_util import index
from datetime import datetime


class Ocean():
    """ Class for handling ocean data requests.

        TODO: Implement loading of temp, salinity and wave data.

        Args:
            lats: 1d numpy array
                Latitude coordinates used for interpolation. 
            lons: 1d numpy array
                Longitude coordinates used for interpolation.
            depths: 1d numpy array
                Depth coordinates used for interpolation.

        Attributes: 

    """
    def __init__(self, bathy=None, temp=None, salinity=None, wave=None, wave_var=None, water_density=1.0):
        self.water_density = water_density
        self.data_source = {'bathy': bathy, 'temp': temp, 'salinity': salinity, 'wave': wave}

        self.bathy_data = None
        self.bathy_interp = None
        self.temp_data = None
        self.temp_interp = None
        self.salinity_data = None
        self.salinity_interp = None
        self.wave_data = None
        self.wave_interp = None

        self.set_origin(0,0)

        self.SW = LatLon(-90, -180)
        self.NE = LatLon(90, 180)

        if not isinstance(bathy, str):
            self.load_bathy()

        if not isinstance(temp, str):
            self.load_temp()

        if not isinstance(salinity, str):
            self.load_salinity()

        if not isinstance(wave, str):
            self.load_wave()


    def load(self, south=-90, north=90, west=-180, east=180, 
            start=datetime(2019, 1, 1), end=datetime(2019, 1, 1, 1), 
            top=0, bottom=5000):

        # save south-west and north-east corners as class attributes
        self.SW = LatLon(south, west)
        self.NE = LatLon(north, east)

        # origo of x-y coordinate system
        lat_ref = 0.5 * (south + north)
        lon_ref = 0.5 * (west + east)
        self.set_origin(lat_ref, lon_ref)

        # load data and initialize interpolation tables
        self.load_bathy(south, north, west, east)
        self.load_temp(south, north, west, east, start, end, top, bottom)
        self.load_salinity(south, north, west, east, start, end, top, bottom)
        self.load_wave(south, north, west, east, start, end)


    def set_origin(self, lat_ref, lon_ref):

        self.origin = LatLon(lat_ref, lon_ref)

        for interp in [
                self.bathy_interp, 
                self.temp_interp, 
                self.salinity_interp, 
                self.wave_interp
                ]:
            if interp is not None: interp.origin = self.origin

        """
        if self.bathy_interp is not None:
            self.bathy_interp.origin = self.origin
        if self.temp_interp is not None:
            self.temp_interp.origin = self.origin
        if self.salinity_interp is not None:
            self.salinity_interp.origin = self.origin
        if self.wave_interp is not None:
            self.wave_interp.origin = self.origin
        """


    def load_bathy(self, south=None, north=None, west=None, east=None, storage=None):
        bathy = self.data_source['bathy']

        if bathy is None:
            self.bathy_data = None
            self.bathy_interp = None

        elif isinstance(bathy, str):
    
            if bathy == "CHS":
                # load data and make coordinate grid for interpolation
                self.bathy_data = Chs().load_bathymetry(south, north, west, east)
                num_lats = int(np.ceil((north - south) / 0.001)) + 1
                num_lons = int(np.ceil((east - west) / 0.001)) + 1
                lats = np.linspace(south, north, num=num_lats)
                lons = np.linspace(west, east, num=num_lons)

                # interpolate
                self.bathy_interp = Interpolator2D(
                        values=self.bathy_data[0],
                        lats=self.bathy_data[1], lons=self.bathy_data[2],
                        origin=self.origin, method_irreg='regular', 
                        lats_reg=lats, lons_reg=lons)       

            elif bathy == "GEBCO":
                self.bathy_data = gebco.load(south, north, west, east)
                self.bathy_interp = Interpolator2D(
                        values=self.bathy_data[0],
                        lats=self.bathy_data[1], lons=self.bathy_data[2],
                        origin=self.origin)       

            else: 
                print('Error: Unknown bathymetry source {0}.'.format(bathy))
                exit(1)

        elif isinstance(bathy, (np.ndarray, tuple)):
        #elif isinstance(bathy, Iterable) and (bathy[-1][0] == 'chs' or bathy[-1][0] == 'gebco'):
            self.bathy_data = bathy
            self.bathy_interp = Interpolator2D(
                    values=self.bathy_data[0],
                    lats=self.bathy_data[1], lons=self.bathy_data[2], 
                    origin=self.origin)

        else:
            self.bathy_data = bathy
            self.bathy_interp = Uniform2D(bathy)


    def load_temp(self, south=None, north=None, west=None, east=None,
            start=None, end=None, top=0, bottom=5000):
        temp = self.data_source['temp']

        if temp is None:
            self.temp_data = None
            self.temp_interp = None

        elif isinstance(temp, str):
            if temp == "HYCOM":
                # TODO:
                # enable hycom nearest time slice loading
                self.temp_data = Hycom().load_temp(
                        south, north, west, east,
                        start, end, top, bottom)
                """ run this code interactively for testing load_temp and interpolation

                south, west = 46, -60
                north, east = 48, -58
                top, bottom = 0, 5000
                start,  end = datetime(2015, 1, 10), datetime(2015, 1, 10, 12)

                from kadlu.geospatial.data_sources import hycom
                from importlib import reload
                reload(hycom)
                self = Ocean()

                hycom.Hycom().fetch_temp(
                        south=south, north=north, west=west, east=east,
                        start=start, end=end, top=top, bottom=bottom
                        )
                self.temp_data = hycom.Hycom().load_temp(
                        south=south, north=north, west=west, east=east,
                        start=start, end=end, top=top, bottom=bottom
                        )
                """

                # build split index, reduce 4th dimension to 3D avg
                splidx = np.nonzero(self.temp_data[3][1:] > self.temp_data[3][:-1])[0] + 1
                if len(splidx) > 0:  # if fourth dimension exists
                    """
                    assert (self.temp_data[1][splidx[0]:splidx[1]] == self.temp_data[1][splidx[1]:splidx[2]]).all()
                    assert (self.temp_data[2][splidx[0]:splidx[1]] == self.temp_data[2][splidx[1]:splidx[2]]).all()
                    """
                    splidx = np.append(splidx, len(self.temp_data[3]))
                    splarr = np.array([self.temp_data[0][splidx[d]:splidx[d+1]] for d in range(0, len(splidx)-1)])
                    vals = (reduce(np.add, splarr) / len(splarr)).astype(float)
                    x = self.temp_data[2] [splidx[0] : splidx[1]] .astype(float)
                    y = self.temp_data[1] [splidx[0] : splidx[1]] .astype(float)
                    z = self.temp_data[4] [splidx[0] : splidx[1]] .astype(float)
                    warnings.warn("query data has been averaged across the time dimension for 3D interpolation")
                else:   # if data is already 3 dimensional, just leave as is
                    vals = self.temp_data[0].astype(float)
                    x = self.temp_data[2].astype(float)
                    y = self.temp_data[1].astype(float)
                    z = self.temp_data[4].astype(float)

                # get size of each dimension and create lat/lon/depth grid arrays
                # TODO:
                # fix bug here on xdim, ydim
                xdim = reduce(np.subtract, np.nonzero(x[1:] != x[:-1])[0][-1:-3:-1])
                ydim = reduce(np.subtract, np.nonzero(y[1:] < y[:-1])[0][-1:-3:-1])
                if min(z) == max(z): zdim = 1
                else: 
                    zdim = reduce(np.subtract, np.nonzero(z[1:] > z[:-1])[0][-1:-3:-1])
                xgrid = x[0::xdim]
                ygrid = y[0:ydim]
                zgrid = z[0::zdim]

                # rows of averaged datapoints with lat/lon/depth coords
                rowdata = np.array((vals, y, x, z)).T

                # create empty grid with default values, then populate with data
                """
                gridspace = np.full((
                            len(np.unique(self.temp_data.T[:,1])),  # lat
                            len(np.unique(self.temp_data.T[:,2])),  # lon
                            len(np.unique(self.temp_data.T[:,4]))   # depth
                        ), 
                        fill_value=-30000
                    )
                """

                gridspace = np.full((
                            len(np.unique(rowdata[:,1])),  # lat
                            len(np.unique(rowdata[:,2])),  # lon
                            len(np.unique(rowdata[:,3]))   # depth
                        ), 
                        fill_value=-30000
                    )

                for row in rowdata:
                    y_ix = index(row[1], ygrid)
                    x_ix = index(row[2], xgrid)
                    z_ix = index(row[3], zgrid)
                    gridspace[y_ix, x_ix, z_ix] = row[0]

                self.temp_interp = Interpolator3D(
                        values=gridspace,
                        lats=ygrid, lons=xgrid,
                        depths=zgrid, origin=self.origin,
                        method='linear')

            else: 
                print('Error: Unknown temperature source {0}.'.format(temp))
                exit(1)

        #elif isinstance(temp, tuple):
        #elif isinstance(temp, Iterable) and temp[-1][0] == 'hycom':
        elif isinstance(temp, (np.ndarray, tuple)):
            self.temp_data = temp
            self.temp_interp = Interpolator3D(
                    values=self.temp_data[0],
                    lats=self.temp_data[1],     lons=self.temp_data[2],
                    depths=self.temp_data[4],   origin=self.origin,
                    method='linear')

        else:
            self.temp_data = temp
            self.temp_interp = Uniform3D(temp)


    def load_salinity(self, south=None, north=None, west=None, east=None, 
            start=None, end=None, top=0, bottom=5000):
        salinity = self.data_source['salinity']

        if salinity is None:
            self.salinity_data = None
            self.salinity_interp = None

        elif isinstance(salinity, str):
            if salinity == "HYCOM":
                self.salinity_data = Hycom().load_salinity(
                        south, north, west, east,
                        start, end, top, bottom)
                """
                num_lats = int(np.ceil((north - south) / 0.001)) + 1
                num_lons = int(np.ceil((east - west) / 0.001)) + 1
                lats = np.linspace(south, north, num=num_lats)
                lons = np.linspace(west, east, num=num_lons)
                """

                # unflatten array into cube

                self.salinity_interp = Interpolator3D(
                        values=self.salinity_data[0], 
                        lats=self.salinity_data[1], lons=self.salinity_data[2], 
                        depths=self.salinity_data[4], origin=self.origin,
                        method='linear')       

            else: 
                print('Error: Unknown salinity source {0}.'.format(salinity))
                exit(1)

        elif isinstance(salinity, (np.ndarray, tuple)):
        #elif isinstance(salinity, Iterable) and salinity[-1][0] == 'hycom':
            self.salinity_data = salinity
            self.salinity_interp = Interpolator3D(
                    values=self.salinity_data[0], 
                    lats=self.salinity_data[1], lons=self.salinity_data[2], 
                    depths=self.salinity_data[4], origin=self.origin,
                    method='linear')

        else:
            self.salinity_data = salinity
            self.salinity_interp = Uniform3D(salinity)


    def load_wave(self, south=None, north=None, west=None, east=None, start=None, end=None):

        wave = self.data_source['wave']

        if wave is None:
            self.wave_data = None
            self.wave_interp = None

        elif isinstance(wave, str):
            if wave == "ERA5":
                # load data and make coordinate grid for interpolation
                self.wave_data = Era5().load_windwaveswellheight(south, north, west, east, start, end)
                num_lats = int(np.ceil((north - south) / 0.001)) + 1
                num_lons = int(np.ceil((east - west) / 0.001)) + 1
                lats = np.linspace(south, north, num=num_lats)
                lons = np.linspace(west, east, num=num_lons)

                # interpolate
                self.wave_interp = Interpolator2D(
                        values=self.wave_data[0],
                        lats=self.wave_data[1], lons=self.wave_data[2],
                        origin=self.origin, method_irreg='regular',
                        lats_reg=lats, lons_reg=lons)       

            elif wave == "RDWPS":
                # load data and make coordinate grid for interpolation
                self.wave_data = Rdwps().load_windwaveswellheight(south, north, west, east, start, end)
                num_lats = int(np.ceil((north - south) / 0.001)) + 1
                num_lons = int(np.ceil((east - west) / 0.001)) + 1
                lats = np.linspace(south, north, num=num_lats)
                lons = np.linspace(west, east, num=num_lons)

                # interpolate
                self.wave_interp = Interpolator2D(
                        values=self.wave_data[0],
                        lats=self.wave_data[1], lons=self.wave_data[2],
                        origin=self.origin, method_irreg='regular',
                        lats_reg=lats, lons_reg=lons)       

            elif wave == "WWIII":
                # load data and make coordinate grid for interpolation
                self.wave_data = Wwiii().load_windwaveheight(south, north, west, east, start, end)
                num_lats = int(np.ceil((north - south) / 0.001)) + 1
                num_lons = int(np.ceil((east - west) / 0.001)) + 1
                lats = np.linspace(south, north, num=num_lats)
                lons = np.linspace(west, east, num=num_lons)

                # interpolate
                self.wave_interp = Interpolator2D(
                        values=self.wave_data[0],
                        lats=self.wave_data[1], lons=self.wave_data[2],
                        origin=self.origin, method_irreg='regular',
                        lats_reg=lats, lons_reg=lons)       

            else: 
                print('Error: Unknown wave source {0}.'.format(wave))
                exit(1)

        elif isinstance(wave, (np.ndarray, tuple)):
        #elif isinstance(wave, Iterable) and len(wave) == 3:
            self.wave_data = wave
            self.wave_interp = Interpolator2D(
                    values=self.wave_data[0], 
                    lats=self.wave_data[1], lons=self.wave_data[2], 
                    origin=self.origin)       

        else:
            self.wave_data = wave
            self.wave_interp = Uniform2D(wave)


    def bathy(self, x=None, y=None, grid=False, geometry='planar'):
        """ Evaluate interpolated bathymetry in spherical (lat-lon) or  
            planar (x-y) geometry.

            x and y can be floats or arrays.

            If grid is set to False, the bathymetry will be evaluated at 
            the positions (x_i, y_i), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_N). Note that in this case, x and y must have 
            the same length.

            If grid is set to True, the bathymetry will be evaluated at 
            all combinations (x_i, y_j), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_M). Note that in this case, the lengths of x 
            and y do not have to be the same.

            If x and y are not specified, the method returns the underlying 
            bathymetric data on which the interpolation is performed, either 
            as a (bathy,lat,lon) tuple, or as a float if the bathymetry is 
            the same everywhere.

            Args: 
                x: float or array
                    x-coordinate(s) or longitude(s)
                y: float or array
                    y-coordinate(s) or latitude(s)
                grid: bool
                    Specify how to combine elements of x and y. If x and y have different
                    lengths, specifying grid has no effect as it is automatically set to True.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'

            Returns:
                z: Interpolated bathymetry values
        """
        assert self.bathy_data is not None, "Bathymetric data have not been loaded" 

        if x is None and y is None:
            z = self.bathy_data
        
        else:
            if geometry == 'planar':
                z = self.bathy_interp.eval_xy(x=x, y=y, grid=grid)

            elif geometry == 'spherical':
                z = self.bathy_interp.eval_ll(lat=y, lon=x, grid=grid)

        return z


    def bathy_gradient(self, x=None, y=None, axis='x', grid=False, geometry='planar'): 
        """ Evaluate interpolated bathymetry gradient in spherical (lat-lon) or  
            planar (x-y) geometry, along the specified axis.

            x and y can be floats or arrays.

            If grid is set to False, the bathymetry will be evaluated at 
            the positions (x_i, y_i), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_N). Note that in this case, x and y must have 
            the same length.

            If grid is set to True, the bathymetry will be evaluated at 
            all combinations (x_i, y_j), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_M). Note that in this case, the lengths of x 
            and y do not have to be the same.

            Args: 
                x: float or array
                   x-coordinate(s) or longitude(s)
                y: float or array
                   y-coordinate(s) or latitude(s)
                axis: str
                    Axis along which gradient is computed. Can be either 'x' (default) or 'y'
                grid: bool
                   Specify how to combine elements of x and y.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'

            Returns:
                grad: Interpolated bathymetry gradient values
        """
        assert self.bathy_data is not None, "Bathymetric data have not been loaded" 

        deriv_order = [(axis=='x'), (axis!='x')]

        if geometry == 'planar':                
            grad = self.bathy_interp.eval_xy(x=x, y=y, grid=grid, x_deriv_order=deriv_order[0], y_deriv_order=deriv_order[1])

        elif geometry == 'spherical':
            grad = self.bathy_interp.eval_ll(lat=y, lon=x, grid=grid, lat_deriv_order=deriv_order[1], lon_deriv_order=deriv_order[0])

        return grad


    def temp(self, x=None, y=None, z=None, grid=False, geometry='planar'):
        """ Evaluate interpolated temperature in spherical (lat-lon) or  
            planar (x-y) geometry.

            x,y,z can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the positions (x_i, y_i, z_i), where x=(x_1,...,x_N),  
            y=(y_1,...,y_N), and z=(z_1,...,z_N). Note that in this case, 
            x,y,z must have the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (x_i, y_j, z_k), where x=(x_1,...,x_N), 
            y=(y_1,...,y_M), and z=(z_1,...,z_K). Note that in this case, the 
            lengths of x,y,z do not have to be the same.

            If x,y,z are not specified, the method returns the underlying 
            temperature data on which the interpolation is performed, either 
            as a (temp,lat,lon,z) tuple, or as a float if the temperature is 
            the same everywhere.

            Args: 
                x: float or array
                    x-coordinate(s) or longitude(s)
                y: float or array
                    y-coordinate(s) or latitude(s)
                z: float or array
                    depth(s)
                grid: bool
                   Specify how to combine elements of x,y,z.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'

            Returns:
                t: Interpolated temperature values
        """
        assert self.temp_data is not None, "Temperature data have not been loaded" 

        if x is None and y is None and z is None:
            t = self.temp_data

        else:
            if geometry == 'planar':
                t = self.temp_interp.eval_xy(x=x, y=y, z=z, grid=grid)

            elif geometry == 'spherical':
                t = self.temp_interp.eval_ll(lat=y, lon=x, z=z, grid=grid)

        return t


    def salinity(self, x=None, y=None, z=None, grid=False, geometry='planar'):
        """ Evaluate interpolated salinity in spherical (lat-lon) or  
            planar (x-y) geometry.

            x,y,z can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the positions (x_i, y_i, z_i), where x=(x_1,...,x_N),  
            y=(y_1,...,y_N), and z=(z_1,...,z_N). Note that in this case, 
            x,y,z must have the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (x_i, y_j, z_k), where x=(x_1,...,x_N), 
            y=(y_1,...,y_M), and z=(z_1,...,z_K). Note that in this case, the 
            lengths of x,y,z do not have to be the same.

            If x,y,z are not specified, the method returns the underlying 
            salinity data on which the interpolation is performed, either 
            as a (salinity,lat,lon,z) tuple, or as a float if the salinity is 
            the same everywhere.

            Args: 
                x: float or array
                    x-coordinate(s) or longitude(s)
                y: float or array
                    y-coordinate(s) or latitude(s)
                z: float or array
                    depth(s)
                grid: bool
                   Specify how to combine elements of x,y,z.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'

            Returns:
                s: Interpolated salinity values
        """
        assert self.salinity_data is not None, "Salinity data have not been loaded" 

        if x is None and y is None and z is None:
            s = self.salinity_data

        if geometry == 'planar':
            s = self.salinity_interp.eval_xy(x=x, y=y, z=z, grid=grid)

        elif geometry == 'spherical':
            s = self.salinity_interp.eval_ll(lat=y, lon=x, z=z, grid=grid)

        return s


    def wave(self, x, y, grid=False, geometry='planar'):
        """ Evaluate interpolated wave data in spherical (lat-lon) or  
            planar (x-y) geometry.

            x and y can be floats or arrays.

            If grid is set to False, the interpolation will be evaluated at 
            the positions (x_i, y_i), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_N). Note that in this case, x and y must have 
            the same length.

            If grid is set to True, the interpolation will be evaluated at 
            all combinations (x_i, y_j), where x=(x_1,...,x_N) and 
            y=(y_1,...,y_M). Note that in this case, the lengths of x 
            and y do not have to be the same.

            If x and y are not specified, the method returns the underlying 
            wave data on which the interpolation is performed, either 
            as a (wave,lat,lon) tuple, or as a float if the wave data is the same 
            everywhere.

            Args: 
                x: float or array
                    x-coordinate(s) or longitude(s)
                y: float or array
                    y-coordinate(s) or latitude(s)
                grid: bool
                    Specify how to combine elements of x and y. If x and y have different
                    lengths, specifying grid has no effect as it is automatically set to True.
                geometry: str
                    Can be either 'planar' (default) or 'spherical'

            Returns:
                w: Interpolated wave data
        """
        assert self.wave_data is not None, "Wave data have not been loaded" 

        if x is None and y is None:
            w = self.wave_data

        else:
            if geometry == 'planar':
                w = self.wave_interp.eval_xy(x=x, y=y, grid=grid)

            elif geometry == 'spherical':
                w = self.wave_interp.eval_ll(lat=y, lon=x, grid=grid)

        return w


