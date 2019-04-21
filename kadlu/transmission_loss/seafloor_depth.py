import numpy as np
from kadlu.bathy_reader import BathyReader, LatLon
from kadlu.bathy_interpolator import BathyInterpolator

class SeafloorDepth():

    def __init__(self, bathy_interpolator=None, flat_bottom_depth=None):

#        assert

        self.interp = bathy_interpolator
        self.flat_bottom_depth = flat_bottom_depth

#            filename = 'kadlu/tests/assets/BathyData_Mariana_500kmx500km.mat'
#            reader = BathyReader(input=filename, lat_name='latgrat', lon_name='longrat', bathy_name='mat')
#            refloc = LatLon(11.2, 142.1)  # reference location at 11.2 deg N and 142.1 deg E
#            self.interp = BathyInterpolator(bathy_reader=reader, origin=refloc)


    def get_depth(self, x, y):
        # Input:
        #  x,y coordinates (cartesian)
        # returns:
        #  1) water depth in all (n) angular bins
        #  2) derivative of depth wrt angle in all angular bins

        if self.flat_bottom_depth:

            n = len(x)
            depth = self.flat_bottom_depth * np.ones(n)
            gradient = np.zeros(n)

        else:

            depth = self.interp.eval_xy(x=x, y=y)
            depth *= (-1.)

            gradient = self.interp.eval_xy(x=x, y=y, y_deriv_order=1)
            gradient *= (-1.)


        if False:
            sigma_x = 10000.
            sigma_y = 30000.
            xc = 0
            yc = 0
            exponent = (x-xc)**2 / (2*sigma_x**2) + (y-yc)**2 / (2*sigma_y**2)
            val = np.exp(-exponent)
            depth = val * 10000.
            gradient = depth * 0

###        print(depth[90], gradient[90], x[90], y[90])

        return depth, gradient