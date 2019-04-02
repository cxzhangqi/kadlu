import numpy as np
from kadlu.bathy_reader import BathyReader, LatLon
from kadlu.bathy_interpolator import BathyInterpolator

class SeafloorDepth():

    def __init__(self, n):

        self.n = n

        filename = 'kadlu/tests/assets/BathyData_Mariana_500kmx500km.mat'
        reader = BathyReader(input=filename, lat_name='latgrat', lon_name='longrat', bathy_name='mat')
        refloc = LatLon(11.2, 142.1)  # reference location at 9 deg N and 140 deg E
        self.interp = BathyInterpolator(bathy_reader=reader, origin=refloc)


    def get_depth(self, x, y):
        # Input:
        #  x,y coordinates (cartesian)
        # returns:
        #  1) water depth in all (n) angular bins
        #  2) derivative of depth wrt angle in all angular bins

#        x = -x

        depth = self.interp.eval_xy(x=x, y=y)
        depth *= (-1.)

        # flat seafloor with depth of 1000 m:
###        depth = 10000 * np.ones(self.n)
        gradient = np.zeros(self.n)

###        print(depth[0], x[0], y[0])

        return depth, gradient