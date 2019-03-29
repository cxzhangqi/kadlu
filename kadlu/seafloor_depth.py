import numpy as np


class SeafloorDepth():

    def __init__(self, n):

        self.n = n

    def get_depth(self, x, y):
        # Input:
        #  x,y coordinates (cartesian)
        # returns:
        #  1) water depth in all (n) angular bins
        #  2) derivative of depth wrt angle in all angular bins

        # flat seafloor with depth of 1000 m:
        depth = 1000 * np.ones(self.n)
        gradient = np.zeros(self.n)

        return depth, gradient