import numpy as np


class RefractiveIndex():

    def __init__(self, m):

        self.m = m

    def get_nsq(self, x, YZ):
        # Input:
        #  x: radial distance
        #  YZ: YZ grid
        # returns:
        #  1) refractive index squared on YZ grid a distance=x

        nsq = np.ones(shape=YZ.shape)

        return nsq