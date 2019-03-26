import numpy as np


class RefractiveIndex():

    def __init__(self, m):

        self.m = m

    def get_nsq(self, x, YZ):

        return np.ones(shape=YZ.shape)