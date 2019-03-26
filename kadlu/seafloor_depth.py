import numpy as np


class SeafloorDepth():

    def __init__(self, n):

        self.n = n

    def get_depth(self, x, y):

        return np.ones(self.n), np.ones(self.n)