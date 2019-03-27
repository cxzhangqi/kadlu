import numpy as np

class ModelOutput():

    def __init__(self, Y, Z):

        # which points to be output?
        self.yso = np.argwhere(Y[0,:] >= 0)[0]
        self.nzhalf = int(Z.shape[0] / 2)

    def get_output(self, dista):

        return 0