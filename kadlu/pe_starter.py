import numpy as np
from enum import Enum
from kadlu.utils import get_member


class PEStarterMethod(Enum):
    GAUSSIAN = 1
    GREENE = 2
    THOMSON = 3


class PEStarter():

    def __init__(self, method='THOMSON', aperture=90):

        self.method = get_member(PEStarterMethod, method)
        self.aperture = aperture

    def eval(self, k0,  zs, Z, kz):
        """ Evaluate PE starter
            
            Args:
                k0: float
                    Reference wavenumber
                zs: float
                    Source depth
                Z: 2d numpy array
                    Depth at each grid point
                kz: 1d numpy array
                    ???
        """

        if self.method is PEStarterMethod.GAUSSIAN:
            psi = self._gaussian(k0, zs, Z)
        elif self.method is PEStarterMethod.GREENE:
            psi = self._greene(k0, zs, Z)
        elif self.method is PEStarterMethod.THOMSON:
            psi = self._thomson(k0, zs, Z, kz)
            
        return psi

    def _gaussian(self, k0, zs, Z):
        psi = np.sqrt(k0) * np.exp( -0.5*k0**2 *(Z-zs)**2 )
        psi = psi - ( np.sqrt(k0) * np.exp( -0.5*k0**2 *(Z+zs)**2 ))
        psi = np.fft.fft(psi)
        return psi

    def _greene(self, k0, zs, Z):
        a = 1.4467
        b = .04201
        c = 3.0512
        psi = np.sqrt(k0) * (a - b * k0**2 * (Z - zs)**2) * np.exp(-(k0**2 * (Z - zs)**2) / c )
        psi = psi - (np.sqrt(k0) * (a - b * k0**2 * (Z + zs)**2) * np.exp(-(k0**2 * (Z + zs)**2) / c ))
        psi = np.fft.fft(psi)
        return psi

    def _thomson(self, k0, zs, z, kz):
        psi = np.exp(-1j * np.pi / 4.) * 2 * np.sqrt(2 * np.pi) * np.sin(kz * zs) / np.sqrt(np.sqrt(k0**2 - kz**2))
        # normalize the starter
        psi = psi / (Z[1] - Z[0])
        psi[Z.shape[0] / 2 + 1, :] = 0 
        # taper the spectrum to obtain desired angle using Turkey window
        kcut1 = k0 * np.sin(self.aperture / 180 * np.pi) 
        kcut0 = k0 * np.sin((self.aperture - 1.5) / 180 * np.pi)
        W = 0.5 * (1 + np.cos(np.pi / (kcut1 - kcut0) * (np.abs(kz) - kcut0)))
        W[np.abs(kz) >= kcut1] = 0
        W[np.abs(kz) <= kcut0] = 1
        psi = psi * W
        psi[np.abs(kz) >= kcut1] = 0
        return psi

