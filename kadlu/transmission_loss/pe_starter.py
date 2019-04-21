import numpy as np
from numpy.lib import scimath
from enum import Enum
from kadlu.utils import get_member


class PEStarterMethod(Enum):
    GAUSSIAN = 1
    GREENE = 2
    THOMSON = 3


class PEStarter():
    """ Computes starting field for Parabolic-Equation propagator.
        
        Args:
            ref_wavenumber: float or numpy array
                Reference wavenumber(s)
            grid: PEGrid
                Computational grid
            method: str
                Options are: GAUSSIAN, GREENE, THOMSON
            aperture: float
                Aperture in degrees
    """
    def __init__(self, ref_wavenumber, grid, method='THOMSON', aperture=90):

        self.method = get_member(PEStarterMethod, method)
        self.aperture = aperture
        self.k0 = ref_wavenumber
        self.grid = grid

    def eval(self, zs):
        """ Evaluate PE starter at specified source depth
            
            Args:
                zs: float
                    Source depth in meters
        """
        if self.method is PEStarterMethod.GAUSSIAN:
            psi = self._gaussian(zs)
        elif self.method is PEStarterMethod.GREENE:
            psi = self._greene(zs)
        elif self.method is PEStarterMethod.THOMSON:
            psi = self._thomson(zs)
            
        return psi

    def _gaussian(self, zs):
        k0 = self.k0
        Z = self.grid.Z
        # compute psi
        psi = np.sqrt(k0) * np.exp( -0.5*k0**2 *(Z-zs)**2 )
        psi = psi - ( np.sqrt(k0) * np.exp( -0.5*k0**2 *(Z+zs)**2 ))
        psi = np.fft.fft(psi)
        return psi

    def _greene(self, zs):
        a = 1.4467
        b = .04201
        c = 3.0512
        k0 = self.k0
        Z = self.grid.Z
        # compute psi        
        psi = np.sqrt(k0) * (a - b * k0**2 * (Z - zs)**2) * np.exp(-(k0**2 * (Z - zs)**2) / c )
        psi = psi - (np.sqrt(k0) * (a - b * k0**2 * (Z + zs)**2) * np.exp(-(k0**2 * (Z + zs)**2) / c ))
        psi = np.fft.fft(psi)
        return psi

    def _thomson(self, zs):
        k0 = self.k0
        kz = self.grid.kz
        dz = self.grid.dz
        Nz = self.grid.Nz
        # compute psi
        psi = np.exp(-1j * np.pi / 4.) * 2 * scimath.sqrt(2 * np.pi) * np.sin(kz * zs) / scimath.sqrt(scimath.sqrt(k0**2 - kz**2))
        # normalize the starter
        psi = psi / dz
        psi[int(Nz/2)] = 0 
        # taper the spectrum to obtain desired angle using Turkey window
        kcut1 = k0 * np.sin(self.aperture / 180 * np.pi) 
        kcut0 = k0 * np.sin((self.aperture - 1.5) / 180 * np.pi)
        W = 0.5 * (1 + np.cos(np.pi / (kcut1 - kcut0) * (np.abs(kz) - kcut0)))
        W[np.abs(kz) >= kcut1] = 0
        W[np.abs(kz) <= kcut0] = 1
        psi = psi * W
        psi[np.abs(kz) >= kcut1] = 0
        psi = psi[:, np.newaxis]
        return psi

