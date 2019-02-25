import numpy as np
from enum import Enum


class PEStarterMethod(Enum):
    GAUSSIAN = 1
    GREENE = 2
    THOMSON = 3


def get_member(cls, member_name):
    for name, member in cls.__members__.items():
        if member_name == name:
            return member

    s = ", ".join(name for name, _ in cls.__members__.items())
    raise ValueError("Unknown value \'{0}\'. Select between: {1}".format(member_name, s))


class PEStarter():

    def __init__(self, method, aperture):

        self.method = get_member(PEStarterMethod, method)
        self.aperture = aperture

    def eval(self, k0, kz, Z, zs):

        # Z is an array containing the z coordinates 
        # of the grid used by the PE solver
        # see l. 64 of propNx2DWAPE.m

        if self.method is PEStarterMethod.GAUSSIAN:
            psi = self._gaussian(k0, Z, zs)
        elif self.method is PEStarterMethod.GREENE:
            psi = self._greene(k0, Z, zs)
        elif self.method is PEStarterMethod.THOMSON:
            psi = self._thomson(k0, kz, Z, zs)
            
        return psi

    def _gaussian(self, k0, Z, zs):
        psi = np.sqrt(k0) * np.exp( -0.5*k0**2 *(Z-zs)**2 )
        psi = psi - ( np.sqrt(k0) * np.exp( -0.5*k0**2 *(Z+zs)**2 ))
        psi = np.fft.fft(psi)
        return psi

    def _greene(self, k0, Z, zs):
        a = 1.4467
        b = .04201
        c = 3.0512
        psi = np.sqrt(k0) * (a - b * k0**2 * (Z - zs)**2) * np.exp(-(k0**2 * (Z - zs)**2) / c )
        psi = psi - (np.sqrt(k0) * (a - b * k0**2 * (Z + zs)**2) * np.exp(-(k0**2 * (Z + zs)**2) / c ))
        psi = np.fft.fft(psi)
        return psi

    def _thomson(self, k0, kz, Z, zs):
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

