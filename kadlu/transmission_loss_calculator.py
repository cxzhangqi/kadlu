import numpy as np
from kadlu.pe_starter import PEStarter


class TransmissionLossCalculator():

    def __init__(self, bathymetry, sound_speed):

        self.bathy = bathymetry
        self.sound_speed = sound_speed

        self.c0 = 1500      # reference sound speed in m/s
        self.rhow = 1       # water density
        self.cb = 1700      # homogeneous bottom sound speed
        self.bloss = 0.5    # homogeneous bottom attenuation db/lambda
        self.rhob = 1.5     # homogeneous bottom density

        self.nsq = 1 # = (self.c0 / self.sound_speed)^2  refractive index squared


    def run(self, frequency, source_depth):

        # frequency in Hz
        freq = frequency  
        
        # source position in meters
        xs = 0
        ys = 0
        zs = source_depth

        # reference wavelenght and wavenumber
        lambda0 = self.c0 / freq        
        k0 = 2 * np.pi * freq / self.c0  

        # create grids
        Z = self._create_grids(lambda0)

        # PE starter
        starter = PEStarter(method='THOMSON', aperture=88)


    def _create_grids(self, lambda0):

        # ------ r (x) ------ 
        steplength = lambda0 / 2         # step size in meters
        rmax = 50e3                      # radial range in meters
        numstep = round(rmax / steplength) # number of radial marching steps
        dx = steplength
        x = np.arange(numstep+1, dtype=float)
        x *= dx

        # ------ theta (y) ------ 
        Ltheta = 2 * np.pi # angular domain in radians
        dtheta = 1/180 * np.pi # angular bin size in radians
        ntheta = int(np.ceil(Ltheta/dtheta)) # number of angular bins
        if ntheta%2==1: ntheta = ntheta+1 # ensure even number of angular bins
        dy = dtheta
        ny = ntheta
        Ly = dy * ny   # rough estimate of angular aperature
        y_pos = np.arange(ny/2, dtype=float)
        y_neg = np.arange(-ny/2-1, step=-1, dtype=float)
        y = np.concatenate((y_pos, y_neg)) * dy

        # ------ z ------ 
        dz = 10 # vertical bin size in meters
        ThinknessOfArtificialAbsorpLayer_ratio_z = 6 # ?
        nz = 2 * 12e3 / ThinknessOfArtificialAbsorpLayer_ratio_z * (ThinknessOfArtificialAbsorpLayer_ratio_z + 1) / dz # number of vertical bins
        nz = round(nz/2) * 2  # ensure even number of vertical bins
        Lz = nz * dz  
        z_pos = np.arange(nz/2, dtype=float)
        z_neg = np.arange(-nz/2-1, step=-1, dtype=float)
        z = np.concatenate((z_pos, z_neg)) * dz
        kz = z * 2 * np.pi / (Lz * dz)


        Y, Z = np.meshgrid(y, z)

        return Z