import numpy as np


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

        freq = frequency  # frequency in Hz
        
        # source position in meters
        xs = 0
        ys = 0
        zs = source_depth

        lambda0 = self.c0 / freq  # reference wavelength

        # --- set parameters for PE solver ---

        # r
        steplength = lambda0 / 2         # step size in meters
        rmax = 50e3                      # radial range in meters
        ndxout = 1                       # ?
        numstep = round(rmax/steplength) # number of radial marching steps

        # theta
        Ltheta = 2 * np.pi # angular domain in radians
        dtheta = 1/180 * np.pi # angular bin size in radians
        ntheta = int(np.ceil(Ltheta/dtheta)) # number of angular bins
        if ntheta%2==1: ntheta = ntheta+1 # ensure even number of angular bins
        ndy_3DSliceout=1 # ?

        # z
        dz = 10 # vertical bin size in meters
        ThinknessOfArtificialAbsorpLayer_ratio_z = 6 # ?
        nz = 2 * 12e3 / ThinknessOfArtificialAbsorpLayer_ratio_z * (ThinknessOfArtificialAbsorpLayer_ratio_z + 1) / dz # number of vertical bins
        nz = round(nz/2) * 2  # ensure even number of vertical bins
        ndz_3DSliceout = 1 # ?

        # PE starter
        pe_starter_type = 'Thomson''s'
        pe_starter_aperature = 88 # degress
