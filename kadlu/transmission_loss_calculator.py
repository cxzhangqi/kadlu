import numpy as np
from numpy.lib import scimath
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

        # reference wavelength and wavenumber
        lambda0 = self.c0 / freq        
        k0 = 2 * np.pi * freq / self.c0  

        smoothing_length_rho = self.c0 / freq / 4
        smoothing_length_ssp = np.finfo(float).eps  # machine epsilon

        # radial step size
        dr = lambda0 / 2

        # vertical range
        ThinknessOfArtificialAbsorpLayer_ratio_z = 6
        vertical_range = 2 * 12e3 / ThinknessOfArtificialAbsorpLayer_ratio_z * (ThinknessOfArtificialAbsorpLayer_ratio_z + 1)

        # create grids
        Z, x, y, z, kz, nx, ny, nz = self._create_grids(radial_step=dr, radial_range=50e3,\
                azimuthal_step=1/180*np.pi, azimuthal_range=2*np.pi,\
                vertical_step=10, vertical_range=vertical_range)

        print('\nnx,ny,nz: ', nx, ny, nz)
        print('x.shape: ', x.shape)
        print('y.shape: ', y.shape)
        print('z.shape: ', z.shape)

        # PE starter
        starter = PEStarter(method='THOMSON', aperture=88)

        # boudary conditions on z
        ArtificialAbsorpCoeff =  1. / np.log10(np.e) / lambda0 * 2 / k0

        ThinknessOfArtificialAbsorpLayer = np.max(np.abs(z)) / ThinknessOfArtificialAbsorpLayer_ratio_z
        D = ThinknessOfArtificialAbsorpLayer / 3
        atten0 = 1j * ArtificialAbsorpCoeff * np.exp(-(np.abs(Z[:]) - np.max(np.abs(z)))**2 / D**2)

        # initial field
        psi = starter.eval(k0=k0, zs=zs, Z=Z, kz=kz) * np.ones(shape=(1,ny))

        print('psi.shape: ', psi.shape)

        # initial Nx2D free propagator
        fr_half, fr_full = self._free_propagator(dr, k0, kz, ny)

        print('fr_half:', fr_half.shape)
        print('fr_full:', fr_full.shape)


    def _create_grids(self, radial_step, radial_range,\
            azimuthal_step, azimuthal_range, vertical_step, vertical_range):
        """ Create grids for PE solver
            
            Args:
                radial_step: float
                    Radial step size (in meters)
                radial_range: float
                    Radial range (in meters)
                azimuthal_step: float
                    Angular step size (in radians)
                azimuthal_range: float
                    Angular domain (in radians)
                vertical_step: float
                    Vertical step size (in meters)
                vertical_range: float
                    Vertical range (in meters)

            Returns:

        """
        dx = radial_step
        xmax = radial_range
        dy = azimuthal_step
        ymax = azimuthal_range
        dz = vertical_step
        zmax = vertical_range

        # ------ x (radial) ------ 
        nx = round(xmax / dx) # number of radial bins
        x = np.arange(nx+1, dtype=float)
        x *= dx

        # ------ y (azimuthal) ------ 
        ny = int(np.ceil(ymax / dy)) # number of angular bins
        if ny%2==1: ny = ny + 1 # ensure even number of angular bins
        y_pos = np.arange(ny/2, dtype=float)
        y_neg = np.arange(-ny/2-1, step=-1, dtype=float)
        y = np.concatenate((y_pos, y_neg)) 
        y *= dy

        # ------ z (vertical) ------ 
        nz = zmax / dz # number of vertical bins
        nz = round(nz / 2) * 2  # ensure even number of vertical bins
        Lz = nz * dz  
        z_pos = np.arange(nz/2, dtype=float)
        z_neg = np.arange(-nz/2-1, step=-1, dtype=float)
        z = np.concatenate((z_pos, z_neg))
        z *= dz

        # ------ kz (wavenumber) ------
        kz = z * 2 * np.pi / (Lz * dz)


        Y, Z = np.meshgrid(y, z)

        return Z, x, y, z, kz, nx, ny, nz

    def _free_propagator(self, dr, k0, kz, ny):
        fr_half = np.exp(1j * dr / 2 * (scimath.sqrt(k0**2 - kz**2) - k0))
        fr_half = fr_half[:,np.newaxis]
        fr_half = fr_half * np.ones(shape=(1,ny))
        fr_full = np.exp(1j * dr * (scimath.sqrt(k0**2 - kz**2) - k0)) 
        fr_full = fr_full[:,np.newaxis]        
        fr_full = fr_full * np.ones(shape=(1,ny))
        return fr_half, fr_full
