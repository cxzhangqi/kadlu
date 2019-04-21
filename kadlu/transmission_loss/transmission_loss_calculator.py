import numpy as np
from numpy.lib import scimath
from kadlu.transmission_loss.pe_starter import PEStarter
from kadlu.transmission_loss.pe_propagator import PEPropagator
from kadlu.transmission_loss.environment_input import EnvironmentInput
from kadlu.transmission_loss.env_input import EnvInput
import math


class TransmissionLossCalculator():

    def __init__(self, bathymetry, sound_speed, step_size=None, range=50e3,\
            angular_bin_size=1, vertical_bin_size=10, max_depth=12e3):

        self.bathy = bathymetry
        self.sound_speed = sound_speed

        self.starter_method = 'THOMSON'
        self.starter_aperture = 88

        self.c0 = 1500      # reference sound speed in m/s
        self.rhow = 1       # water density
        self.cb = 1700      # homogeneous bottom sound speed
        self.bloss = 0.5    # homogeneous bottom attenuation db/lambda
        self.rhob = 1.5     # homogeneous bottom density
        self.ndx_ChangeWD = 3 #1  how often to update bathymetry

        self.ndx_ChangeNSQ = math.inf  # how often to update water column (inf => range-independent SSP)

        self.nsq = 1 # = (self.c0 / self.sound_speed)^2  refractive index squared

        self.step_size = step_size
        self.range = range
        self.angular_bin_size = angular_bin_size
        self.vertical_bin_size = vertical_bin_size
        self.max_depth = max_depth


    def run(self, frequency, source_depth, vertical_slice=True, depths=[.1]):


        import time
        start = time.time()

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
        if self.step_size is None:
            dr = lambda0 / 2
        else:
            dr = self.step_size

        # azimuthal step size
        azimuthal_step = self.angular_bin_size / 180 * np.pi

        # vertical range
        ThinknessOfArtificialAbsorpLayer_ratio_z = 6
        vertical_range = 2 * self.max_depth / ThinknessOfArtificialAbsorpLayer_ratio_z * (ThinknessOfArtificialAbsorpLayer_ratio_z + 1)
        vertical_step = self.vertical_bin_size

        grid = PEGrid(radial_step=dr, radial_range=self.range,\
                azimuthal_step=azimuthal_step, azimuthal_range=2*np.pi,\
                vertical_step=vertical_step, vertical_range=vertical_range)

        print('\nnx,ny,nz: ', grid.Nr, grid.Nq, grid.Nz)
        print('x.shape: ', grid.r.shape)
        print('y.shape: ', grid.q.shape)
        print('z.shape: ', grid.z.shape)
        print('Y.shape: ', grid.Q.shape)
        print('Z.shape: ', grid.Z.shape)
        print('y: {0:.3f},{1:.3f}...{2:.3f},{3:.3f},{4:.3f}...,{5:.3f},{6:.3f}'.format(grid.q[0],grid.q[1],grid.q[int(grid.Nq/2)-1],grid.q[int(grid.Nq/2)],grid.q[int(grid.Nq/2)+1],grid.q[-2],grid.q[-1]))
        print('dx: {0:.1f} m'.format(grid.dr))

        # PE starter
        starter = PEStarter(ref_wavenumber=k0, grid=grid, method=self.starter_method, aperture=self.starter_aperture)

        # initial field
        psi = starter.eval(zs) * np.ones(shape=(1,grid.Nq))
        print('psi.shape: ', psi.shape)

        # module handling updates of environmental input
        env_input = EnvironmentInput(ref_wavenumber=k0, grid=grid, xs=xs, ys=ys,\
            freq=freq, ndx_ChangeWD=self.ndx_ChangeWD, ndx_ChangeNSQ=self.ndx_ChangeNSQ, c0=self.c0,\
            cb=self.cb, bloss=self.bloss, rhob=self.rhob, rhow=self.rhow,\
            smoothing_length_ssp=smoothing_length_ssp, smoothing_length_rho=smoothing_length_rho,
            ThinknessOfArtificialAbsorpLayer_ratio_z=ThinknessOfArtificialAbsorpLayer_ratio_z)

        # PE propagator
        propagator = PEPropagator(ref_wavenumber=k0, grid=grid, env_input=env_input)

        output = propagator.run(psi=psi, vertical_slice=vertical_slice, depths=depths)


        # ------- output ------- #

        # take only first 1/2 of z axis (?)
        z = grid.z[:int(grid.Nz/2)]

        # rearrange y axis (azimuthal) so values are increasing order
        y = np.fft.fftshift(grid.q)

        import matplotlib.pyplot as plt
        plot_r = grid.r[1:]
        plot_theta = np.fft.fftshift(np.squeeze(output.Ez_y))

        end = time.time()
        print('Calculation completed in {0:.1f} seconds'.format(end - start))

        SPfield = np.fft.fftshift(np.squeeze(output.Ez[0,:,1:]), axes=0)
        SPfield = 20 * np.log10(np.abs(SPfield))

        if False:
            R, TH = np.meshgrid(plot_r, plot_theta)
            XX = R * np.cos(TH)
            YY = R * np.sin(TH)
            fig = plt.contourf(XX, YY, SPfield, 100)
            plt.colorbar(fig)
            plt.show()

        if False:
            fig=plt.figure()
            ZZ = 20 * np.log10(np.abs(output.Af[:,:]))
            XX, YY = np.meshgrid(grid.r, z)
            fig = plt.contourf(XX, YY, ZZ, 100)
            plt.colorbar(fig)
            ax = plt.gca()
            ax.invert_yaxis()
            plt.show()

        return SPfield


class PEGrid():
    """ Grid for Parabolic Equation solver
        
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
    """
    def __init__(self, radial_step, radial_range,\
            azimuthal_step, azimuthal_range, vertical_step, vertical_range):

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
        y_pos = np.arange(start=0, stop=ny/2, step=1, dtype=float)
        y_neg = np.arange(start=-ny/2, stop=0, step=1, dtype=float)
        y = np.concatenate((y_pos, y_neg)) 
        y *= dy

        # ------ z (vertical) ------ 
        nz = zmax / dz # number of vertical bins
        nz = round(nz / 2) * 2  # ensure even number of vertical bins
        Lz = nz * dz  
        z_pos = np.arange(start=0, stop=nz/2, step=1, dtype=float)
#        z_neg = np.arange(start=-1, stop=-nz/2-1, step=-1, dtype=float)
        z_neg = np.arange(start=-nz/2, stop=0, step=1, dtype=float)
        z = np.concatenate((z_pos, z_neg))
        z *= dz

        # ------ kz (wavenumber) ------
        kz = z * 2 * np.pi / (Lz * dz)

        Y, Z = np.meshgrid(y, z)

        self.r = x
        self.q = y
        self.z = z
        self.kz = kz
        self.Nr = nx
        self.Nq = ny
        self.Nz = nz
        self.Q = Y
        self.Z = Z
        self.dr = dx
        self.dq = dy
        self.dz = dz
