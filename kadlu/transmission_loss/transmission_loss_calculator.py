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
    """ Grid for Parabolic Equation solver.

        Creates a regular cylindrical grid, where

            r: radial distance
            q: azimuthal angle
            z: vertical depth

        The radial (r) and vertical coordinates are in meters
        while the angular coordinate (q) is in radians.

        Args:
            radial_step: float
                Radial step size
            radial_range: float
                Radial range
            azimuthal_step: float
                Angular bin size
            azimuthal_range: float
                Angular domain
            vertical_step: float
                Vertical bin size
            vertical_range: float
                Vertical range

        Attributes:
            r: 1d numpy array
                Radial grid values
            dr: float
                Radial step size
            Nr: int
                Number of radial steps
            q: 1d numpy array
                Angular grid values
            dq: float
                Angular bin size
            Nq: int
                Number of angular bins
            z: 1d numpy array
                Vertical grid values
            dz: float
                Vertical bin size
            Nz: int
                Number of vertical bins
            z: 1d numpy array
                Wavenumber grid values
            Q: 2d numpy array
                Angular values of azimuthal-vertical matrix; has shape (Nz,Nq).
            Z: 2d numpy array
                Vertical values of azimuthal-vertical matrix; has shape (Nz,Nq).
    """
    def __init__(self, radial_step, radial_range,\
            azimuthal_step, azimuthal_range, vertical_step, vertical_range):

        # radial
        self.r, self.dr, self.Nr = self.__make_radial_grid__(radial_step, radial_range)

        # azimuthal
        self.q, self.dq, self.Nq = self.__make_azimuthal_grid__(azimuthal_step, azimuthal_range)

        # vertical
        self.z, self.dz, self.Nz = self.__make_vertical_grid__(vertical_step, vertical_range)

        # wavenumber
        L = self.Nz * self.dz  
        self.kz = self.z * 2 * np.pi / (L * self.dz)

        # mesh-grid
        self.Q, self.Z = np.meshgrid(self.q, self.z)

    def __make_radial_grid__(self, dr, rmax):
        N = round(rmax / dr)
        r = np.arange(N+1, dtype=float)
        r *= dr
        return r, dr, N

    def __make_azimuthal_grid__(self, dq, qmax):
        N = int(np.ceil(qmax / dq))
        if N%2 == 1: 
            N = N + 1 # ensure even number of angular bins

        q_pos = np.arange(start=0, stop=N/2, step=1, dtype=float)
        q_neg = np.arange(start=-N/2, stop=0, step=1, dtype=float)
        q = np.concatenate((q_pos, q_neg)) 
        q *= dq
        return q, dq, N

    def __make_vertical_grid__(self, dz, zmax):
        N = zmax / dz
        N = round(N / 2) * 2  # ensure even number of vertical bins
        z_pos = np.arange(start=0, stop=N/2, step=1, dtype=float)
        z_neg = np.arange(start=-N/2, stop=0, step=1, dtype=float)
        z = np.concatenate((z_pos, z_neg))
        z *= dz
        return z, dz, N