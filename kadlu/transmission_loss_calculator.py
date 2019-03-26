import numpy as np
from numpy.lib import scimath
from kadlu.pe_starter import PEStarter
from kadlu.env_input import EnvInput
import math


class TransmissionLossCalculator():

    def __init__(self, bathymetry, sound_speed):

        self.bathy = bathymetry
        self.sound_speed = sound_speed

        self.c0 = 1500      # reference sound speed in m/s
        self.rhow = 1       # water density
        self.cb = 1700      # homogeneous bottom sound speed
        self.bloss = 0.5    # homogeneous bottom attenuation db/lambda
        self.rhob = 1.5     # homogeneous bottom density
        self.ndx_ChangeWD = 3 #1  how often to update bathymetry

        self.ndx_ChangeNSQ = math.inf  # how often to update water column (inf => range-independent SSP)

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
        Y, Z, x, y, z, kz, nx, ny, nz = self._create_grids(radial_step=dr, radial_range=50e3,\
                azimuthal_step=1/180*np.pi, azimuthal_range=2*np.pi,\
                vertical_step=10, vertical_range=vertical_range)

        print('\nnx,ny,nz: ', nx, ny, nz)
        print('x.shape: ', x.shape)
        print('y.shape: ', y.shape)
        print('z.shape: ', z.shape)
        print('y: {0:.3f},{1:.3f}...{2:.3f},{3:.3f},{4:.3f}...,{5:.3f},{6:.3f}'.format(y[0],y[1],y[int(ny/2)-1],y[int(ny/2)],y[int(ny/2)+1],y[-2],y[-1]))
        print('dx: {0:.1f} m'.format(dr))

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

        # allocate memory
        Af = np.empty(shape=(int(nz/2), len(x)))
        U = np.zeros(shape=Z.shape)
        print('Af.shape: ', Af.shape)
        print('U.shape: ', U.shape)

        # module handling updates of environmental input
        env = EnvInput(Y=Y, Z=Z, xs=xs, ys=ys, dx=dr, nx=nx, ny=ny,\
            freq=freq, ndx_ChangeWD=self.ndx_ChangeWD, ndx_ChangeNSQ=self.ndx_ChangeNSQ, c0=self.c0,\
            cb=self.cb, bloss=self.bloss, rhob=self.rhob, rhow=self.rhow)

        # output field at 0
        nfft = 0
        dista = 0
        iNextOutput = 0
        Af[:,iNextOutput] = self._model_output() #icase,psi,dista,ndy_3DSliceout,ndz_3DSliceout,YZSlice_output_folder,isplot);
        iNextOutput += 1

        # PE marching starts here
        is_halfstep = True    # initially, half step to dx/2
        for jj in range(10): #range(1, nx+1):
    
            print("loop: {0}/{1}".format(jj, nx+1), end="\r")

            # (1) x --> x+dx/2 free propagation
            if is_halfstep:
                psi = fr_half * psi
    
            # (2) do phase adjustment at x+dx/2
            isnewscreen, isupdate = env.get_input(dista=dista) #(freq,c0,dista+dx/2,dx,xs,ys,smoothing_length_rho,smoothing_length_ssp);
            
#            if isnewscreen:
#                U[:, isupdate] = np.exp(1j * dr * k0 * (-1 + scimath.sqrt( n2in[:,isupdate] + atten0[:,isupdate] +\
#                    1/2 / k0**2 * (d2denin[:,isupdate] / denin[:,isupdate] - 3/2 * (ddenin[:,isupdate] / denin[:,isupdate])**2))))
            
            psi = np.fft.fft(U * np.fft.ifft(psi))    
            nfft += 2 
            
            # (3) x+dx/2 --> x+dx free propagation
            dista = dista + dr
        
            # output field?
            if np.any(x == dista):
                is_halfstep = True        # output field at x+dx, so half step from x+dx/2
                psi = fr_half * psi
                Af[:, iNextOutput] = self._model_output() #icase,psi,dista,ndy_3DSliceout,ndz_3DSliceout,YZSlice_output_folder,isplot);
                iNextOutput += 1

            else:                        # if not output filed, full step to x+dx/2+dx
                is_halfstep = False
                psi = fr_full * psi



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
        y_pos = np.arange(start=0, stop=ny/2, step=1, dtype=float)
        y_neg = np.arange(start=-1, stop=-ny/2-1, step=-1, dtype=float)
        y = np.concatenate((y_pos, y_neg)) 
        y *= dy

        # ------ z (vertical) ------ 
        nz = zmax / dz # number of vertical bins
        nz = round(nz / 2) * 2  # ensure even number of vertical bins
        Lz = nz * dz  
        z_pos = np.arange(start=0, stop=nz/2, step=1, dtype=float)
        z_neg = np.arange(start=-1, stop=-nz/2-1, step=-1, dtype=float)
        z = np.concatenate((z_pos, z_neg))
        z *= dz

        # ------ kz (wavenumber) ------
        kz = z * 2 * np.pi / (Lz * dz)


        Y, Z = np.meshgrid(y, z)

        return Y, Z, x, y, z, kz, nx, ny, nz

    def _free_propagator(self, dr, k0, kz, ny):
        fr_half = np.exp(1j * dr / 2 * (scimath.sqrt(k0**2 - kz**2) - k0))
        fr_half = fr_half[:,np.newaxis]
        fr_half = fr_half * np.ones(shape=(1,ny))
        fr_full = np.exp(1j * dr * (scimath.sqrt(k0**2 - kz**2) - k0)) 
        fr_full = fr_full[:,np.newaxis]        
        fr_full = fr_full * np.ones(shape=(1,ny))
        return fr_half, fr_full


    def _model_output(self):
        return 0
