import numpy as np
from numpy.lib import scimath
from kadlu.transmission_loss.model_output import OutputCollector


class PEPropagator():

    def __init__(self, ref_wavenumber, grid, env_input):

        self.k0 = ref_wavenumber
        self.grid = grid
        self.env_input = env_input


    def run(self, psi, vertical_slice=True, depths=[.1]):

        # output collector
        output = OutputCollector(ref_wavenumber=self.k0, grid=self.grid, env_input=self.env_input,\
            vertical_slice=vertical_slice, depths=depths)

        print('psi.shape: ', psi.shape)

        # initial Nx2D free propagator
        fr_half, fr_full = self.__free_propagator__()

        print('fr_half.shape:', fr_half.shape)
        print('fr_full.shape:', fr_full.shape)

        # output field at 0
        output.collect(dist=0, psi=psi)

        # PE marching starts here
        dist = 0
        Nr = self.grid.Nr
        dr = self.grid.dr
        for step in range(1, Nr+1):
    
            print("loop: {0}/{1}".format(step, Nr+1), end="\r")

            # (1) r --> r + dr/2 free propagation
            psi = fr_half * psi

            # (2) phase adjustment at r + dr/2
            U = self.env_input.update(dist + dr/2)            
            psi = np.fft.fft(U * np.fft.ifft(psi, axis=0), axis=0)    

            # (3) r + dr/2 --> r + dr free propagation
            psi = fr_half * psi

            # increment distance
            dist = dist + dr

            # collect output
            output.collect(dist=dist, psi=psi)

#        psi_final = np.fft.ifft(psi, axis=0) * np.exp(1j * self.k0 * dist) / np.sqrt(dist) * np.sqrt(env.denin)
#        psi_final = np.fft.fftshift(psi_final[:int(self.grid.Nz/2),:], axes=(1,))

        return output


    def __free_propagator__(self):

        k0 = self.k0
        dr = self.grid.dr
        kz = self.grid.kz
        ny = self.grid.Nq

        fr_half = np.exp(1j * dr / 2 * (scimath.sqrt(k0**2 - kz**2) - k0))
        fr_half = fr_half[:,np.newaxis]
        fr_half = fr_half * np.ones(shape=(1,ny))

        fr_full = np.exp(1j * dr * (scimath.sqrt(k0**2 - kz**2) - k0)) 
        fr_full = fr_full[:,np.newaxis]        
        fr_full = fr_full * np.ones(shape=(1,ny))

        return fr_half, fr_full


