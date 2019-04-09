import numpy as np

class ModelOutput():

    def __init__(self, Y, Z, kz, ny, x, k0):

        self.k0 = k0
        self.Z = Z

        # which points to be output?
        self.yso = np.argwhere(Y[0,:] >= 0)[0]
        self.nzhalf = int(Z.shape[0] / 2)

        # this are used for storing the calculated 3d field values
        self.Ez_z = np.array([.1]) 
#        self.Ez_z = (np.arange(99, dtype=float) + 0.5) * 10. 
        self.Ez_z = self.Ez_z[:, np.newaxis] 
        self.Ez_y = Y[0,:]  # y values (azimuthal)
        self.Ez = np.empty(shape=(len(self.Ez_z), ny, len(x)), dtype=complex)  # sound intensity values

        self.Ez_ifft_kernel = np.exp(1j * np.matmul(self.Ez_z, kz[np.newaxis,:])) / len(kz)

        self.iout_Ez = 0


    def get_output(self, dista, psi, denin, nfft):

        self.iout_Ez += 1

        if dista != 0:

            dz = self.Z[1,0] - self.Z[0,0]

            idx = np.squeeze(np.round(self.Ez_z/dz).astype(int))
            
            if np.ndim(idx) == 0:
                idx = np.array([idx])

            A, B = self._calc_A_B_matrices(psi, denin, idx)

            self.Ez[:,:,self.iout_Ez-1] = self._calc_Ez(A, B, dista)

            psi = self._calc_psi(psi, dista, denin)

            nfft += 1 

        else:
            psi = np.fft.ifft(psi, axis=0)
            nfft += 1 

        # sound intensity on the vertical x-z plane crossing the source
        Af = psi[:self.nzhalf, self.yso]
        Af = np.squeeze(Af)

        return Af, nfft

    def _calc_A_B_matrices(self, psi, denin, idx):
        A = np.matmul(self.Ez_ifft_kernel, psi)
        B = np.sqrt(denin[idx])
        return A, B

    def _calc_Ez(self, A, B, dista):
        y = A * B * np.exp(1j * self.k0 * dista) / np.sqrt(dista) 
        return y

    def _calc_psi(self, psi, dista, denin):
        psi = np.fft.ifft(psi, axis=0) * self._calc_g(dista) * _calc_f(denin)
        return psi

    def _calc_g(self, dista):
        return np.exp(1j * self.k0 * dista) / np.sqrt(dista)

def _calc_f(denin):
        return np.sqrt(denin)

