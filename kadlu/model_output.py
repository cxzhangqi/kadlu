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
        self.Ez_y = Y[0,:]  # y values (azimuthal)
        self.Ez = np.empty(shape=(len(self.Ez_z), ny, len(x)), dtype=complex)  # sound intensity values

        self.Ez_ifft_kernel = np.exp(1j * self.Ez_z * kz) / len(kz)
        self.Ez_ifft_kernel = self.Ez_ifft_kernel[np.newaxis, :]        
        self.iout_Ez = 0


    def get_output(self, dista, psi, denin, nfft):

        self.iout_Ez += 1

        if dista != 0:
            dz = self.Z[1] - self.Z[0]

            idx = np.round(self.Ez_z/dz).astype(int)

            A = np.matmul(self.Ez_ifft_kernel, psi)
            B = np.sqrt(denin[idx,:])

            self.Ez[:,:,self.iout_Ez] = np.matmul(A, B) * np.exp(1j * self.k0 * dista) / np.sqrt(dista) 

            psi = np.fft.ifft(psi) * np.exp(1j * self.k0 * dista) / np.sqrt(dista) * np.sqrt(denin)
            nfft += 1 

        else:
            psi = np.fft.ifft(psi)
            nfft += 1 

        # sound intensity on the vertical x-z plane crossing the source
        Af = psi[:self.nzhalf, self.yso]
        Af = np.squeeze(Af)

        return Af, nfft