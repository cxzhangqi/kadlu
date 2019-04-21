import numpy as np
from numpy.lib import scimath
from kadlu.transmission_loss.seafloor_depth import SeafloorDepth
from kadlu.transmission_loss.refractive_index import RefractiveIndex

class EnvInput():

    def __init__(self, Y, Z, xs, ys, dx, nx, ny, freq, ndx_ChangeWD,\
                    ndx_ChangeNSQ, c0, cb, bloss, rhob, rhow,\
                    smoothing_length_ssp, smoothing_length_rho):

        self.smoothing_length_ssp = smoothing_length_ssp
        self.smoothing_length_rho = smoothing_length_rho

        self.Z = Z

        self.dx = dx

        self.xs = xs
        self.ys = ys

        self._range_independent = (nx == 1) and (ny == 1) # range independent environment

        # allocate memory
        m = Z.shape
        n = Y.shape[1]

        self.seafloor = SeafloorDepth(n)
        self.refractive_index = RefractiveIndex(m)

        self._unchanged = np.zeros(shape=n, dtype=bool)

        self.wd_x_next = 0

        self.wd_old = np.empty(shape=(1,n))
        self.wd_mask = np.empty(m)
        self.DwdDy = np.empty(m)
        self.Z_sub_wd = np.empty(m)

        self.NSQ_x_next = 0

        self.n2w = np.zeros(m)
        self.n2w_new = np.zeros(m)
        
        # updated in _update_phase_screen()
        self.H_c = np.empty(m, dtype=complex)
        self.n2in = np.empty(m, dtype=complex)
        self.H_rho = np.empty(m, dtype=float)
        self.ddenin = np.empty(m, dtype=complex)   
        self.d2denin = np.empty(m, dtype=complex)
        self.denin = np.empty(m, dtype=float)


        self.costheta = np.cos(Y[0,:])
        self.sintheta = np.sin(Y[0,:])

        self.ndx_ChangeWD = ndx_ChangeWD
        self.ndx_ChangeNSQ = ndx_ChangeNSQ
        
        # complex bottom sound speed
        self.rhob = rhob
        kbi = bloss / (cb / freq) / 20. / np.log10(np.e)   # note: np.e is faster than np.exp()
        betab = kbi / 2 / np.pi / freq
        cbi = np.roots([betab, -1, betab*cb**2])  # roots of polynomial p[0]*x^n+...p[n]
        cbi = cbi[np.imag(cbi) == 0] 
        cbi = cbi[np.logical_and(cbi >= 0, cbi < cb)]
        cb = cb - 1j * cbi
        self.n2b = (c0 / cb)**2

        self.rhow = rhow


    def get_input(self, dista):

        new_env = self._new_env(dista)[0]
        new_env_any = np.any(new_env)

        if new_env_any:
            self._update_phase_screen(dista, new_env)

        return new_env_any, new_env


    def _new_bathymetry(self, x, y, dista):

        new_bathy = self._unchanged

        if (dista == self.dx/2) or (dista >= self.wd_x_next):
            
            # next distance at which to update bathymetry
            self.wd_x_next = dista + self.ndx_ChangeWD * self.dx
            
            # get bathymetry
            wd_new, DwdDy_new = self.seafloor.get_depth(x=x, y=y)
            wd_new = wd_new[np.newaxis,:]
            DwdDy_new = DwdDy_new[np.newaxis,:]
            #MATLAB: [wd_new(:),DwdDy_new(:)]  = sub_SeafloorDepth(x,y);

            new_bathy = np.logical_and(self.wd_old != wd_new, np.logical_not(np.isnan(wd_new)))

###            print('wd_mask.shape: ', self.wd_mask.shape)
###            print('wd_new.shape: ', wd_new.shape)
###            print('new_bathy.shape: ', new_bathy.shape)

            self.wd_old = wd_new

            k = self.Z.shape

            self.wd_mask[:,new_bathy[0]] = np.ones((k[0],1)) * wd_new[new_bathy]  # water depth mask
            self.DwdDy[:,new_bathy[0]] = np.ones((k[0],1)) * DwdDy_new[new_bathy]
            self.Z_sub_wd[:,new_bathy[0]] = np.abs(self.Z[:,new_bathy[0]]) - self.wd_mask[:,new_bathy[0]]
            
            print('Updating bathymetry at {0:.2f} m'.format(dista))

        return new_bathy


    def _new_refractive_index(self, x, dista): # NSQ

        new_refr = self._unchanged

        if (dista == self.dx/2) or ((dista >= self.NSQ_x_next) and not self._range_independent):  # update water column

            # next distance at which to update water column
            self.NSQ_x_next = dista + self.ndx_ChangeNSQ * self.dx

            # interpolate water nsq
            self.n2w_new = np.copy(self.n2w)
            IDZ = (self.Z <= 0)
            idy = np.nonzero(IDZ)[1]
            IDZ = np.nonzero(IDZ)
            self.n2w_new[IDZ] = self.refractive_index.get_nsq(x=x[idy], YZ=self.Z[IDZ])   #sub_NSQ(x(idy).',y(idy).',Z(IDZ));

            k = self.Z.shape
            one = np.array([0], dtype=int)
            indeces = np.arange(start=k[0]-1, step=-1, stop=k[0]/2, dtype=int)
            indeces = np.concatenate([one, indeces])

            new_refr = np.any(self.n2w[indeces,:] - self.n2w_new[indeces,:] != 0, axis=0)

            indeces2 = np.arange(start=1, step=1, stop=k[0]/2, dtype=int)

            self.n2w_new[np.ix_(indeces2, new_refr)] = self.n2w_new[np.ix_(indeces[1:], new_refr)]

            self.n2w[:,new_refr] = self.n2w_new[:,new_refr]       

            print('Updating water column at {0:.2f} m'.format(dista))

        return new_refr


    def _new_env(self, dista):
        x = self.xs + self.costheta * dista
        y = self.ys + self.sintheta * dista
        new = np.logical_or(self._new_bathymetry(x, y, dista), self._new_refractive_index(x, dista))
        return new


    def _update_phase_screen(self, dista, new_env):

        # smooth ssp
        self.H_c[:,new_env] = (1 + np.tanh(self.Z_sub_wd[:,new_env] / self.smoothing_length_ssp / 2)) / 2
        self.n2in[:,new_env] = self.n2w[:,new_env] + (self.n2b - self.n2w[:,new_env]) * self.H_c[:,new_env]
        itmp = (self.wd_mask[0,:] == 0) 
        if np.any(itmp):
            self.n2in[0,itmp] = self.n2b
        
        # smooth density
        TANH = np.tanh(self.Z_sub_wd[:,new_env] / self.smoothing_length_rho / 2)
        self.H_rho[:,new_env] = (1 + TANH) / 2
        self.denin[:,new_env] = self.rhow + (self.rhob - self.rhow) * self.H_rho[:,new_env]
        
        SECH2 = 1 / np.cosh(self.Z_sub_wd[:,new_env] / self.smoothing_length_rho / 2)
        SECH2 = SECH2 * SECH2
        self.ddenin[:,new_env] =  SECH2 / self.smoothing_length_rho / 2 * np.sqrt(1 + self.DwdDy[:,new_env]**2)
        self.ddenin[:,new_env] =  (self.rhob - self.rhow) / 2 * self.ddenin[:,new_env]

        self.d2denin[:,new_env] = -SECH2 / self.smoothing_length_rho / 2 * (TANH / self.smoothing_length_rho * (1+self.DwdDy[:,new_env]**2))
        self.d2denin[:,new_env] =  (self.rhob - self.rhow) / 2 * self.d2denin[:,new_env]
        
        print('Updating phase screen at {0:.2f} m'.format(dista))


