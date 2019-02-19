# fit mcvine-simulated point spread function to an empirical form

import os, numpy as np
import histogram.hdf as hh, histogram as H
import scipy.optimize as sopt


class PSF_Affine_Model(object):

    def __init__(self, Ei, E):
        self.Ei = Ei
        self.E = E
        self.E_profile = None
        self.q_profile = None
        self.dq_over_dE = 0
        return

    def set_qE_profile(self, dq_over_dE):
        self.qE_profile = self.make_qE_profile_function(dq_over_dE)
        return self.qE_profile

    def make_qE_profile_cost(self, res_x, res_y, res_z):
        xgrid,ygrid = np.meshgrid(res_x, res_y)
        dx = res_x[1]-res_x[0]
        dy = res_y[1]-res_y[0]
        def cost(x):
            dq_over_dE, = x
            f = self.make_qE_profile_function(dq_over_dE)
            zfit = f(xgrid, ygrid)
            # print "sum of area", zfit.sum()*dx*dy
            zfit /= zfit.sum()*dx*dy
            c =  np.sum((zfit-res_z)**2)/res_z.size
            # print "cost", c
            return c
        return cost

    def make_qE_profile_function(self, dq_over_dE):
        def _(qgrid, Egrid):
            return ellipsoid_affine(qgrid, Egrid, dq_over_dE, self.q_profile, self.E_profile)
        return _

    def set_E_profile(self, *args, **kwds):
        raise NotImplementedError("set_E_profile")

    def set_q_profile(self, sigma):
        self.q_profile = Qprofile(sigma)
        return self.q_profile
    
    def make_E_profile_function(self, a, b, R, sigma=None, t0=None):
        raise NotImplementedError("make_E_profile_function")

    def make_E_profile_cost(self, res_x, res_y):
        def cost(x):
            f = self.make_E_profile_function(*x)
            yfit = f(res_x)
            dx = res_x[1]-res_x[0]
            yfit /= np.sum(yfit) * dx
            return np.sum((yfit-res_y)**2)/res_y.size
        return cost


class Qprofile:

    def __init__(self, sigma):
        self.sigma = sigma
        self._h = 1./np.sqrt(2*np.pi)/np.abs(sigma)

    def __call__(self, x):
        return gaus(x, self._h, 0., self.sigma)


def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def ellipsoid_affine(q, E, dq_over_dE, q_profile, E_profile):
    """q and E are grids
    I(q) are gaussians. their centers are linearly related to E
    """
    pE = E_profile(E)
    q1 = E*dq_over_dE
    pq = q_profile(q-q1)
    return pq*pE

