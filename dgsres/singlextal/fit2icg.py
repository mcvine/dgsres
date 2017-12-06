# fit mcvine-simulated point spread function to an empirical form

import os, numpy as np
import histogram.hdf as hh, histogram as H
from .. import icg
import scipy.optimize as sopt


# cncs_geom = icg.Geom(l1=6.413, l2=36.2-6.413, l3=3.5)

from . import fit_2d_psf
class Fit2ICG(fit_2d_psf.Fit):

    def __init__(self, sim_path, instrument_geom, qaxis=(-0.3, 0.3, 0.01), Eaxis=(-2, 1, 0.05), Ei=None, E=None):
        self.instrument_geom = instrument_geom
        super(Fit2ICG, self).__init__(sim_path, qaxis, Eaxis, Ei, E)
        return

    def createModel(self):
        self.model = ICGAffineModel(self.instrument_geom, self.Ei, self.E)
        return

    def fit_E_profile(self, a=None, b=None, R=None, sigma=None, t0=None):
        # fit res(E)
        res_E = self.mcvine_psf_E
        res_y = res_E.I.copy()
        res_x = res_E.E.copy()
        dx = res_x[1] - res_x[0]; res_y /= res_y.sum()*dx # normalize. density function
        if isinstance(a, tuple):
            # inputs are bounds
            bounds = [a, b, R, sigma, t0]
            x0 = [np.average(b) for b in bounds]
        else:
            # inputs are single numbers
            x0 = a,b,R,sigma,t0
            bounds = np.array([
                [a/3, b/3, 0.01, sigma/10., 0.],
                [a*3, b*3, 0.99, sigma*10., t0*10],
            ]).T
        cost = self.model.make_E_profile_cost(res_x, res_y)
        fitres = sopt.minimize(cost, x0=x0, bounds=bounds, method='L-BFGS-B')
        a,b,R,sigma,t0 = fitres.x
        E_profile = self.model.set_E_profile(a,b,R,sigma,t0)
        yfit = E_profile(res_x); yfit /= yfit.sum()*dx
        # fit parameters, x values, y values, fit y values, E_profile fitted
        return fitres.x, res_x, res_y, yfit, E_profile


class ICGAffineModel(fit_2d_psf.PSF_Affine_Model):

    """Energy profile is ICG
    q profile is gaussian
    Center of I(q) gaussians at different q is linearly associated with E
    """

    def __init__(self, instrument_geom, Ei, E):
        self.instrument_geom = instrument_geom
        super(ICGAffineModel, self).__init__(Ei, E)
        return

    def set_E_profile(self, a, b, R, sigma=None, t0=None):
        self.E_profile = self.make_E_profile_function(a,b,R,sigma,t0)
        return self.E_profile

    def make_E_profile_function(self, a, b, R, sigma=None, t0=None):
        E = self.E
        Ei = self.Ei
        return ICGModel(Ei, E, a, b, R, sigma, t0, self.instrument_geom)


class ICGModel:

    def __init__(self, Ei, E, a, b, R, sigma=None, t0=None, geom=None):
        self.Ei = Ei
        self.E = E
        self.a = a
        self.b = b
        self.R = R
        self.sigma = sigma
        self.t0 = t0
        self.geom = geom
        return
    
    def __call__(self, x):
        a,b,R,sigma,t0 = self.a, self.b, self.R, self.sigma, self.t0
        return icg.resolution(self.E+x, self.Ei, self.E, a, b, R, sigma=sigma, t0=t0, geom=self.geom)

