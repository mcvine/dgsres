# fit mcvine-simulated point spread function to an empirical form

import os, numpy as np
import histogram.hdf as hh, histogram as H
import scipy.optimize as sopt


from . import fit_2d_psf
class Fit2EE(fit_2d_psf.Fit):

    """Eenrgy profile is sum of two exponential X error_func
    """

    def createModel(self):
        self.model = EEAffineModel(self.Ei, self.E)
        return

    def fit_E_profile(self, sigma_left, sigma_right, weight_left, ef_width, E_offset):
        # fit res(E)
        res_E = self.mcvine_psf_E
        res_y = res_E.I.copy()
        res_x = res_E.E.copy()
        dx = res_x[1] - res_x[0]; res_y /= res_y.sum()*dx # normalize. density function
        if isinstance(sigma_left, tuple):
            # inputs are bounds
            bounds = [sigma_left, sigma_right, weight_left, ef_width, E_offset]
            x0 = [np.average(b) for b in bounds]
        else:
            # inputs are single numbers
            x0 = [sigma_left, sigma_right, weight_left, ef_width, E_offset]
            tmpE = np.max(np.abs([sigma_left, sigma_right, ef_width, E_offset]))
            bounds = np.array([
                [sigma_left/3, sigma_right/3, 0, ef_width/10., -tmpE*2],
                [sigma_left*3, sigma_right*3, 1., ef_width*10., tmpE*2]
            ]).T
        cost = self.model.make_E_profile_cost(res_x, res_y)
        fitres = sopt.minimize(cost, x0=x0, bounds=bounds, method='L-BFGS-B')
        E_profile = self.model.set_E_profile(*fitres.x)
        yfit = E_profile(res_x); yfit /= yfit.sum()*dx
        # fit parameters, x values, y values, fit y values, E_profile fitted
        return fitres.x, res_x, res_y, yfit, E_profile


class EEAffineModel(fit_2d_psf.PSF_Affine_Model):

    """Energy profile is sum of two exponential X error_func
    q profile is gaussian
    Center of I(q) gaussians at different q is linearly associated with E
    """

    def update_parameters(self, sigma_left=None, sigma_right=None, weight_left=None, ef_width=None, E_offset=None, q_sigma=None, dq_over_dE=None):
        self.set_E_profile(sigma_left, sigma_right, weight_left, ef_width, E_offset)
        self.set_q_profile(q_sigma)
        self.set_qE_profile(dq_over_dE)
        return

    def set_E_profile(self, sigma_left, sigma_right, weight_left, ef_width, E_offset):
        self.E_profile = self.make_E_profile_function(sigma_left, sigma_right, weight_left, ef_width, E_offset)
        return self.E_profile

    def make_E_profile_function(self, sigma_left, sigma_right, weight_left, ef_width, E_offset):
        return EEModel(sigma_left, sigma_right, weight_left, ef_width, E_offset)


class EEModel:

    def __init__(self, sigma_left, sigma_right, weight_left, ef_width, E_offset):
        self.sigma_left = sigma_left
        self.sigma_right = sigma_right
        self.weight_left = weight_left
        self.ef_width = ef_width
        self.E_offset = E_offset
        return
    
    def __call__(self, x):
        from scipy.special import erfc
        weight_left = self.weight_left
        weight_right = 1.0 - weight_left
        sigma_left = self.sigma_left
        sigma_right = self.sigma_right
        ef_width = self.ef_width
        E_offset = self.E_offset
        x = x-E_offset
        t1 = fit_2d_psf.gaus(x, 1./np.sqrt(2*np.pi)/np.abs(sigma_left), 0., sigma_left) * (1-erfc(x/ef_width)/2.)
        t2 = fit_2d_psf.gaus(x, 1./np.sqrt(2*np.pi)/np.abs(sigma_right), 0., sigma_right) * erfc(x/ef_width)/2.
        return t1*weight_left + t2*weight_right

