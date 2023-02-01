# -*- Python -*-
#
# Jiao Lin <jiao.lin@gmail.com>
#
# See https://user-images.githubusercontent.com/1796155/52989676-8809c080-33d2-11e9-8f47-1c0687be9ba0.jpeg
#
# This is more sophisticated than fit_2d_psf


import os, numpy as np, lmfit
from .PSF_Affine_Model import gaus
from .fit2ee import EEModel

import copy
from .fit_2d_psf import Fit as base
from multiprocessing import Process, Pipe

class Fit(base):

    def fit(self, rounds=None, gaussian2d_threshold=0.5, alpha_bounds=(-np.pi/2, np.pi/2), **kwds):
        """
        the fitting is done by looking at two axes, one along the main axis of the covariance matrix,
        without considering the tail, the other along the tail direction.

        alpha_bounds: 

            alpha (and beta) values can be 180 degree off.
            to make interpolation easy it is better to get all alpha values in a certain range
            if alpha_bounds is given, the code here will try to keep alpha within bounds.
            beta should have bounds too; but it is right now somewhat hard-coded.

        gaussian2d_threshold:

            The threshold here determins the data included in calculating the main axis of the 
            covariance matrix.
        """
        qgrid, Egrid = self.qEgrids
        self.res_z = res_z = self.get_res_z()
        return fit(qgrid, Egrid, res_z, rounds=rounds, gaussian2d_threshold=gaussian2d_threshold, alpha_bounds=alpha_bounds, **kwds)

    def get_res_z(self):
        reshist = self.mcvine_psf_qE
        res_x = reshist.q
        res_y = reshist.E
        res_z = reshist.I.T.copy()
        dx = res_x[1]-res_x[0]
        dy = res_y[1]-res_y[0]
        res_z /= np.sum(res_z)*dx*dy
        return res_z

    def createModel(self):
        # model is in the ellipsoid method below
        return


class FitResult: pass


class InterpModel:

    """interpolated resolution model

    Given parameter values at some (q,E) points, this class provides convenient methods
    to calculate the parameters at a given (q',E') point, and generate a model.
    """
    def __init__(self, qE_points, param_values, qrange, Erange):
        self.qE_points = qE_points
        self.param_values = param_values
        self.qrange = qrange
        self.Erange = Erange
        return

    def getModel(self, q, E):
        from scipy.interpolate import griddata
        kwds = dict()
        keys = list(self.param_values.keys())
        for k in keys:
            vals = self.param_values[k]
            kwds[k] = griddata(self.qE_points, vals, [[q, E]]) [0]
            continue
        return Model(self.qrange, self.Erange, **kwds)


class Model:

    def __init__(
            self, qrange, Erange,
            alpha, beta,
            xp_center, xp_sigma,
            y_sigma_left, y_sigma_right, y_weight_left, y_ef_width, y_offset,
            scale):
        self.qrange = qrange
        self.Erange = Erange
        self.alpha, self.beta = alpha, beta
        self.xp_center, self.xp_sigma = xp_center, xp_sigma
        self.y_sigma_left, self.y_sigma_right, self.y_weight_left, self.y_ef_width, self.y_offset\
            = y_sigma_left, y_sigma_right, y_weight_left, y_ef_width, y_offset
        self.scale = scale
        return

    def ongrid(self, qgrid, Egrid):
        u = qgrid/self.qrange
        v = Egrid/self.Erange
        ret = ellipsoid(
            u, v, self.alpha, self.beta,
            self.xp_center, self.xp_sigma,
            self.y_sigma_left, self.y_sigma_right, self.y_weight_left, self.y_ef_width, self.y_offset,
            self.scale)
        ret/=ret.sum()
        return ret


def qE2uv_grid(qgrid, Egrid):
    qrange, Erange = qEgrid2range(qgrid, Egrid)
    return qgrid/qrange, Egrid/Erange

def qEgrid2range(qgrid, Egrid):
    qrange = qgrid[0][-1] - qgrid[0][0]
    Erange = Egrid[:, 0][-1] - Egrid[:, 0][0]
    return qrange, Erange

class _fitcore(object):
    # helper class for `fit` method below
    def __init__(self, ugrid,vgrid,z,model,fit_rounds=None):
        self.ugrid = ugrid
        self.vgrid = vgrid
        self.z = z
        self.model = model
        self.fit_rounds = fit_rounds
        self.results = None

    def fit(self,conn=None):
        if not self.fit_rounds: self.fit_rounds = 3
        self.results = []
        for i in range(self.fit_rounds):
            print(" -- Fitting round %s" % i)
            result = self.model.fit(self.z, u=self.ugrid.flatten(), v=self.vgrid.flatten(), method='differential_evolution')
            # print result.fit_report()
            self.results.append(result)
            print("    chisq=%s" % result.chisqr)
        #print(self.results)
        if conn: conn.send(self)

def fit(qgrid, Egrid, I, rounds=None, gaussian2d_threshold=0.5, alpha_bounds=None, return_all_results=False,multiprocess=True):
    if not rounds: rounds = 3
    # convert to unitless
    ugrid, vgrid = qE2uv_grid(qgrid, Egrid)
    z = I.copy().flatten()/I.max()
    model, guess_results = guessModel(qgrid, Egrid, I, gaussian2d_threshold, alpha_bounds=alpha_bounds)
    print("Established model:", model)
    model.print_param_hints(colwidth=12)
    print("Start fitting...")
    fit_0 =_fitcore(ugrid,vgrid,z,model,fit_rounds=rounds)
    if multiprocess:
        conn0_s,conn0_r = Pipe()
        p0 = Process(target=fit_0.fit,args=(conn0_s,)) 
    else:
        fit_0.fit()
    # alpha may be 90 degrees off

    print("Start fitting with alpha_guess+90degree...")
    alpha, beta = guess_results[:2]
    model, guess_results = guessModel(qgrid, Egrid, I, gaussian2d_threshold, alpha=alpha+np.pi/2, beta=beta, alpha_bounds=alpha_bounds)
    fit_90 = _fitcore(ugrid,vgrid,z,model,fit_rounds=rounds)
    if multiprocess:
        conn90_s,conn90_r = Pipe()
        p90 =Process(target=fit_90.fit,args=(conn90_s,))
    else:
        fit_90.fit()
    if multiprocess:
        p0.start()
        p90.start()
        fit_0 = conn0_r.recv()
        fit_90 = conn90_r.recv()
        p0.join()
        p90.join()
    print(fit_0.results)
    print(fit_90.results)
    results = fit_0.results + fit_90.results
    results.sort(key=lambda x: x.chisqr)
    if return_all_results:
        return results
    result = results[0]
    # fix alpha and refit
    alpha1 = result.best_values['alpha']
    alpha = adjust_to_within_bounds(alpha1, alpha_bounds, period=np.pi)
    if alpha != alpha1:
        model.set_param_hint('alpha', value=alpha, vary=False)
        result = model.fit(z, u=ugrid.flatten(), v=vgrid.flatten(), method='differential_evolution')
    return result

def guessModel(qgrid, Egrid, I, gaussian2d_threshold, alpha=None, beta=None, alpha_bounds=None):
    alpha, beta, xp_bc, Ixp, y_bc, Iy, xp_center, xp_sigma, y_center, y_sigma = guess_results = fitguess(qgrid, Egrid, I, gaussian2d_threshold, alpha=alpha, beta=beta, alpha_bounds=alpha_bounds)
    model = lmfit.Model(ellipsoid, independent_vars=['u', 'v'])
    model.set_param_hint('alpha', value=alpha, min=alpha-np.pi/10, max=alpha+np.pi/10)
    model.set_param_hint('beta', value=beta, min=beta-np.pi/10, max=beta+np.pi/10)
    model.set_param_hint('xp_center', value=xp_center, min=xp_center-1, max=xp_center+1)
    model.set_param_hint('xp_sigma', value=xp_sigma, min=xp_sigma/10, max=xp_sigma*10)
    model.set_param_hint('y_sigma_left', value=y_sigma/2., min=y_sigma/10, max=y_sigma*10)
    model.set_param_hint('y_sigma_right', value=y_sigma/2., min=y_sigma/10, max=y_sigma*10)
    model.set_param_hint('y_weight_left', value=.5, min=0.01, max=.99)
    model.set_param_hint('y_ef_width', value=y_sigma/10, min=y_sigma/100, max=y_sigma*2)
    model.set_param_hint('y_offset', value=0, min=-y_sigma, max=y_sigma)
    model.set_param_hint('scale', value=1., min=0.1, max=10)
    return model, guess_results
    
def ellipsoid(
        u, v, alpha, beta,
        xp_center, xp_sigma,
        y_sigma_left, y_sigma_right, y_weight_left, y_ef_width, y_offset,
        scale):
    "ellipsoid function. this one has tail"
    cos_alpha= np.cos(alpha); sin_alpha = np.sin(alpha)
    x = u * cos_alpha + v * sin_alpha
    y = -u * sin_alpha + v * cos_alpha
    gamma = beta - alpha
    x_prime = x - y/np.tan(gamma)
    eem = EEModel(y_sigma_left, y_sigma_right, y_weight_left, y_ef_width, y_offset)
    return gaus(x_prime, 1., xp_center, xp_sigma) * eem(y) * scale

def ellipsoid_notail(u, v, alpha, beta, xp_center, xp_sigma, y_center, y_sigma, scale):
    "ellipsoid function. this is a simple 2D gaussian without tail"
    cos_alpha= np.cos(alpha); sin_alpha = np.sin(alpha)
    x = u * cos_alpha + v * sin_alpha
    y = -u * sin_alpha + v * cos_alpha
    gamma = beta - alpha
    x_prime = x - y/np.tan(gamma)
    return gaus(x_prime, 1., xp_center, xp_sigma) * gaus(y, 1., y_center, y_sigma) * scale
                            
# use cov matrix
def getPrimaryAxisAngle(ugrid, vgrid, I):
    cov = np.cov(
        np.array([ugrid.flatten(), vgrid.flatten()]),
        aweights=I.flatten())
    vals, vecs = np.linalg.eig(cov)
    return np.arctan2(*vecs[0])

def getAlpha(ugrid, vgrid, I, threshold=.5):
    "guess alpha angle. the bright part"
    # median = np.nanmedian(I[I>0])
    max = np.nanmax(I)
    # large=np.where(I>median)
    large = np.where(I>max*threshold)
    # print "Alpha. max=", max*.5
    return getPrimaryAxisAngle(ugrid[large], vgrid[large], I[large])

def getBeta(ugrid, vgrid, I):
    "guess beta angle. everything same weight"
    positive=np.where(I>0)
    return getPrimaryAxisAngle(ugrid[positive], vgrid[positive], (I*0+1.)[positive])

def transform_grid(ugrid, vgrid, alpha, beta):
    "convert u,v to x' and y"
    cos_alpha= np.cos(alpha); sin_alpha = np.sin(alpha)
    xgrid = ugrid * cos_alpha + vgrid * sin_alpha
    ygrid = -ugrid * sin_alpha + vgrid * cos_alpha
    gamma = beta-alpha
    x_prime_grid = xgrid - ygrid/np.tan(gamma)
    return x_prime_grid, ygrid

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.
    
    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

def fitguess(qgrid, Egrid, I, gaussian2d_threshold=0.5, alpha=None, beta=None, alpha_bounds=None, beta_bounds=None):
    """return guess parameters for fitting

    alpha and beta if provided, will not be guessed
    """
    # convert to unitless
    ugrid, vgrid = qE2uv_grid(qgrid, Egrid)
    # initial guess of parameters
    if alpha is None:
        alpha = getAlpha(ugrid, vgrid, I, gaussian2d_threshold)
    if alpha_bounds is not None:
        # alpha (and beta) values can be 180 degree off.
        # to make interpolation easy it is better to get all alpha values in a certain range
        # if alpha_bounds is given, the code here will try to keep alpha within bounds
        alpha = adjust_to_within_bounds(alpha, alpha_bounds, period=np.pi)
    if beta is None:
        beta = getBeta(ugrid, vgrid, I)
    beta_bounds = beta_bounds or (0, np.pi)
    beta = adjust_to_within_bounds(beta, beta_bounds, period=np.pi)
    print("Guessed alpha,beta=", np.rad2deg([alpha, beta]))
    # initial guess for x and y profiles
    # x is gaussian, y is asymmetric
    xpg, yg = transform_grid(ugrid, vgrid, alpha, beta)
    # return xpg, yg
    Imin = I.max()/50
    xpmin = np.min(xpg[I>Imin]); xpmax = np.max(xpg[I>Imin])
    ymin = np.min(yg[I>Imin]); ymax = np.max(yg[I>Imin])
    # histogram I(xp)
    Nbins = 200
    Ixp, xp_bb = np.histogram(xpg.flatten(), bins=Nbins, weights=I.flatten())
    norm_xp, xp_bb = np.histogram(xpg.flatten(), bins=Nbins)
    Ixp /= norm_xp; Ixp[Ixp!=Ixp] = 0
    xp_bc = (xp_bb[1:] + xp_bb[:-1])/2
    # histogram I(y)
    Iy, y_bb = np.histogram(yg.flatten(), bins=Nbins, weights=I.flatten())
    norm_y, y_bb = np.histogram(yg.flatten(), bins=Nbins)
    Iy /= norm_y; Iy[Iy!=Iy] = 0
    y_bc = (y_bb[1:] + y_bb[:-1])/2
    # return xpg, yg, xp_bc, Ixp, y_bc, Iy
    # guess of gaussian width and center
    xp_center, xp_sigma = weighted_avg_and_std(xp_bc, Ixp)
    y_center, y_sigma = weighted_avg_and_std(y_bc, Iy)
    print("Guessed x' center, sigma: ", xp_center, xp_sigma)
    print("Guessed y center, sigma: ", y_center, y_sigma)
    return alpha, beta, xp_bc, Ixp, y_bc, Iy, xp_center, xp_sigma, y_center, y_sigma

def adjust_to_within_bounds(x, bounds, period):
    min, max = bounds
    min1 = min % period
    x1 = x % period
    xdiff = x-x1; mindiff = min-min1
    x += mindiff - xdiff
    for i in range(0,2):
        x += i*period
        if x <= max and x >= min: return x
        continue
    raise RuntimeError("failed to adjust %s to be withint %s. period=%s" % (x, bounds, period))


def test_adjust_to_within_bounds():
    assert adjust_to_within_bounds(3, (185.,365.), 180.) == 363
    assert adjust_to_within_bounds(-3, (185.,365.), 180.) == 357
    assert adjust_to_within_bounds(187, (185.,365.), 180.) == 187
    for i in range(-5, 10):
        assert adjust_to_within_bounds(187+i*180, (185.,365.), 180.) == 187
    return

# End of file 
