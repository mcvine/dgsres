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


from .fit_2d_psf import Fit as base

class Fit(base):

    def fit(self, rounds=None):
        qgrid, Egrid = self.qEgrids
        reshist = self.mcvine_psf_qE
        res_x = reshist.q
        res_y = reshist.E
        res_z = reshist.I.T.copy()
        dx = res_x[1]-res_x[0]
        dy = res_y[1]-res_y[0]
        res_z /= np.sum(res_z)*dx*dy
        self.res_z = res_z
        return fit(qgrid, Egrid, res_z, rounds=rounds)


    def createModel(self):
        # model is in the ellipsoid method below
        return
    
    
def fit(qgrid, Egrid, I, rounds=None):
    if not rounds: rounds = 3
    # convert to unitless
    qrange = qgrid[0][-1] - qgrid[0][0]
    Erange = Egrid[:, 0][-1] - Egrid[:, 0][0]
    ugrid = qgrid/qrange; vgrid = Egrid/Erange
    alpha, beta, xp_bc, Ixp, y_bc, Iy, xp_center, xp_sigma, y_center, y_sigma = fitguess(qgrid, Egrid, I)
    z = I.copy().flatten()/I.max()
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
    print "Established model:", model
    model.print_param_hints(colwidth=12)
    print "Start fitting..."
    results = []
    for i in range(rounds):
        print " -- Fitting round %s" % i
        result = model.fit(z, u=ugrid.flatten(), v=vgrid.flatten(), method='differential_evolution')
        # print result.fit_report()
        results.append(result)
        print "    chisq=%s" % result.chisqr
        # print
    results.sort(key=lambda x: x.chisqr)
    return results[0]


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

def getAlpha(ugrid, vgrid, I):
    "guess alpha angle. the bright part"
    # median = np.nanmedian(I[I>0])
    max = np.nanmax(I)
    # large=np.where(I>median)
    large = np.where(I>max*.5)
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

def fitguess(qgrid, Egrid, I):
    "return guess parameters for fitting"
    # convert to unitless
    qrange = qgrid[0][-1] - qgrid[0][0]
    Erange = Egrid[:, 0][-1] - Egrid[:, 0][0]
    ugrid = qgrid/qrange; vgrid = Egrid/Erange
    # initial guess of parameters
    alpha = getAlpha(ugrid, vgrid, I)
    beta = getBeta(ugrid, vgrid, I)
    print "Guessed alpha,beta=", np.rad2deg([alpha, beta])
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
    Ixp /= norm_xp
    xp_bc = (xp_bb[1:] + xp_bb[:-1])/2
    # histogram I(y)
    Iy, y_bb = np.histogram(yg.flatten(), bins=Nbins, weights=I.flatten())
    norm_y, y_bb = np.histogram(yg.flatten(), bins=Nbins)
    Iy /= norm_y
    y_bc = (y_bb[1:] + y_bb[:-1])/2
    # return xpg, yg, xp_bc, Ixp, y_bc, Iy
    # guess of gaussian width and center
    xp_center, xp_sigma = weighted_avg_and_std(xp_bc, Ixp)
    y_center, y_sigma = weighted_avg_and_std(y_bc, Iy)
    print "Guessed x' center, sigma: ", xp_center, xp_sigma
    print "Guessed y center, sigma: ", y_center, y_sigma
    return alpha, beta, xp_bc, Ixp, y_bc, Iy, xp_center, xp_sigma, y_center, y_sigma

# End of file 