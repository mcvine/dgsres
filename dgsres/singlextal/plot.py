import os, numpy as np
from . import use_covmat
from .pointcloud import PointCloud
from .simdata import McvineResolutionData
from .violini import VioliniModel
# backward compatible
AnalyticalModel = VioliniModel

# obsolete
def createARCSAnalyticalModel(
    tau_P, tau_M, pix_r, pix_h,
    sample_thickness, 
    sigma_thetai, sigma_phii,
    Ei, sample_yml, psi_scan
):
    """create analytical resolution model for ARCS. cf. Violini et al.
    """
    tofwidths = use_covmat.tofwidths(P=tau_P, M=tau_M)
    beamdivs = use_covmat.beamdivs(theta=sigma_thetai, phi=sigma_phii)
    samplethickness = sample_thickness
    instrument = use_covmat.instrument(
        name = 'ARCS',
        detsys_radius = "3.*meter",
        L_m2s = "13.6*meter",
        L_m2fc = "11.61*meter",
        offset_sample2beam = "-0.15*meter" # offset from sample to saved beam
        )
    pixel = use_covmat.pixel(
        radius = "%s*inch" % pix_r,
        height = "meter*%s" % pix_h,
        pressure = "10*atm",
        )
    return VioliniModel(
        instrument, pixel, tofwidths, beamdivs, 
        sample_yml, samplethickness,
        Ei, psi_scan)


def plot_qE_ellipse(covmat, q, symbol='.', **kwds):
    """plot 2d ellipsoid along a particular hkl direction, given the cov matrix.
    """
    invcm = np.linalg.inv(covmat)
    _, u = compute_qE_ellipse(invcm, q)
    from matplotlib import pyplot as plt
    plt.plot(u[:,0], u[:,1], symbol, **kwds)
    return


def compute_qE_ellipse(InvCov4D, q):
    "compute ellipsoid along a q direction"
    qE2qE = np.array(
        [np.hstack([q, [0]]),
         [0,0,0,1]])
    inv_cov_qE = np.dot(qE2qE, np.dot(InvCov4D, qE2qE.T))
    # print inv_cov_hE
    r = np.linalg.eig(inv_cov_qE)
    mR = r[1]; lambdas = r[0]
    RR = 2*np.log(2)
    theta = np.arange(0, 360, 1.)*np.pi/180
    u1p = np.sqrt(RR/lambdas[0])*np.cos(theta)
    u2p = np.sqrt(RR/lambdas[1])*np.sin(theta)
    up = np.array([u1p, u2p]).T
    u = np.dot(up, mR.T)
    return inv_cov_qE, u

def compute_qE_ellipses(InvCov4D, directions=None):
    if directions is None:
        directions = np.eye(3, dtype=float)
    return [(q, compute_qE_ellipse(InvCov4D, q)) for q in directions]

