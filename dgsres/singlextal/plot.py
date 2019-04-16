import os, numpy as np
from . import use_covmat
from .pointcloud import PointCloud
from .simdata import McvineResolutionData

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
    return AnalyticalModel(
        instrument, pixel, tofwidths, beamdivs, 
        sample_yml, samplethickness,
        Ei, psi_scan)


class AnalyticalModel:

    """analytical resolution model based on paper by Violini et al.
    The main calculation is done in module .use_covmat.
    """
    
    def __init__(
        self, instrument, pixel, tofwidths, beamdivs, 
        sample_yml, samplethickness,
        Ei, psi_scan):
        self.instrument = instrument
        self.pixel = pixel
        self.tofwidths = tofwidths
        self.beamdivs = beamdivs
        self.sample_yml = sample_yml
        self.samplethickness = samplethickness
        self.Ei = Ei
        self.psi_scan = psi_scan
        return
    
    def computePointCloud(self, hkl, E, N=int(1e6)):
        covmat = self.computeCovMat(hkl, E)
        events = np.random.multivariate_normal(np.zeros(4), covmat, size=N)
        dhs, dks, dls, dEs = events.T
        ws = np.ones(dhs.shape)
        return PointCloud(dhs, dks, dls, dEs, ws)
    
    def computeCovMat(self, hkl, E):
        class dynamics:
            hkl_dir = np.array([1.,0.,0.])
            dq = 0
        dynamics.hkl0 = hkl
        dynamics.E = E
        cm_res = use_covmat.compute(
            self.sample_yml, self.Ei, 
            dynamics, 
            self.psi_scan,
            self.instrument, self.pixel,
            self.tofwidths, self.beamdivs, self.samplethickness,
            plot=False)
        # ellipsoid_trace = cm_res['u']
        InvCov4D = cm_res['hklE_inv_cov']
        return np.linalg.inv(InvCov4D)/2.355
    
def plotEllipsoid(covmat, q, symbol='.'):
    """plot 2d ellipsoid along a particular hkl direction, given the cov matrix.
    """
    invcm = np.linalg.inv(covmat)
    _, u = computeEllipsoid(invcm, q)
    from matplotlib import pyplot as plt
    plt.plot(u[:,0], u[:,1], symbol)
    return

def computeEllipsoid(InvCov4D, q):
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

def computeEllipsoids(InvCov4D, directions=None):
    if directions is None:
        directions = np.eye(3, dtype=float)
    return [(q, computeEllipsoid(InvCov4D, q)) for q in directions]

