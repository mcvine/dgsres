import os, numpy as np
from . import use_covmat
from .pointcloud import PointCloud

class VioliniModel:

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
    
