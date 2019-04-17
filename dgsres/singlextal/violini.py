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
    
    def pdf4D(self, hkl, E):
        "return a pdf function pdf((h,k,l,E)) gives value of pdf at that point"
        sigma2 = self.computeCovMat(hkl, E)
        return create_pdf(sigma2)

    def pdf2D(self, hkl, E, eu, ev, ew):
        """return a 2D pdf function along uaxis and Eaxis at the point h,k,l,E
        the pdf is 4D pdf integrated over v and w direction

        eu, ev, ew are unit vectors along uvw

        see https://jupyter.sns.gov/user/{UID}/notebooks/data/SNS/SEQ/IPTS-21411/shared/resolution/violini-04162019.ipynb
        for an example
        """
        # compute covmat
        cm = self.computeCovMat(hkl, E)
        cm = np.matrix(cm)
        # inverse
        ic = cm.I
        # transform the coordinate system
        A = np.array( (eu,ev,ew) ).T
        A = np.block( [[A, np.zeros((3,1))],
                       [np.zeros((1,3)), 1]])
        A = np.matrix(A)
        cm0 = cm
        # print cm0
        cm = A.T * cm * A # u,v,w,E
        # print cm
        # integration along v and w
        cm2d = [[cm[0,0] - (cm[0,1]+cm[1,0])**2/4./cm[1,1] - (cm[0,2]+cm[2,0])**2/4./cm[2,2],
                 cm[0,3]],
                [cm[3,0],
                 cm[3,3] - (cm[3,1]+cm[1,3])**2/4./cm[1,1] - (cm[3,2]+cm[2,3])**2/4./cm[2,2]
                ]]
        # print cm2d
        return create_pdf(cm2d)


def create_pdf(covmat):
    "create probability distribution function out of covariance matrix"
    sigma2 = np.matrix(covmat)
    det = np.linalg.det(sigma2)
    # det = np.abs(det)
    if det == 0:
        raise ValueError("The covariance matrix can't be singular")
    size, size2 = sigma2.shape
    assert size == size2
    # print "det", det
    norm_const = 1.0/ ( np.power((2*np.pi),float(size)/2) * np.power(det,1.0/2) )
    # print "norm_const", norm_const
    inv = sigma2.I
    # print "inv", inv
    def _(x):
        x = np.matrix(x)
        result = np.exp( -0.5 * (x * inv * x.T) )
        return norm_const * result
    return _
