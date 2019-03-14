import os, numpy as np
from . import use_covmat

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


class PointCloud:
    
    """resolution point cloud. can be regarded as "events" with dh,dk,dl,dE,weight data.
    These data could be generated from mcvine simulation of the resolution function,
    or generated from a simple normal distribution using a 4D cov matrix (see AnalyticalModel).
    """
    
    def __init__(self, dhs, dks, dls, dEs, weights):
        self.dhs = dhs
        self.dks = dks
        self.dls = dls
        self.dEs = dEs
        self.weights = weights
        return
    
    def getThinSlice(
        self, 
        axis1=('h', -0.1,0.1, 0.002), axis2=('E', -3, 3., 0.04),
        axis3=('k', -0.006, 0.006), axis4=('l', -0.006, 0.006)
    ):
        axis1_name = axis1[0]; axis1_ticks = np.arange(*axis1[1:])
        axis2_name = axis2[0]; axis2_ticks = np.arange(*axis2[1:])
        condition = True
        for ax in [axis3, axis4]:
            name = ax[0]
            min, max = ax[1:]
            arr = self._evts(name)
            condition *= (arr<max) * (arr>min)
            continue
        Ixy, xedges, yedges = np.histogram2d(
            self._evts(axis1_name)[condition], self._evts(axis2_name)[condition],
            bins=(axis1_ticks, axis2_ticks), weights=self.weights[condition]) 
        xbc = (xedges[:-1] + xedges[1:])/2
        ybc = (yedges[:-1] + yedges[1:])/2
        xg, yg = np.meshgrid(xbc, ybc)
        return xg,yg,Ixy
            
    def _evts(self, name):
        return getattr(self, 'd%ss' % name)
    

class McvineResolutionData:

    """resolution data simulated by mcvine. 
    The data was simulated earlier and saved to disk. 
    This class handles locating the data files and reading the data, 
    and convert the data to different forms: point cloud, cov matrix, etc.
    """
    
    def __init__(self, parent_dir, dirname_template='E%s_hkl%s'):
        self.parent_dir = parent_dir
        self.dirname_template = dirname_template
        return
    
    def path(self, hkl, E):
        return os.path.join(self.parent_dir, self.dirname_template % (E, '%s,%s,%s' % tuple(hkl)))
    
    def loadData(self, hkl, E):
        p = self.path(hkl, E)
        return self._loadData(p)
    
    def loadPointCloud(self, hkl, E):
        dhs, dks, dls, dEs, probs = self.loadData(hkl, E)
        return PointCloud(dhs, dks, dls, dEs, probs)
    
    def _loadData(self, outdir1):
        dhkls = np.load('%s/dhkls.npy' % outdir1)
        dEs = np.load('%s/dEs.npy' % outdir1)
        probs = np.load('%s/probs.npy' % outdir1)
        dhs,dks,dls = dhkls.T
        # there might be unreasonable data points with unreasonable weights
        mask = (dhs> -2.)*(dhs<2.) \
            * (dks> -2.)*(dks<2.) \
            * (dls> -2.)*(dls<2.) 
        return np.array([dhs[mask], dks[mask], dls[mask], dEs[mask], probs[mask]])
    
    def computeCovMat(self, hkl, E):
        data = self.loadData(hkl, E)
        Data = data[:4]; probs = data[-1]
        return np.cov(Data, aweights=probs)

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

