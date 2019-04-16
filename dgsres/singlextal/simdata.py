import os, numpy as np
from .pointcloud import PointCloud

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

