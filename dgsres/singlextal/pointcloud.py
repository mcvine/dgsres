import os, numpy as np

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
    

