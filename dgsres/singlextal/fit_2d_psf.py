# fit mcvine-simulated point spread function to an empirical form

import os, numpy as np
import histogram.hdf as hh, histogram as H
import scipy.optimize as sopt


class Fit(object):

    def __init__(self, sim_path, qaxis=(-0.3, 0.3, 0.01), Eaxis=(-2, 1, 0.05), Ei=None, E=None):
        self.sim_path = sim_path
        self.qaxis = qaxis
        self.Eaxis = Eaxis
        self.Ei = Ei
        self.E = E
        self.createModel()
        return

    def createModel():
        raise NotImplementedError("createModel")
        self.model = None
        return

    def load_mcvine_psf_qE(self):
        path = self.sim_path
        qaxis = self.qaxis; Eaxis = self.Eaxis
        dEs = np.load(os.path.join(path, 'dEs.npy'))
        dqs = np.load(os.path.join(path, 'dxs.npy'))
        probs = np.load(os.path.join(path, 'probs.npy'))
        hist, qedges, Eedges = np.histogram2d(dqs, dEs, bins=(np.arange(*qaxis), np.arange(*Eaxis)), weights=probs)
        qaxis = H.axis('q', boundaries=qedges)
        Eaxis = H.axis('E', boundaries=Eedges)
        # remove outliers
        median = np.median(hist[hist>0])
        normal = hist[hist<median*100]; normal_max = np.max(normal)
        hist[hist>median*100] = normal_max
        # res(q, E) histogram
        self.mcvine_psf_qE = reshist  = H.histogram('res', (qaxis, Eaxis), data=hist)
        self.qEgrids = np.meshgrid(reshist.q, reshist.E)
        # res(E) histogram
        self.mcvine_psf_E = reshist.sum('q')
        # res(q) histogram
        self.mcvine_psf_q = reshist.sum('E')
        return

    def fit_E_profile(self, a=None, b=None, R=None, sigma=None, t0=None):
        raise NotImplementedError("fit_E_profile")

    def fit_q_profile(self, sigma=None):
        res_q = self.mcvine_psf_q
        res_x = res_q.q.copy()
        res_y = res_q.I.copy()
        dx = res_x[1] - res_x[0];  res_y /= res_y.sum()*dx
        p0 = [1./sigma, 0., sigma]
        popt,pcov = sopt.curve_fit(gaus,res_x,res_y,p0=p0)
        q_profile = self.model.set_q_profile(popt[-1])
        yfit = gaus(res_x, *popt)
        return popt, res_x, res_y, yfit, q_profile

    def fit_qE_profile(self, dq_over_dE):
        reshist = self.mcvine_psf_qE
        res_x = reshist.q
        res_y = reshist.E
        res_z = reshist.I.T.copy()
        dx = res_x[1]-res_x[0]
        dy = res_y[1]-res_y[0]
        res_z /= np.sum(res_z)*dx*dy
        # dq_over_dE could be a single value or a tuple (bounds)
        try:
            # bound
            l, u = dq_over_dE
            input_type = 'bounds'
        except:
            input_type = 'single value'
            pass

        if input_type == 'bounds':
            bounds = np.array([[l,u]])
            x0 = ((l+u)/2., )
        else:
            # single value
            x0 =  (dq_over_dE,)
            if dq_over_dE>0:
                bounds = np.array(
                    [(0.,),
                     (dq_over_dE*3,)]
                ).T
            elif dq_over_dE < 0:
                bounds = np.array(
                    [(dq_over_dE*3,),
                     (0.,),]
                ).T
            else:
                raise ValueError("dq_over_dE estimate must not be zero")
        cost = self.model.make_qE_profile_cost(res_x, res_y, res_z)
        fitres = sopt.minimize(
            cost, x0=x0, bounds=bounds)#, method='L-BFGS-B')
        qE_profile = self.model.set_qE_profile(fitres.x[0])
        qgrid,Egrid = np.meshgrid(res_x, res_y)
        zfit = qE_profile(qgrid, Egrid); zfit /= np.sum(zfit)*dx*dy
        return fitres.x, qgrid, Egrid, res_z, zfit, qE_profile


class PSF_Affine_Model(object):

    def __init__(self, Ei, E):
        self.Ei = Ei
        self.E = E
        self.E_profile = None
        self.q_profile = None
        return

    def set_qE_profile(self, dq_over_dE):
        self.qE_profile = self.make_qE_profile_function(dq_over_dE)
        return self.qE_profile

    def make_qE_profile_cost(self, res_x, res_y, res_z):
        xgrid,ygrid = np.meshgrid(res_x, res_y)
        dx = res_x[1]-res_x[0]
        dy = res_y[1]-res_y[0]
        def cost(x):
            dq_over_dE, = x
            f = self.make_qE_profile_function(dq_over_dE)
            zfit = f(xgrid, ygrid)
            # print "sum of area", zfit.sum()*dx*dy
            zfit /= zfit.sum()*dx*dy
            c =  np.sum((zfit-res_z)**2)/res_z.size
            # print "cost", c
            return c
        return cost

    def make_qE_profile_function(self, dq_over_dE):
        def _(qgrid, Egrid):
            return ellipsoid_affine(qgrid, Egrid, dq_over_dE, self.q_profile, self.E_profile)
        return _

    def set_E_profile(self, *args, **kwds):
        raise NotImplementedError("set_E_profile")

    def set_q_profile(self, sigma):
        def _(x):
            h = 1./np.sqrt(2*np.pi)/np.abs(sigma)
            return gaus(x, h, 0., sigma)
        self.q_profile = _
        return self.q_profile
    
    def make_E_profile_function(self, a, b, R, sigma=None, t0=None):
        raise NotImplementedError("make_E_profile_function")

    def make_E_profile_cost(self, res_x, res_y):
        def cost(x):
            f = self.make_E_profile_function(*x)
            yfit = f(res_x)
            dx = res_x[1]-res_x[0]
            yfit /= np.sum(yfit) * dx
            return np.sum((yfit-res_y)**2)/res_y.size
        return cost


def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def ellipsoid_affine(q, E, dq_over_dE, q_profile, E_profile):
    """q and E are grids
    I(q) are gaussians. their centers are linearly related to E
    """
    pE = E_profile(E)
    q1 = E*dq_over_dE
    pq = q_profile(q-q1)
    return pq*pE

