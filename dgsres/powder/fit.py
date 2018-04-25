import numpy as np
import histogram.hdf as hh, histogram as H
import lmfit
from dgsres import icg
from .use_ConstantQEKernel import result_path, list_results

class Fit_IE:

    def __init__(self, simout, Ei, params=None, geom=None):
        self.simout = simout
        self.Ei = Ei
        # fitting parameters
        if params is None:
            params = lmfit.Parameters()
            params.add('a', min=0., max=1.)
            params.add('b', min=0., max=.3)
            params.add('R', value=0.3, vary=False)
            params.add('sigma', min=0., max=20.)
            params.add('t0', min=0., max=100.)
        self.params = params
        self.geom = geom
        return

    def getData(self, Q, E):
        simout = self.simout
        # get mcvine sim result
        res_E = hh.load(result_path(simout, Q, E, 'I(E)'))
        x = res_E.E
        y0 = res_E.I
        scale = y0.sum()
        y0/=scale
        y0err = res_E.E2**.5
        y0err/=scale
        return x, y0, y0err
    
    def fit(self, Q, E):
        # get mcvine sim result
        x, y0, y0err = self.getData(Q,E)
        # prepare fitting data
        data = y0
        eps_data = 1.
        #
        geom = self.geom
        # fit
        out = lmfit.minimize(
            residual, self.params,
            args=(x, data, eps_data, geom, self.Ei, E),
            method='differential_evolution')
        return out

    def fit_all(self):
        QE_pairs = list_results(self.simout, 'I(E)')
        results = {}
        for Q,E in QE_pairs:
            results[(Q, E)] = self.fit(Q, E)
            continue
        return results

    pass

def residual(params, x, data, eps_data, geom, Ei, E):
    a = params['a']
    b = params['b']
    R = params['R']
    sigma = params['sigma']
    t0 = params['t0']
    model = icg.resolution(x, Ei=Ei, E0=E, a=a, b=b, R=R, sigma=sigma, t0=t0, geom=geom)
    sum = model.sum()
    if np.abs(sum)<1e-10: sum = 1e-10
    model/=sum
    model[model!=model] = np.nanmax(model)
    # plt.plot(data, )
    # plt.plot(model)
    ret = (data-model) / eps_data
    # print np.sum(ret**2)/ret.size
    return ret


