# helper classes and methods for interpolation

import numpy as np

class InterpedParameters:

    def __init__(self, coordinates2parameters, **interp_options):
        """
        * coordinates2parameters: map coordinates to parameters

            For example:

                {(-2.0, -6.0): {'dq_over_dE': 0.041900057195662116,
                  'ef_width': 0.51000000000000001,
                  'q_sigma': 0.027102026809677951,
                  'sigma_left': 0.5,
                  'sigma_right': 0.5,
                  'weight_left': 0.51261838147713479},
                  ...
                }

        * interp_options: interpolation options. for example: kind='linear'
        """
        self.c2p = c2p = coordinates2parameters

        self.coordinates = keys = c2p.keys()
        k0 = keys[0]
        assert len(k0) == 2, "Only support 2D"

        # columns of coordinates 
        xs, ys = np.array(keys).T

        #
        v0 = c2p[k0]
        interps = dict()
        self.param_names = v0.keys()
        import scipy.interpolate as si
        for name in self.param_names:
            interps[name] = si.interp2d(xs, ys, self.build_z(name), **interp_options)
            continue
        self.interps = interps
        return


    def build_z(self, name):
        return [self.c2p[c][name] for c in self.coordinates]


    def at(self, *x):
        v = dict()
        for name in self.param_names:
            v[name] = self.interps[name](*x)
            continue
        return v
