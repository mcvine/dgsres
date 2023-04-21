import os


class config_cls(object):

    def __init__(self, Ei, beampth, sim_Nrounds_beam=1):
        self.slices = []
        self.Ei = Ei
        self.beam = beampth
        self.sim_Nrounds_beam = sim_Nrounds_beam
        self.thisdir = os.path.abspath(os.path.dirname(__file__) or '.')
        self.psi_scan = None
        self.sample_yaml = None
        self.fitting = None
        self.instrument = None

    def simdir(self, q, E, slice):
        return 'sim-{}-q_{:0.3f},E_{:0.3f}'.format(slice.name, q, E)