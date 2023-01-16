import numpy as np, histogram.hdf as hh, histogram as H
import os, sys

class Sim:
    
    scatterer_xml = "sampleassembly/V-scatterer.xml"

    def __init__(
            self,
            instrument = 'ARCS',
            workdir = ".",
            beamdir = "/SNS/users/lj7/simulations/ARCS/beam/300meV-n3e9",
            Ei = 300,
            Q = 12.,
            dQ_axis=(-1, 1, 0.02),
            dE_axis=(-60, 40, 1.),
            ncount = 1e7,
            nodes = 20
            ):
        """create a simulation working directory for simulating the point spread function
        """
        workdir = os.path.abspath(workdir)
        self.workdir = workdir
        self.beamdir = beamdir
        pwd = os.path.abspath(".")
        outdir = self.outdir = os.path.join(workdir, 'out')
        # workdir is the root of everything 
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)
        # outdir: output
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        # simdir: mcvine sim
        simdir = os.path.join(workdir, 'res-sim')
        if not os.path.exists(simdir):
            # init
            if sys.version_info < (3,0):
                cmd = u'mcvine workflow powder --instrument=%s --sample=V --workdir=%s' % (instrument, simdir)
            else:
                cmd = 'mcvine workflow powder --instrument=%s --sample=V --workdir=%s' % (instrument, simdir)
            os.system(cmd)
            os.chdir(simdir)
            os.system("rm -rf beam")
            os.system('ln -s %s beam' % beamdir)
        else:
            os.chdir(simdir)
        self.simdir = simdir
        self.Ei, self.Q = Ei, Q
        self.dQ_axis, self.dE_axis = dQ_axis, dE_axis
        self.ncount = ncount
        self.nodes = nodes
        self.instrument = instrument
        os.chdir(pwd)
        return

    def run(self, E, Q=None, dQ_axis=None, dE_axis=None):
        "run one simulation at the requested energy transfer and momentum transfer and reduce"
        self.simulate(E, Q)
        self.reduce(dQ_axis, dE_axis)
        return


    def simulate(self, E, Q=None):
        "run simulation at the requested energy transfer and momentum transfer"
        Q = self.Q if Q is None else Q
        self.E = E; self.Q = Q
        pwd = os.path.abspath(".")
        simdir = self.simdir
        os.chdir(simdir)
        # update scatterer
        text = scatterer_template.format(Q=Q, E=E)
        open(self.scatterer_xml, 'wt').write(text)
        # run
        _exec('make clean')
        def printlog(path):
            if os.path.exists(path):
                print(open(path).read())
        def onerror():
            printlog('log.scatter')
            printlog('log.create-nxs')
            raise RuntimeError("simulation failed")
        _exec(
            'time make NCOUNT={ncount} NODES={nodes}'.format(ncount=self.ncount, nodes=self.nodes),
            onerror=onerror
        )
        os.chdir(pwd)
        return


    def reduce(self, dQ_axis=None, dE_axis=None):
        "reduce simulation"
        Q, E = self.Q, self.E
        dQ_axis = self.dQ_axis if dQ_axis is None else dQ_axis
        dE_axis = self.dE_axis if dE_axis is None else dE_axis
        dQ_min,dQ_max,ddQ=dQ_axis
        dE_min,dE_max,ddE=dE_axis
        pwd = os.path.abspath(".")
        simdir = self.simdir
        os.chdir(simdir)
        cmd = ("time mcvine instruments {instrument} nxs reduce sim.nxs "
               "--qaxis {Qmin} {Qmax} {dQ} --eaxis {Emin} {Emax} {dE} --tof2E >log.reduce2 2>&1").format(
                   instrument = self.instrument.lower(),
                   Qmin=Q+dQ_min, Qmax=Q+dQ_max, dQ=ddQ,
                   Emin=E+dE_min, Emax=E+dE_max, dE=ddE)
        _exec(cmd)
        iqe = hh.load("iqe.h5")
        iqe.I[iqe.I!=iqe.I] = 0
        ie = iqe.sum('Q')
        hh.dump(ie, result_path(self.outdir, Q, E, 'I(E)'))
        hh.dump(iqe, result_path(self.outdir, Q, E, 'I(Q,E)'))
        os.chdir(pwd)
        return


kind2prefix = {
    'I(E)': 'ie',
    'I(Q,E)': 'iqe',
    }
def result_path(outdir, Q, E, kind='I(E)'):
    Q = float(Q); E = float(E)
    identifier = 'Q_%s-E_%s' % (Q,E)
    pre = kind2prefix[kind]
    return os.path.join(outdir, '%s-%s.h5' % (pre, identifier))


def list_results(outdir, kind):
    pre = kind2prefix[kind]
    pattern = os.path.join(outdir, '%s-Q_*-E_*.h5' % pre)
    import glob
    files = glob.glob(pattern)
    l = []
    for fn in files:
        b = os.path.basename(fn)
        _, Q, E = b.rstrip('.h5').split('-')
        Q, E = list(map(float, (Q[2:], E[2:])))
        l.append((Q,E))
        continue
    return l


scatterer_template = """<?xml version="1.0"?>

<!DOCTYPE scatterer>

<!-- weights: absorption, scattering, transmission -->
<homogeneous_scatterer mcweights="0, 1, 0.2">

  <ConstantQEKernel momentum-transfer="{Q}/angstrom" energy-transfer="{E}*meV">
  </ConstantQEKernel>

</homogeneous_scatterer>
"""

def _exec(cmd, onerror=None):
    if os.system(cmd):
        if onerror:
            onerror()
        else:
            raise RuntimeError("%s failed" % cmd)
    return


