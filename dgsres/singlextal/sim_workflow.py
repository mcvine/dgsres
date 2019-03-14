import os, shutil, subprocess as sp, time
import numpy as np
from . import use_res_comps

def plotDynRange(hkl0, hkl_projection, qaxis, Erange, config):
    from mcvine.workflow.singlextal import dynrange
    from mcvine.workflow.sample import loadSampleYml

    sample = loadSampleYml(config.sample_yaml)
    psi_scan = config.psi_scan
    psilist = np.arange(psi_scan.min, psi_scan.max, psi_scan.step)

    dynrange.plotDynRangeOfSlice(
        sample, psi_scan.ticks(), config.Ei, hkl0, hkl_projection, qaxis,
        config.instrument.scattering_angle_constraints,
        Erange=Erange)
    return 


def simulate(q, E, slice, outdir, config, Nrounds_beam=1):
    hkl0 = slice.hkl0
    hkl_projection = slice.hkl_projection
    hkl = hkl0 + hkl_projection*q
    # setup
    use_res_comps.setup(
        outdir,
        config.sample_yaml, config.beam, E, hkl, hkl_projection,
        config.psi_scan, config.instrument.instrument, config.instrument.pixel)

    # more configuration
    open(os.path.join(outdir, 'mc_params.yml'), 'wt').write("Nrounds_beam: %s"%Nrounds_beam)

    # run
    cmd = "python run.py"
    start = time.time()
    out = sp.check_output(cmd, shell=True, cwd=outdir)
    end = time.time()
    duration = end - start
    print "* simulation took %s seconds" % duration
    return out
