"""
see https://jupyter.sns.gov/user/{USER}/notebooks/data/SNS/SEQ/IPTS-16800/shared/resolution/resolution%20simulations%20-%20improve%20workflow.ipynb#
"""

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


def simulate_all_grid_points(slice, config, Nrounds_beam=1):
    failed = []; outputs = {}
    outdir_pattern = config.simdir_pattern
    for q in slice.grid.qaxis.ticks():
        for E in slice.grid.Eaxis.ticks():
            simdir = outdir_pattern % (q,E)
            try:
                outputs[(q,E)] = simulate(q=q, E=E, slice=slice, outdir=simdir, config=config, Nrounds_beam=Nrounds_beam)
            except:
                failed.append( (q,E) )
        continue
    return outputs, failed


def plot_resolution_on_grid(slice, config, figsize=(10, 7)):
    from matplotlib import pyplot as plt
    qs = slice.grid.qaxis.ticks()
    Es = slice.grid.Eaxis.ticks()
    ncols = len(qs)
    nrows = len(Es)
    res_2d_grid = slice.res_2d_grid
    hkl_projection = slice.hkl_projection
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    for irow in range(nrows):
        for icol in range(ncols):
            q = qs[icol]
            E = Es[irow]
            simdir = config.simdir_pattern % (q,E)
            try:
                probs = np.load('%s/probs.npy' % simdir)
            except IOError:
                continue
            dEs = np.load('%s/dEs.npy' % simdir)
            dhkls = np.load('%s/dhkls.npy' % simdir)
            dqs = np.dot(dhkls, hkl_projection)/np.linalg.norm(hkl_projection)**2
            I, dqedges, dEedges = np.histogram2d(
                bins=(res_2d_grid.qaxis.ticks(), res_2d_grid.Eaxis.ticks()), weights=probs, x=dqs, y=dEs )
            dqg, dEg = np.meshgrid(dqedges, dEedges)
            ax1 = axes[irow][icol]
            median = np.nanmedian(I[I>0])
            good = I[I<median*200]
            goodmax = good.max()
            I[I>median*100] = goodmax
            ax1.set_title("q=%.2f, E=%.2f" % (q, E))
            ax1.pcolormesh(dqg, dEg, I.T)
    plt.tight_layout()
    return
