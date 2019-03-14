"""
see
* https://jupyter.sns.gov/user/{USER}/notebooks/data/SNS/SEQ/IPTS-16800/shared/resolution/resolution%20simulations%20-%20improve%20workflow.ipynb
* https://jupyter.sns.gov/user/lj7/notebooks/data/SNS/SEQ/IPTS-16800/shared/resolution/resolution%20fit%20-%20improve%20workflow.ipynb
"""

import os, shutil, subprocess as sp, time
import numpy as np
from matplotlib import pyplot as plt
from . import use_res_comps, fit_ellipsoid

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
    for q in slice.grid.qaxis.ticks():
        for E in slice.grid.Eaxis.ticks():
            simdir = config.simdir(q,E, slice)
            try:
                outputs[(q,E)] = simulate(q=q, E=E, slice=slice, outdir=simdir, config=config, Nrounds_beam=Nrounds_beam)
            except:
                failed.append( (q,E) )
        continue
    return outputs, failed


def plot_resolution_on_grid(slice, config, figsize=(10, 7)):
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
            simdir = config.simdir(q,E, slice)
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


def fit(q, E, slice, config):
    from dgsres.singlextal import fit_ellipsoid
    datadir = config.simdir(q,E,slice)
    qaxis = slice.res_2d_grid.qaxis
    Eaxis = slice.res_2d_grid.Eaxis
    fitter = fit_ellipsoid.Fit(
        datadir,
        qaxis=(qaxis.min, qaxis.max, qaxis.step),
        Eaxis=(Eaxis.min, Eaxis.max, Eaxis.step),
        Ei=config.Ei, E=E
    )
    fitter.load_mcvine_psf_qE(adjust_energy_center=True)
    fitting_params = dict([(k,v) for k,v in slice.fitting.__dict__.items() if not k.startswith('_')])
    fitter.fit_result = fitter.fit(**fitting_params)
    return fitter


def plot_compare_fit_to_data(fitter, figsize=(8,4)):
    res_z = fitter.res_z
    qgrid, Egrid = fitter.qEgrids
    result = fitter.fit_result
    plt.figure(figsize=figsize)
    plt.subplot(1,2,1)
    plt.pcolormesh(qgrid, Egrid, res_z)
    plt.colorbar()
    plt.subplot(1,2,2)
    scale = res_z.sum()/result.best_fit.sum()
    plt.pcolormesh(qgrid, Egrid, result.best_fit.reshape(res_z.shape)*scale)
    plt.colorbar()
    return

def fit_all_grid_points(slice, config):
    qE2fitter = dict()
    nofit = []
    for q in slice.grid.qaxis.ticks():
        for E in slice.grid.Eaxis.ticks():
            print q, E
            try:
                qE2fitter[(q,E)] = fit(q, E, slice, config)
            except:
                nofit.append((q,E))
        continue
    # qEranges is an import parameter that need to be rememberd
    fitter1 = qE2fitter.values()[0]
    slice.res_2d_grid.qEranges = fit_ellipsoid.qEgrid2range(*fitter1.qEgrids)
    for fitter in qE2fitter.values()[1:]:
        assert slice.res_2d_grid.qEranges == fit_ellipsoid.qEgrid2range(*fitter.qEgrids)
    return qE2fitter, nofit

def plot_resfits_on_grid(qE2fitter, slice, config, figsize=(10, 7)):
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
            fitter = qE2fitter.get((q,E))
            if fitter is None: continue
            dqgrid, dEgrid = fitter.qEgrids
            result = fitter.fit_result
            ax1 = axes[irow][icol]
            ax1.set_title("q=%.2f, E=%.2f" % (q, E))
            ax1.pcolormesh(dqgrid, dEgrid, result.best_fit.reshape(dqgrid.shape))
    plt.tight_layout()
    return

def save_fits_as_pickle(qE2fitter, path):
    qE2fitres_tosave = dict()

    for qe in qE2fitter.keys():
        fr = qE2fitter[qe].fit_result
        fr_tosave = fit_ellipsoid.FitResult()
        fr_tosave.best_values = fr.best_values
        qE2fitres_tosave[qe] = fr_tosave
        continue

    import pickle as pkl
    pkl.dump(qE2fitres_tosave, open(path, 'w'))
    return

def print_parameter_table(qE2fitres):
    keys = qE2fitres.values()[0].best_values.keys()
    print "%5s %5s" % ('q','E'),
    for k in keys: print ('%8s' % k[:8]),
    print
    qEs = qE2fitres.keys()
    for q, E in qEs:
        if not (q,E) in qE2fitres: continue
        result = qE2fitres[(q,E)]
        print "%5.1f %5.1f" % (q,E),
        for k in keys:
            v = result.best_values[k]
            print ('%8.4f' % v),
        print
    return


def createInterpModel(qE2fitres, slice):
    # Get parameters as lists, ready for interpolation
    keys = qE2fitres.values()[0].best_values.keys()
    qEs_all = qE2fitres.keys()
    qE_points = []
    param_values = dict()
    for k in keys:
        param_values[k] = []
    for q,E in qEs_all:
        if (q,E) not in qE2fitres: continue
        result = qE2fitres[(q,E)]
        bv = result.best_values
        qE_points.append((q,E))
        for k in keys:
            vals = param_values[k]
            v = bv[k]
            vals.append(v)
        continue
    qrange, Erange = slice.res_2d_grid.qEranges
    return fit_ellipsoid.InterpModel(qE_points, param_values, qrange, Erange)
