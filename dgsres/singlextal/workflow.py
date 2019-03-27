"""
see
* https://jupyter.sns.gov/user/{USER}/notebooks/data/SNS/SEQ/IPTS-16800/shared/resolution/resolution%20simulations%20-%20improve%20workflow.ipynb
* https://jupyter.sns.gov/user/lj7/notebooks/data/SNS/SEQ/IPTS-16800/shared/resolution/resolution%20fit%20-%20improve%20workflow.ipynb
"""

import os, sys, shutil, subprocess as sp, time
import numpy as np
from matplotlib import pyplot as plt
from . import use_res_comps, fit_ellipsoid, _workflow_pdf_helpers as _wph
from mcvine.workflow import singlextal as sx


def simulate_all_in_one(config):
    "simulate all grid points, and compose PDF reports"
    import pylatex
    Ei = config.Ei
    Erange = (-0.3*Ei, .95*Ei)
    for sl in config.slices:
        doc = _wph.initReportDoc("%s-sim-report" % sl.name) # report document
        # info
        _wph.slice_info_section(sl, doc)
        
        qaxis = sl.grid.qaxis; Eaxis = sl.grid.Eaxis

        # dyn range plot
        # larger q range for a broader view 
        ratio = 1.
        expanded_qaxis = sx.axis(
            min=qaxis.min-(qaxis.max-qaxis.min)*ratio/2,
            max=qaxis.max+(qaxis.max-qaxis.min)*ratio/2,
            step=qaxis.step
        ).ticks()
        width = r'1\textwidth'
        with doc.create(pylatex.Section('Dynamical range')):
            with doc.create(pylatex.Figure(position='htbp')) as plot:
                plt.figure()
                plotDynRange(
                    sl.hkl0, sl.hkl_projection,
                    qaxis= expanded_qaxis, Erange=Erange,
                    config=config)
                plot.add_plot(width=pylatex.NoEscape(width))
                plot.add_caption('Dynamical range for slice %s' % sl.name)
                plt.close()

        # simulate
        with doc.create(pylatex.Section('Simulated resolution functions on a grid')):
            outputs, failed = simulate_all_grid_points(
                slice=sl, config=config, Nrounds_beam=config.sim_Nrounds_beam, overwrite=False)

            if failed:
                # this seems unecessary as what is missing is clear in the plot
                """
                doc.append("Failed to calculate resolution functions for the following (Q,E) pairs:")
                with doc.create(pylatex.Itemize()) as itemize:
                    for f in failed:
                        itemize.add_item(str(f))
                """
                pass
            # plot
            with doc.create(pylatex.Figure(position='htbp')) as plot:
                plt.figure()
                plot_resolution_on_grid(sl, config, figsize=(10, 10))
                plot.add_plot(width=pylatex.NoEscape(width))
                plot.add_caption('Simulated resolution functions for %s' % sl.name)
                plt.close()
        # save pdf
        doc.generate_pdf(clean_tex=False)
        continue
    return

def fit_all_in_one(config):
    "fit all grid points, and compose PDF reports"
    import pylatex, dill
    Ei = config.Ei
    Erange = (-0.3*Ei, .95*Ei)
    width = r'1\textwidth'
    for sl in config.slices:
        doc = _wph.initReportDoc("%s-fit-report" % sl.name) # report document
        qaxis = sl.grid.qaxis; Eaxis = sl.grid.Eaxis
        # info
        _wph.slice_info_section(sl, doc)
        
        # fit
        with doc.create(pylatex.Section('Fit resolution functions on grid')):
            # path to saved result
            path = '%s-fit_all_grid_points.dill' % sl.name
            if os.path.exists(path):
                qE2fitter, nofit = dill.load(open(path))
            else:
                qE2fitter, nofit = fit_all_grid_points(sl, config, use_cache=True)
                dill.dump((qE2fitter, nofit), open(path, 'w'), recurse=True)
            # plot
            with doc.create(pylatex.Figure(position='htbp')) as plot:
                plt.figure()
                plot_resfits_on_grid(qE2fitter, sl, config, figsize=(10,10))
                plot.add_plot(width=pylatex.NoEscape(width))
                plot.add_caption('Fitted resolution functions for %s' % sl.name)
                plt.close()
            doc.append(pylatex.utils.NoEscape(r"\clearpage"))
        # save
        pklfile = '%s-fit_results.pkl' % sl.name
        save_fits_as_pickle(qE2fitter, pklfile)
        import pickle as pkl
        qE2fitres = pkl.load(open(pklfile))
        
        # parameters
        with doc.create(pylatex.Subsection('Fitted parameters')):
            s = format_parameter_table(qE2fitres)
            doc.append(_wph.verbatim(s))

        # interpolated model
        with doc.create(pylatex.Subsection('Interpolated model')):
            imodel = get_interped_resolution_model(sl)
            qs = (qaxis.ticks() + qaxis.step/2.)[:-1]
            Es = (Eaxis.ticks() + Eaxis.step/2.)[:-1]
            dqgrid, dEgrid = qE2fitter.values()[0].qEgrids
            # plot
            with doc.create(pylatex.Figure(position='htbp')) as plot:
                plt.figure()
                plot_interpolated_resolution_on_grid(imodel, qs, Es, dqgrid, dEgrid, figsize=(10,10))
                plot.add_plot(width=pylatex.NoEscape(width))
                plot.add_caption('Interpolated resolution functions for %s' % sl.name)
                plt.close()
            doc.append(pylatex.utils.NoEscape(r"\clearpage"))
            
        # one by one comparison plots
        with doc.create(pylatex.Section('Comparing fits to mcvine simulations')):
            for qE, fitter in qE2fitter.items():
                with doc.create(pylatex.Figure(position='htbp')) as plot:
                    plt.figure()
                    plot_compare_fit_to_data(fitter)
                    plot.add_plot(width=pylatex.NoEscape(width))
                    plot.add_caption('Resolution at q=%s, E=%s' % qE)
                    plt.close()
                doc.append(pylatex.utils.NoEscape(r"\clearpage")) # otherwise latex complain about "too many floats"
                
        # save PDF
        doc.generate_pdf(clean_tex=False)
        continue
    return

def create_convolution_calculator(slice, resolution_range=None):
    """Create a "convoler", instance of .convolve2d.Convolver from slice convolution specs
    
    Example slice convolution specification:

        %%file convolution_config.py
        import os, numpy as np, imp
        narr = np.array
        from dgsres import axis
        here = os.path.abspath(os.path.dirname(__file__) or '.')
        # load the resolution calculation configration
        rwc = imp.load_source('rwc', os.path.join(here, './mcvinesim/resolution_workflow_config.py'))
        resolution_datadir = os.path.join(here, 'mcvinesim')

        class Slice_00L(rwc.Slice_00L):

            class expdata:
                "The experimental data to model"
                class grid:
                    qaxis = axis(min=0, max=4.+1e-5, step=0.5/20)
                    Eaxis = axis(min=10., max=27.+1e-5, step=0.2)
                perp_hkl_directions = narr([[1.,0.,0.], [0.,1.,0.]])
                dh = 0.1
                perp_hkl_range = narr([[-dh, dh], [-dh, dh]])

            class convolution:
                expansion_ratio = 0.1
                N_subpixels = 5,5
    """
    from dgsres import axis
    from . import convolve2d as cvv2

    if not _qEgrid_bigger_than(slice.grid, slice.expdata.grid):
        sys.stderr.write("slice %s: Resolution calculation grid smaller than exp data grid" % slice.name)

    grid = cvv2.Grid(slice.expdata.grid.qaxis, slice.expdata.grid.Eaxis)
    expansion_ratio = slice.convolution.expansion_ratio
    qticks = slice.expdata.grid.qaxis.ticks()
    hkl_start = slice.hkl0 + slice.hkl_projection*qticks[0]
    hkl_end =  slice.hkl0 + slice.hkl_projection*qticks[-1]
    perp_hkl_directions = slice.expdata.perp_hkl_directions
    perp_hkl_range = slice.expdata.perp_hkl_range
    def res_func(q, E):
        import os
        pwd = os.path.abspath(os.curdir)
        os.chdir(slice.resolution_datadir)
        imodel = get_interped_resolution_model(slice)
        ret = imodel.getModel(q,E)
        os.chdir(pwd)
        return ret
    N_subpixels = slice.convolution.N_subpixels
    if resolution_range is None:
        res_grid = slice.res_2d_grid
        dqticks = res_grid.qaxis.ticks()
        dEticks = res_grid.Eaxis.ticks()
        resolution_range = dqticks[-1]-dqticks[0], dEticks[-1]-dEticks[0]
    convolver = cvv2.Convolver(grid, hkl_start, hkl_end, expansion_ratio, N_subpixels, res_func, resolution_range, transpose_res_matrix=False)
    slice.convolution.calculator = convolver
    if not _qEgrid_bigger_than(slice.grid, slice.convolution.calculator.finer_expanded_grid):
        sys.stderr.write("slice %s: Resolution calculation grid smaller than convolution expanded grid" % slice.name)
    return slice

def _qEgrid_bigger_than(grid1, grid2):
    qticks1 = grid1.qaxis.ticks()
    qticks2 = (getattr(grid2, 'qaxis', None) or grid2.xaxis).ticks()
    Eticks1 = grid1.Eaxis.ticks()
    Eticks2 = (getattr(grid2, 'Eaxis', None) or grid2.yaxis).ticks()
    return qticks1[0]<qticks2[0] and qticks1[-1]>qticks2[-1] \
        and Eticks1[0]<Eticks2[0] and Eticks1[-1]>Eticks2[-1]


def get_interped_resolution_model(sl):
    import pickle as pkl
    qE2fitres = pkl.load(open('%s-fit_results.pkl' % sl.name))
    qE2fitres = fill_in_blanks_for_fits(qE2fitres, sl)
    return create_interp_model(qE2fitres, sl)

def fill_in_blanks_for_fits(qE2fitter, sl):
    """due to limit of dynamical range measured, some grid points no data is available. 
    has to fill in blanks at corners"""
    # find corners
    qticks = sl.grid.qaxis.ticks()
    qmin = qticks[0]; qmax = qticks[-1]
    Eticks = sl.grid.Eaxis.ticks()
    Emin = Eticks[0]; Emax = Eticks[-1]
    corners = [(qmin, Emin), (qmax, Emin), (qmin, Emax), (qmax, Emax)]
    #
    for corner in corners:
        if corner in qE2fitter: continue
        q1, E1 = corner
        Esearch = Eticks
        if E1==Eticks[-1]: Esearch = Esearch[::-1]
        qsearch = qticks
        if q1==qticks[-1]: qsearch = qsearch[::-1]
        found = False
        for E2 in Esearch:
            for q2 in qsearch:
                if (q2,E2) in qE2fitter:
                    found = True
                    break
            if found: break
        if found:
            qE2fitter[corner] = qE2fitter[(q2, E2)]
    return qE2fitter
        
                    
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


def simulate_all_grid_points(slice, config, Nrounds_beam=1, overwrite=False):
    failed = []; outputs = {}
    for q in slice.grid.qaxis.ticks():
        for E in slice.grid.Eaxis.ticks():
            simdir = config.simdir(q,E, slice)
            if not overwrite and os.path.exists(simdir): continue
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
            E = Es[nrows-irow-1]
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


def fit(q, E, slice, config, use_cache=False, extra_fitting_params=None):
    if use_cache:
        import dill
        path = '%s-q_%.3f-E_%.3f-fitter.dill' % (slice.name, q, E)
        if os.path.exists(path):
            return dill.load(open(path))
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
    if extra_fitting_params: fitting_params.update(extra_fitting_params)
    fitter.fit_result = fitter.fit(**fitting_params)
    if use_cache:
        dill.dump(fitter, open(path, 'w'))
    return fitter


def plot_compare_fit_to_data(fitter, figsize=(8,8)):
    res_z = fitter.res_z
    qgrid, Egrid = fitter.qEgrids
    result = fitter.fit_result
    plt.figure(figsize=figsize)
    plt.subplot(2,2,1)
    plt.title('mcvine sim')
    plt.pcolormesh(qgrid, Egrid, res_z)
    plt.colorbar()
    plt.subplot(2,2,2)
    plt.title('fit')
    scale = res_z.sum()/result.best_fit.sum()
    iqe_fit = result.best_fit.reshape(res_z.shape)*scale
    plt.pcolormesh(qgrid, Egrid, iqe_fit)
    plt.colorbar()

    qs = qgrid[0]; Es = Egrid[:,0]
    plt.subplot(2,2,3)
    plt.plot(qs, res_z.sum(0), label='mcvine sim')
    plt.plot(qs, iqe_fit.sum(0), label='fit')
    plt.legend()
    plt.subplot(2,2,4)
    plt.plot(Es, res_z.sum(1), label='mcvine sim')
    plt.plot(Es, iqe_fit.sum(1), label='fit')
    plt.legend()
    
    return

def fit_all_grid_points(slice, config, use_cache=False):
    qE2fitter = dict()
    nofit = []
    for q in slice.grid.qaxis.ticks():
        for E in slice.grid.Eaxis.ticks():
            print q, E
            try:
                qE2fitter[(q,E)] = fit(q, E, slice, config, use_cache=use_cache)
            except:
                # import traceback as tb
                # tb.print_exc()
                nofit.append((q,E))
        continue
    # correct alpha
    # in some cases, candidate values of alpha can be off by pi/2 but give very similar
    # agreements between the model and the data.
    # therefore, we need to go through the alpha values of all data and find the outliers,
    # and change them by 90 degrees
    outliers = []
    # create an "image" of alpha values
    qticks = slice.grid.qaxis.ticks()
    Eticks = slice.grid.Eaxis.ticks()
    alpha_image = np.ones( (len(qticks), len(Eticks)) ) * np.nan
    for iq, q1 in enumerate(qticks):
        for iE, E1 in enumerate(Eticks):
            this = q1, E1
            if this not in qE2fitter: continue
            alpha_image[iq, iE] = qE2fitter[this].fit_result.best_values['alpha']
        continue
    # find outliers
    N_alphas = (alpha_image==alpha_image).sum()
    outliers, median = _find_outliers(alpha_image, .26*np.pi, Niter=5)
    # for each outlier, change alpha value to something similar to neighbor values
    for iq, q1 in enumerate(qticks):
        for iE, E1 in enumerate(Eticks):
            if not outliers[iq, iE]: continue
            this = q1, E1
            print "* working on outlier %s" % (this,)
            old_alpha = alpha_image[iq, iE]
            m = median[iq, iE]
            alpha_bounds = (m-np.pi/10, m+np.pi/10)
            fitter = fit(
                q1, E1, slice, config, use_cache=False,
                extra_fitting_params=dict(return_all_results=True)
            )
            # choose the best one within the range
            results = fitter.fit_result
            within_bounds  = lambda a: a<alpha_bounds[1] and a>alpha_bounds[0]
            results = [r for r in results if within_bounds(r.best_values['alpha'])]
            fitter.fit_result = results[0]
            qE2fitter[this] = fitter
            print "   old alpha: %s. median alpha: %s. new alpha: %s" % (
                old_alpha, m, fitter.fit_result.best_values['alpha'])
        continue
    #
    # qEranges is an import parameter that need to be rememberd
    fitter1 = qE2fitter.values()[0]
    slice.res_2d_grid.qEranges = fit_ellipsoid.qEgrid2range(*fitter1.qEgrids)
    for fitter in qE2fitter.values()[1:]:
        assert np.allclose(slice.res_2d_grid.qEranges, fit_ellipsoid.qEgrid2range(*fitter.qEgrids)), \
            "qEranges mistmatch: %s vs %s" % (slice.res_2d_grid.qEranges, fit_ellipsoid.qEgrid2range(*fitter.qEgrids))
    return qE2fitter, nofit

def _find_outliers(img, max_deviation, Niter):
    img2 = img.copy()
    for i in range(Niter):
        outliers, median = _find_outliers_1(img2, max_deviation)
        if outliers.sum()==0: break
        img2[outliers] = median[outliers]
    return (img != img2)&(img==img), median

def _find_outliers_1(img, max_deviation):
    """find outliers in img
    
    an outlier is defined as with value significantly different from the median value of its neighborhood
    `>max_deviation` means too different.

    returns outlier map and median map
    """
    # change image value range to 0-240. scikit image median does not work for random floats
    scale = 240; amax = np.nanmax(img); amin = np.nanmin(img)
    img1 = img-amin; img1*=scale/(amax-amin) 
    from skimage.filters import median
    alpha_median = median(img1.astype('uint8'), mask=(img==img))
    return np.abs(img1-alpha_median)>max_deviation*scale/(amax-amin), alpha_median*(amax-amin)/scale+amin

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
            E = Es[nrows-irow-1]
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
        fitter = qE2fitter[qe]
        fr = fitter.fit_result
        fr_tosave = fit_ellipsoid.FitResult()
        fr_tosave.best_values = fr.best_values
        fr_tosave.qEranges = fit_ellipsoid.qEgrid2range(*fitter.qEgrids)
        qE2fitres_tosave[qe] = fr_tosave
        continue

    import pickle as pkl
    pkl.dump(qE2fitres_tosave, open(path, 'w'))
    return

def format_parameter_table(qE2fitres):
    keys = qE2fitres.values()[0].best_values.keys()
    lines = []
    line = "%6s%6s" % ('q','E')
    for k in keys: line += '%8s' % k[:8]
    lines.append(line)
    qEs = list(qE2fitres.keys())
    qEs.sort(key=lambda x: (x[1],x[0]))
    for q, E in qEs:
        if not (q,E) in qE2fitres: continue
        result = qE2fitres[(q,E)]
        line = "%6.1f%6.1f" % (q,E)
        for k in keys:
            v = result.best_values[k]
            line += '%8.4f' % v
        lines.append(line)
    return '\n'.join(lines)

def print_parameter_table(qE2fitres):
    s = format_parameter_table(qE2fitres)
    return s

def create_interp_model(qE2fitres, slice):
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
    qrange, Erange = result.qEranges
    return fit_ellipsoid.InterpModel(qE_points, param_values, qrange, Erange)

def plot_interpolated_resolution_on_grid(model, qs, Es, dqgrid, dEgrid, figsize=(10,10)):
    ncols = len(qs)
    nrows = len(Es)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    for irow in range(nrows):
        for icol in range(ncols):
            q = qs[icol]
            E = Es[nrows-irow-1]
            ax1 = axes[irow][icol]
            ax1.set_title("q=%.2f, E=%.2f" % (q, E))
            z = model.getModel(q=q, E=E).ongrid(dqgrid, dEgrid)
            ax1.pcolormesh(dqgrid, dEgrid, z)
            continue
        continue
    plt.tight_layout()
    return
