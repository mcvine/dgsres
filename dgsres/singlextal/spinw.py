# use spinw to calculate S(q, w)
#
# The methods here require a data structure for a slice
# The code here was derived from exploratory code in https://jupyter.sns.gov/user/{UID}/notebooks/data/SNS/SEQ/IPTS-21411/shared/resolution/spinw.ipynb#
#
# The methods here require the "convolver" is already created by calling workflow.create_convolution_calculator(slice)

import os, numpy as np, tqdm
import matlab.engine, matlab

def get_dispersions_along_slice_using_spinw(ml_slice_func, slice, Nq_disp=500):
    convolver = slice.convolution.calculator
    spec = ml_slice_func(
        matlab.double(list(convolver.expanded_hkl_start)),
        matlab.double(list(convolver.expanded_hkl_end)), Nq_disp)
    
    qs = convolver.finer_expanded_grid.xaxis.ticks()
    omega = np.asarray(spec['omega'])
    swInt = np.asarray(spec['swInt'])
    # q values for the dispersion data calcd by spinw
    disp_qstep = (qs[-1]-qs[0])/(Nq_disp-1)
    spinw_qs = np.arange(qs[0], qs[-1]+disp_qstep/2., disp_qstep)
    return spinw_qs, omega, swInt

def get_thin_slice_using_spinw(ml_slice_func, slice):
    convolver = slice.convolution.calculator
    qs = convolver.finer_expanded_grid.xaxis.ticks()
    Es = convolver.finer_expanded_grid.yaxis.ticks()
    Nq = len(qs)
    NE = len(Es)
    slice_img = np.zeros((Nq,NE))

    spinw_qs, omega, swInt = get_dispersions_along_slice_using_spinw(ml_slice_func, slice)

    Emin = convolver.finer_expanded_grid.yaxis.min
    Estep = convolver.finer_expanded_grid.yaxis.step

    for branch in range(omega.shape[0]//2):
        disp1_Es = np.interp(qs, spinw_qs, omega[branch])
        disp1_Is = np.interp(qs, spinw_qs, swInt[branch])
        disp1_Ebinindexes = np.array((disp1_Es - Emin)/Estep, dtype=int)
        qindexes = np.arange(Nq, dtype=int)
        good = np.logical_and(disp1_Ebinindexes>=0, disp1_Ebinindexes<NE)
        slice_img[qindexes[good], disp1_Ebinindexes[good]] += disp1_Is[qindexes[good]]
    qmg,Emg = np.meshgrid(qs, Es)
    return qmg, Emg, slice_img


def get_slice_using_spinw(ml_slice_func, slice, Nq_disp=500, Nsample_perp=17):
    convolver = slice.convolution.calculator
    qs = convolver.finer_expanded_grid.xaxis.ticks()
    Es = convolver.finer_expanded_grid.yaxis.ticks()
    Nq = len(qs)
    NE = len(Es)
    slice_img = np.zeros((Nq,NE))

    # q values for the dispersion data calcd by spinw
    disp_qstep = (qs[-1]-qs[0])/(Nq_disp-1)
    spinw_qs = np.arange(qs[0], qs[-1]+disp_qstep/2., disp_qstep)

    Emin = convolver.finer_expanded_grid.yaxis.min
    Estep = convolver.finer_expanded_grid.yaxis.step

    perp1_range, perp2_range = slice.expdata.perp_hkl_range
    perp1_dir, perp2_dir = slice.expdata.perp_hkl_directions
    for q_perp1 in tqdm.tqdm(np.linspace(perp1_range[0], perp1_range[1], Nsample_perp)):
        for q_perp2 in np.linspace(perp2_range[0], perp2_range[1], Nsample_perp):
            # the q offset at the perpendicular directions
            q_offset = q_perp1*perp1_dir+ q_perp2*perp2_dir
            # print q_offset
            # start, end of q points with offset
            q_start = convolver.expanded_hkl_start + q_offset
            q_end = convolver.expanded_hkl_end + q_offset
            # call spinw
            spec = ml_slice_func(matlab.double(list(q_start)), matlab.double(list(q_end)), Nq_disp)
            # get data
            omega = np.asarray(spec['omega'])
            swInt = np.asarray(spec['swInt'])
            # loop over branches
            for branch in range(omega.shape[0]//2):
                disp1_Es = np.interp(qs, spinw_qs, omega[branch])
                disp1_Is = np.interp(qs, spinw_qs, swInt[branch])
                disp1_Ebinindexes = np.array((disp1_Es - Emin)/Estep, dtype=int)
                qindexes = np.arange(Nq, dtype=int)
                good = np.logical_and(disp1_Ebinindexes>=0, disp1_Ebinindexes<NE)
                # print qindexes[good]
                # slice_img[qindexes[good], disp1_Ebinindexes[good]] += 1
                slice_img[qindexes[good], disp1_Ebinindexes[good]] += disp1_Is[qindexes[good]]
    qmg,Emg = np.meshgrid(qs, Es)
    return qmg, Emg, slice_img


def _():
    
    grid = cvv2.Grid(Slice_00L.expdata.grid.qaxis, Slice_00L.expdata.grid.Eaxis)
    expansion_ratio = Slice_00L.convolution.expansion_ratio
    qticks = Slice_00L.expdata.grid.qaxis.ticks()
    hkl_start = Slice_00L.hkl0 + Slice_00L.hkl_projection*qticks[0]
    hkl_end =  Slice_00L.hkl0 + Slice_00L.hkl_projection*qticks[-1]
    perp_hkl_directions = Slice_00L.expdata.perp_hkl_directions
    perp_hkl_range = Slice_00L.expdata.perp_hkl_range
    def res_func(q, E):
        import os
        pwd = os.path.abspath(os.curdir)
        os.chdir('../mcvinesim/')
        imodel = get_interped_resolution_model(sl)
        ret = imodel.getModel(q,E)
        os.chdir(pwd)
        return ret
    N_subpixels = Slice_00L.convolution.N_subpixels
    res_range = 0.1, 10
    convolver = cvv2.Convolver(grid, hkl_start, hkl_end, expansion_ratio, N_subpixels, res_func, res_range)
