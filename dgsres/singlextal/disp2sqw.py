# calculate S(q, w) from dispersion
#
# The methods here require a data structure for a slice
# The code here was derived from exploratory code in https://jupyter.sns.gov/user/{UID}/notebooks/data/SNS/SEQ/IPTS-21411/shared/resolution/spinw.ipynb#
#
# The methods here require the "convolver" is already created by calling workflow.create_convolution_calculator(slice)

import os, numpy as np, tqdm
# import matlab.engine, matlab

def get_dispersions_along_slice(disp_calc, slice1, Nq=500, branches=slice(None, None)):
    """obtain dispersion at the given slice between convolver expanded grid
    
    disp_calc: dispersion calculation function. returns (omega, intensity).
        omega: array. shape: branch,Nq
        intensity: array. shape: branch,Nq
    
    slice1: slice spec. examples see workflow.py

    Nq: number of q points
    """
    convolver = slice1.convolution.calculator
    # call dispersion calculator
    hkl_start, hkl_end = convolver.expanded_hkl_start, convolver.expanded_hkl_end
    omega, intensity = disp_calc(hkl_start, hkl_end, Nq)
    qmin = get_q(hkl_start, slice1); qmax = get_q(hkl_end, slice1)
    qs = np.linspace(qmin, qmax, Nq)
    return qs, omega[branches, :], intensity[branches, :]

def get_q(hkl, slice):
    return np.dot( hkl-slice.hkl0, slice.hkl_projection ) / np.dot(slice.hkl_projection, slice.hkl_projection)

def get_thin_slice(disp_calc, slice, branches=None):
    convolver = slice.convolution.calculator
    qs = convolver.finer_expanded_grid.xaxis.ticks()
    Es = convolver.finer_expanded_grid.yaxis.ticks()
    Nq = len(qs)
    NE = len(Es)
    slice_img = np.zeros((Nq,NE))

    _qs, omega, intensity = get_dispersions_along_slice(disp_calc, slice)
    nbr = omega.shape[0]

    Emin = convolver.finer_expanded_grid.yaxis.min
    Estep = convolver.finer_expanded_grid.yaxis.step

    for branch in branches or range(nbr):
        disp1_Es = np.interp(qs, _qs, omega[branch])
        disp1_Is = np.interp(qs, _qs, intensity[branch])
        disp1_Ebinindexes = np.array((disp1_Es - Emin)/Estep, dtype=int)
        qindexes = np.arange(Nq, dtype=int)
        good = np.logical_and(disp1_Ebinindexes>=0, disp1_Ebinindexes<NE)
        slice_img[qindexes[good], disp1_Ebinindexes[good]] += disp1_Is[qindexes[good]]
    qmg,Emg = np.meshgrid(qs, Es)
    return qmg, Emg, slice_img


def get_slice(disp_calc, slice, Nq_sample=None, Nsample_perp=None, sampling_method='linear', branches=None):
    """get slice integrated over thickness (defined in slice object. see workflow.py for example)
    """
    convolver = slice.convolution.calculator
    qs = convolver.finer_expanded_grid.xaxis.ticks()
    Es = convolver.finer_expanded_grid.yaxis.ticks()
    Nq = len(qs)
    NE = len(Es)
    slice_img = np.zeros((Nq,NE))
    # q values sampled (calculated by disp_calc, to be interpolated later)
    hkl_start, hkl_end = convolver.expanded_hkl_start, convolver.expanded_hkl_end
    qmin = get_q(hkl_start, slice); qmax = get_q(hkl_end, slice)
    Nq_sample = Nq_sample or (Nq*2)
    _qs = np.linspace(qmin, qmax, Nq_sample)
    # energy
    Emin = convolver.finer_expanded_grid.yaxis.min
    Estep = convolver.finer_expanded_grid.yaxis.step
    # 
    perp1_range, perp2_range = slice.expdata.perp_hkl_range
    perp1_dir, perp2_dir = slice.expdata.perp_hkl_directions
    Nsample_perp = Nsample_perp or slice.expdata.Nsample_perp_hkl
    if sampling_method == 'linear':
        perp_qs = [(qp1, qp2)
                   for qp1 in np.linspace(perp1_range[0], perp1_range[1], Nsample_perp)
                   for qp2 in np.linspace(perp2_range[0], perp2_range[1], Nsample_perp)
        ]
    elif sampling_method == 'uniform':
        rs = np.random.random( (Nsample_perp*Nsample_perp, 2) )
        rs[:, 0] *= perp1_range[1]-perp1_range[0]; rs[:, 0] += perp1_range[0]
        rs[:, 1] *= perp2_range[1]-perp2_range[0]; rs[:, 0] += perp2_range[0]
        perp_qs = rs
    else:
        raise ValueError("Unknown sampling method %r" % sampling_method)
    #
    for q_perp1, q_perp2 in tqdm.tqdm(perp_qs):
        # the q offset at the perpendicular directions
        q_offset = q_perp1*perp1_dir+ q_perp2*perp2_dir
        # print q_offset
        # start, end of q points with offset
        q_start = convolver.expanded_hkl_start + q_offset
        q_end = convolver.expanded_hkl_end + q_offset
        # call spinw
        omega, intensity = disp_calc(q_start, q_end, Nq_sample)
        nbr = omega.shape[0]
        # loop over branches
        for branch in branches or range(nbr):
            disp1_Es = np.interp(qs, _qs, omega[branch])
            disp1_Is = np.interp(qs, _qs, intensity[branch])
            disp1_Ebinindexes = np.array((disp1_Es - Emin)/Estep, dtype=int)
            qindexes = np.arange(Nq, dtype=int)
            good = np.logical_and(disp1_Ebinindexes>=0, disp1_Ebinindexes<NE)
            # print qindexes[good]
            # slice_img[qindexes[good], disp1_Ebinindexes[good]] += 1
            slice_img[qindexes[good], disp1_Ebinindexes[good]] += disp1_Is[qindexes[good]]
    qmg,Emg = np.meshgrid(qs, Es)
    return qmg, Emg, slice_img


def get_slice_fast(disp_calc, slice, Nq_sample=None, Nsample_perp=None, sampling_method='linear', branches=None):
    """get slice integrated over thickness (defined in slice object. see workflow.py for example)
    """
    convolver = slice.convolution.calculator
    qs = convolver.finer_expanded_grid.xaxis.ticks()
    Es = convolver.finer_expanded_grid.yaxis.ticks()
    Nq = len(qs)
    NE = len(Es)
    slice_img = np.zeros((Nq,NE))
    # q values sampled (calculated by disp_calc, to be interpolated later)
    hkl_start, hkl_end = convolver.expanded_hkl_start, convolver.expanded_hkl_end
    qmin = get_q(hkl_start, slice); qmax = get_q(hkl_end, slice)
    # energy
    Emin = convolver.finer_expanded_grid.yaxis.min
    Estep = convolver.finer_expanded_grid.yaxis.step
    # 
    perp1_range, perp2_range = slice.expdata.perp_hkl_range
    perp1_dir, perp2_dir = slice.expdata.perp_hkl_directions
    Nsample_perp = Nsample_perp or slice.expdata.Nsample_perp_hkl
    if sampling_method == 'linear':
        perp_qs = [(qp1, qp2)
                   for qp1 in np.linspace(perp1_range[0], perp1_range[1], Nsample_perp)
                   for qp2 in np.linspace(perp2_range[0], perp2_range[1], Nsample_perp)
        ]
    elif sampling_method == 'uniform':
        rs = np.random.random( (Nsample_perp*Nsample_perp, 2) )
        rs[:, 0] *= perp1_range[1]-perp1_range[0]; rs[:, 0] += perp1_range[0]
        rs[:, 1] *= perp2_range[1]-perp2_range[0]; rs[:, 0] += perp2_range[0]
        perp_qs = rs
    else:
        raise ValueError("Unknown sampling method %r" % sampling_method)
    #
    # print "build hkls"
    hkls = [] # 
    for q_perp1, q_perp2 in perp_qs:
        # the q offset at the perpendicular directions
        q_offset = q_perp1*perp1_dir+ q_perp2*perp2_dir
        # print q_offset
        # start, end of q points with offset
        q_start = convolver.expanded_hkl_start + q_offset
        q_end = convolver.expanded_hkl_end + q_offset
        _x = np.linspace(0., 1., Nq)
        hkls1 = q_start +(q_end-q_start)*_x[:, np.newaxis] # Nq, 3
        hkls.append( hkls1 )
        continue
    # print "built hkls"
    hkls = np.array(hkls) # nq_perp1and2, Nq, 3
    hkls.shape = -1, 3  # nq_perp1and2 X Nq, 3
    # call spinw
    # print "call spinw"
    omega, intensity = disp_calc(hkls) # both shapes: nbranch, nhkls
    # print "spinw done"
    #
    nbr = omega.shape[0]
    omega.shape = nbr, -1, Nq # nbr, nq_perp1and2, Nq
    nq_perp1and2 = omega.shape[1]
    intensity.shape = nbr, nq_perp1and2, Nq
    # loop over branches
    qindexes = np.arange(Nq, dtype=int)
    for branch in branches or range(nbr):
        for iq_perp in range(nq_perp1and2):
            disp1_Es = omega[branch, iq_perp]
            disp1_Is = intensity[branch, iq_perp]
            disp1_Ebinindexes = np.array((disp1_Es - Emin)/Estep, dtype=int)
            good = np.logical_and(disp1_Ebinindexes>=0, disp1_Ebinindexes<NE)
            slice_img[qindexes[good], disp1_Ebinindexes[good]] += disp1_Is[qindexes[good]]
    qmg,Emg = np.meshgrid(qs, Es)
    return qmg, Emg, slice_img


"""
Does not work. also not faster

    # make "events"
    omega.shape = intensity.shape =  -1,
    #
    nrepeats = omega.size// Nq
    qevents = np.tile(qs, nrepeats)
    # boundaries
    dq = qs[1]-qs[0]; qbb = np.append(qs-dq/2, qs[-1]+dq/2)
    dE = Es[1]-Es[0]; Ebb = np.append(Es-dE/2, Es[-1]+dE/2)
    img, xedges, yedges = np.histogram2d(omega, qevents, bins=(qbb, Ebb), weights=intensity)
    qmg,Emg = np.meshgrid(qs, Es)
    return qmg, Emg, slice_img
    
"""
