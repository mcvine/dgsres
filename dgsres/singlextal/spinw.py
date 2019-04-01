# use spinw to calculate S(q, w)
#
# The methods here require a data structure for a slice
# The code here was derived from exploratory code in https://jupyter.sns.gov/user/{UID}/notebooks/data/SNS/SEQ/IPTS-21411/shared/resolution/spinw.ipynb#
#
# The methods here require the "convolver" is already created by calling workflow.create_convolution_calculator(slice)


import os, numpy as np, tqdm
import matlab.engine, matlab
from . import disp2sqw as d2s


def get_dispersions_along_slice(ml_slice_func, slice, Nq_disp=500, **kwds):
    return d2s.get_dispersions_along_slice(wrap_spinw_disp_func(ml_slice_func), slice, Nq=Nq_disp, **kwds)


def get_thin_slice(ml_slice_func, slice, **kwds):
    return d2s.get_thin_slice(wrap_spinw_disp_func(ml_slice_func), slice, **kwds)


def get_slice(ml_slice_func, slice, Nq_disp=500, Nsample_perp=None, sampling_method='linear', **kwds):
    return d2s.get_slice(
        wrap_spinw_disp_func(ml_slice_func), slice,
        Nq_sample=Nq_disp, Nsample_perp=Nsample_perp, sampling_method=sampling_method,
        **kwds)


def _backward_compatible_func(f, name):
    def _(*args, **kwds):
        import warnings
        warnings.warn("`%s` is obsolete. Use `%s`" % (name, f.__name__))
        return f(*args, **kwds)
    return _

get_dispersions_along_slice_using_spinw = _backward_compatible_func(get_dispersions_along_slice, 'get_dispersions_along_slice_using_spinw')
get_thin_slice_using_spinw = _backward_compatible_func(get_thin_slice, 'get_thin_slice_using_spinw')
get_slice_using_spinw = _backward_compatible_func(get_slice, 'get_slice_using_spinw')


def wrap_spinw_disp_func(f):
    def _(start, end, Nq_disp):
        sp = f(matlab.double(list(start)), matlab.double(list(end)), Nq_disp)
        return sp['omega'], sp['swInt']
    return _
