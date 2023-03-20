import os, numpy as np, imp
narr = np.array
from dgsres import axis
from dgsres.singlextal import convolve2d as cvv2
here = os.path.abspath(os.path.dirname(__file__) or '.')
rwc = imp.load_source('rwc', os.path.join(here, 'resolution_workflow_config.py'))

resolution_datadir = '/SNS/SEQ/IPTS-21411/shared/resolution/mcvinesim/'

class Slice_00L(rwc.Slice_00L):
    
    class expdata:
        class grid:
            qaxis = axis(min=0, max=4.+1e-5, step=0.5/10)
            Eaxis = axis(min=10., max=28.+1e-5, step=0.5)
        perp_hkl_directions = narr([[1.,0.,0.], [-1.,2.,0.]])
        dh = 0.1
        perp_hkl_range = narr([[-dh, dh], [-dh, dh]])
        Nsample_perp_hkl = 20  # this is for sampling the perp directions.
        
    class convolution:
        expansion_ratio = 0.1
        N_subpixels = 10,10

class Slice_H00(rwc.Slice_H00):
    
    class expdata:
        class grid:
            qaxis = axis(min=0, max=2.5+1e-5, step=0.5/10)
            Eaxis = axis(min=10., max=28.+1e-5, step=0.5)
        perp_hkl_directions = narr([[-1.,2.,0.], [0.,0.,1.]])
        dh = 0.1
        perp_hkl_range = narr([[-dh, dh], [-dh, dh]])
        Nsample_perp_hkl = 20  # this is for sampling the perp directions.
    
    class convolution:
        expansion_ratio = 0.1
        N_subpixels = 10,10

class Slice_H01(rwc.Slice_H01):
    
    class expdata:
        class grid:
            qaxis = axis(min=0, max=3.+1e-5, step=0.5/10)
            Eaxis = axis(min=10., max=28.+1e-5, step=0.5)
        perp_hkl_directions = narr([[-1.,2.,0.], [0.,0.,1.]])
        dh = 0.1
        perp_hkl_range = narr([[-dh, dh], [-dh, dh]])
        Nsample_perp_hkl = 20  # this is for sampling the perp directions.
    
    class convolution:
        expansion_ratio = 0.1
        N_subpixels = 10,10

class Slice_H02(rwc.Slice_H02):
    
    class expdata:
        class grid:
            qaxis = axis(min=0, max=3.+1e-5, step=0.5/10)
            Eaxis = axis(min=10., max=28.+1e-5, step=0.5)
        perp_hkl_directions = narr([[0.,1.,0.], [0.,0.,1.]])
        dh = 0.1
        perp_hkl_range = narr([[-dh, dh], [-dh, dh]])
        Nsample_perp_hkl = 20  # this is for sampling the perp directions.
    
    class convolution:
        expansion_ratio = 0.1
        N_subpixels = 10,10

class Slice_H03(rwc.Slice_H03):
    
    class expdata:
        class grid:
            qaxis = axis(min=0, max=3.+1e-5, step=0.5/10)
            Eaxis = axis(min=10., max=28.+1e-5, step=0.5)
        perp_hkl_directions = narr([[-1.,2.,0.], [0.,0.,1.]])
        dh = 0.1
        perp_hkl_range = narr([[-dh, dh], [-dh, dh]])
        Nsample_perp_hkl = 20  # this is for sampling the perp directions.
    
    class convolution:
        expansion_ratio = 0.1
        N_subpixels = 10,10

class Slice_H0H(rwc.Slice_H0H):
    
    class expdata:
        class grid:
            qaxis = axis(min=0, max=3.+1e-5, step=0.5/10)
            Eaxis = axis(min=10., max=28.+1e-5, step=0.5)
        perp_hkl_directions = narr([[-0.18225,0.,1.], [-1.,2.,0.]])
        dh = 0.1
        perp_hkl_range = narr([[-dh, dh], [-dh, dh]])
        Nsample_perp_hkl = 20  # this is for sampling the perp directions.
    
    class convolution:
        expansion_ratio = 0.1
        N_subpixels = 10,10

class Slice_0K1(rwc.Slice_0K1):
    
    class expdata:
        class grid:
            qaxis = axis(min=0, max=2.5+1e-5, step=0.5/10)
            Eaxis = axis(min=10., max=28.+1e-5, step=0.5)
        perp_hkl_directions = narr([[1.,0.,0.], [0.,0.,1.]])
        dh = 0.1
        perp_hkl_range = narr([[-dh, dh], [-dh, dh]])
        Nsample_perp_hkl = 20  # this is for sampling the perp directions.
    
    class convolution:
        expansion_ratio = 0.1
        N_subpixels = 10,10

slices = [Slice_00L, Slice_H00, Slice_H01, Slice_H02, Slice_H03, Slice_H0H, Slice_0K1]
# slices = [Slice_00L]
for sl in slices:
    sl.resolution_datadir = resolution_datadir
