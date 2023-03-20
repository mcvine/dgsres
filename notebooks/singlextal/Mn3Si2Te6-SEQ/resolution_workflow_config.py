import os
thisdir = os.path.abspath(os.path.dirname(__file__) or '.')

import numpy as np

# instrument
from dgsres.instruments import sequoia as instrument

# sample
sample_yaml = os.path.join(thisdir, 'sample.yaml')

# exp condition
beam = os.path.join(thisdir, "beam-60meV-n5e9/")
Ei = 60.47848057039433

from mcvine.workflow import singlextal as sx
psi_scan = sx.axis(min=0., max=180., step=1.)

# sim directory name
def simdir(q, E, slice):
    return 'sim-%s-q_%.3f,E_%.3f' % (slice.name, q, E)
sim_Nrounds_beam = 30

# slice
# 
class Slice_00L:
    name = '00L'
    hkl_projection = np.array([0,0,1.])
    hkl0 = np.array([0,0,0.])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-.5, max=4.5, step=0.6)
        Eaxis = sx.axis(min=5., max=31., step=5.)
GammaA = Slice_00L

# 
class Slice_H00:
    name = 'H00'
    hkl_projection = np.array([1.,0,0.])
    hkl0 = np.array([0.,0,0])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-0.5, max=3., step=0.4)
        Eaxis = sx.axis(min=5., max=31., step=5.)

#
class Slice_H01:
    name = 'H01'
    hkl_projection = np.array([1.,0,0.])
    hkl0 = np.array([0,0,1.])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-0.5, max=3., step=0.2)
        Eaxis = sx.axis(min=5., max=31., step=5.)

#
class Slice_H02:
    name = 'H02'
    hkl_projection = np.array([1.,0,0.])
    hkl0 = np.array([0,0,2.])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-0.5, max=3., step=0.4)
        Eaxis = sx.axis(min=5., max=31., step=5.)
        
# 
class Slice_H03:
    name = 'H03'
    hkl_projection = np.array([1.,0,0])
    hkl0 = np.array([0,0,3.])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-0.3, max=2.6, step=0.4)
        Eaxis = sx.axis(min=5., max=31., step=5.)

# 
class Slice_H0H:
    name = 'H0H'
    hkl_projection = np.array([1,0,1.])
    hkl0 = np.array([0,0,0.])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-0.3, max=3.35, step=0.6)
        Eaxis = sx.axis(min=5., max=31., step=5.)

# 
class Slice_0K1:
    name = '0K1'
    hkl_projection = np.array([0,1.,0.])
    hkl0 = np.array([0,0.,1.])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-0.5, max=3., step=0.4)
        Eaxis = sx.axis(min=5., max=31., step=5.)
        
# 
class Slice_mK2K2:
    name = 'mK2K2'
    hkl_projection = np.array([-1.,2.,0.])
    hkl0 = np.array([0,0.,2.])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-0.5, max=3., step=0.4)
        Eaxis = sx.axis(min=5., max=31., step=5.)
        
# 
class Slice_mK2K3:
    name = 'mK2K3'
    hkl_projection = np.array([-1.,2.,0.])
    hkl0 = np.array([0,0.,3.])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-0.5, max=3., step=0.4)
        Eaxis = sx.axis(min=5., max=31., step=5.)     
           

#
#slices = [Slice_00L, Slice_H00, Slice_H01, Slice_H02, Slice_H03, Slice_H0H, Slice_0K1, Slice_mK2K2, Slice_mK2K3]
slices = [Slice_00L]

class res_2d_grid:
    "resolution data will be histogrammed into this grid"
    qaxis = sx.axis(min=-0.5, max=0.5, step=0.02)
    Eaxis = sx.axis(min=-5., max=5., step=.2)

class fitting:
    rounds = 3
    gaussian2d_threshold = 0.5
    alpha_bounds = (-np.pi/2, np.pi/2)

for sl in slices:
    sl.res_2d_grid = res_2d_grid
    sl.fitting = fitting
