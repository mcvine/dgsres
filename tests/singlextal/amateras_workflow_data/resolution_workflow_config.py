import os
thisdir = os.path.abspath(os.path.dirname(__file__) or '.')

import numpy as np

# instrument
from dgsres.instruments import amateras as instrument   # CHANGE INSTRUMENT HERE

# sample
sample_yaml = os.path.join(thisdir, 'sample.yaml')

# exp condition
beam = os.path.join(thisdir, "beam/")
Ei = 2.63

from mcvine.workflow import singlextal as sx
psi_scan = sx.axis(min=-10., max=181., step=1.)  # Angular Range

# sim directory name
def simdir(q, E, slice):
    return 'sim-%s-q_%.3f,E_%.3f' % (slice.name, q, E)
sim_Nrounds_beam = 10

# slice
#
class Slice_H00:
    name = 'H00'
    hkl_projection = np.array([1.,0,0])
    hkl0 = np.array([0.,0,0])

    class grid:
        "simulations will be done for points on this grid"
        # qaxis = sx.axis(min=-2.0, max=2.51, step=0.5)
        # Eaxis = sx.axis(min=-0.5, max=2.51, step=0.25)
        qaxis = sx.axis(min=-2.0, max=3.01, step=1.0)
        Eaxis = sx.axis(min=-0.5, max=2.51, step=0.5)
GammaA = Slice_H00

#
class Slice_H10:
    name = 'H10'
    hkl_projection = np.array([1.,0,0])
    hkl0 = np.array([0,1.,0])

    class grid:
        "simulations will be done for points on this grid"
        # qaxis = sx.axis(min=-2.0, max=2.51, step=0.5)
        # Eaxis = sx.axis(min=-0.5, max=2.51, step=0.25)
        qaxis = sx.axis(min=-2.0, max=2.51, step=1.5)
        Eaxis = sx.axis(min=-0.5, max=2.51, step=0.7)
        # qaxis = sx.axis(min=-2.0, max=2.51, step=2.0)
        # Eaxis = sx.axis(min=-0.5, max=2.51, step=1.0)
# 
class Slice_H20:
    name = 'H20'
    hkl_projection = np.array([1.,0,0])
    hkl0 = np.array([0,2.,0])

    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-2.0, max=2.51, step=0.5)
        Eaxis = sx.axis(min=-0.5, max=2.51, step=0.25)

# 
class Slice_0K0:
    name = '0K0'
    hkl_projection = np.array([0,1.,0])
    hkl0 = np.array([0,0,0])

    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-4.0, max=4.01, step=0.5)
        Eaxis = sx.axis(min=-0.5, max=2.51, step=0.25)

# 
class Slice_1K0:
    name = '1K0'
    hkl_projection = np.array([0,1.,0])
    hkl0 = np.array([1.,0,0])

    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-3.5, max=3.51, step=0.5)
        Eaxis = sx.axis(min=-0.5, max=2.51, step=0.25)

# 
class Slice_0p5K0:
    name = '0p5K0'
    hkl_projection = np.array([0,1.,0])
    hkl0 = np.array([0.5,0,0])

    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-4.0, max=3.51, step=0.5)
        Eaxis = sx.axis(min=-0.5, max=2.51, step=0.25)
        
#
slices = [Slice_H00, Slice_H10, Slice_H20, Slice_0K0, Slice_1K0, Slice_0p5K0]
slices = [Slice_H10]

class res_2d_grid:
    "resolution data will be histogrammed into this grid"
    qaxis = sx.axis(min=-0.125, max=0.126, step=0.005)
    Eaxis = sx.axis(min=-0.2, max=0.201, step=0.008)

class fitting:
    rounds = 3
    gaussian2d_threshold = 0.5
    alpha_bounds = (-np.pi/2, np.pi/2)

for sl in slices:
    sl.res_2d_grid = res_2d_grid
    sl.fitting = fitting

