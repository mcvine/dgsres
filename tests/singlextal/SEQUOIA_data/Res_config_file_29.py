import os
thisdir = os.path.abspath(os.path.dirname(__file__) or '.')

import numpy as np

# instrument
from dgsres.instruments import sequoia as instrument

# sample
sample_yaml = os.path.join(thisdir, 'sample.yml')

# exp condition
beam = os.path.join(thisdir,"res_files_29","beam/")
Ei = 28.94

from mcvine.workflow import singlextal as sx
psi_scan = sx.axis(min=0., max=180., step=1.)

# sim directory name
def simdir(q, E, slice):
    return 'sim-%s-q_%.3f,E_%.3f' % (slice.name, q, E)
sim_Nrounds_beam = 1

# slice
# 
class Slice_1p5k0:
    name = '1p5k0'
    hkl_projection = np.array([0,1,0])
    hkl0 = np.array([1.5,0,0])
    
    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-2, max=2, step=0.1)
        Eaxis = sx.axis(min=-5, max=29, step=0.5)

slices = [Slice_1p5k0]

class res_2d_grid:
    "resolution data will be histogrammed into this grid"
    qaxis = sx.axis(min=-0.2, max=0.2, step=0.005)
    Eaxis = sx.axis(min=-1, max=1, step=0.05)

class fitting:
    rounds = 3
    gaussian2d_threshold = 0.5
    alpha_bounds = (-np.pi/2, np.pi/2)

for sl in slices:
    sl.res_2d_grid = res_2d_grid
    sl.fitting = fitting
    
