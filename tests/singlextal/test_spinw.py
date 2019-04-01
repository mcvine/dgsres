#!/usr/bin/env python

import os, numpy as np, unittest as ut
from dgsres.singlextal import spinw, convolve2d as cvv2
from mcvine.workflow import singlextal as sx
plot = False
# plot = True

here = os.path.dirname(__file__)
expected_res_dir = os.path.join(here, 'expected-results')

def ml_slice_func(hkl_start, hkl_end, Nq_disp):
    step =( hkl_end-hkl_start)/(Nq_disp-1)
    Q = hkl_start + np.outer( np.arange(Nq_disp), step)
    Qx, Qy, Qz = Q.T
    E = 20*np.sin(np.pi*(Qx+Qy+Qz))**2 + 10
    I = np.ones(Nq_disp)
    return dict(omega=np.array([E,E]), swInt=np.array([I,I]))


# slice
class Slice_00L:
    name = '00L'
    hkl_projection = np.array([0,0,1.])
    hkl0 = np.array([0,0,0.])

    class grid:
        "simulations will be done for points on this grid"
        qaxis = sx.axis(min=-.5, max=4.5, step=0.6)
        Eaxis = sx.axis(min=-4., max=42.1, step=8.)

    class res_2d_grid:
        "resolution data will be histogrammed into this grid"
        qaxis = sx.axis(min=-0.25, max=0.25, step=0.01)
        Eaxis = sx.axis(min=-2., max=2., step=.05)
        
    class expdata:
        "The experimental data to model"
        class grid:
            qaxis = sx.axis(min=0, max=4.+1e-5, step=0.5/10)
            Eaxis = sx.axis(min=0., max=40.+1e-5, step=1.)
        perp_hkl_directions = np.array([[1.,0.,0.], [0.,1.,0.]])
        dh = 0.1
        perp_hkl_range = np.array([[-dh, dh], [-dh, dh]])
        Nsample_perp_hkl = 25
          
class convolution:
    expansion_ratio = 0.1
    N_subpixels = 5,5

    grid = cvv2.Grid(Slice_00L.expdata.grid.qaxis, Slice_00L.expdata.grid.Eaxis)
    qticks = Slice_00L.expdata.grid.qaxis.ticks()
    hkl_start = Slice_00L.hkl0 + Slice_00L.hkl_projection*qticks[0]
    hkl_end =  Slice_00L.hkl0 + Slice_00L.hkl_projection*qticks[-1]
    res_grid = Slice_00L.res_2d_grid
    dqticks = res_grid.qaxis.ticks()
    dEticks = res_grid.Eaxis.ticks()
    resolution_range = dqticks[-1]-dqticks[0], dEticks[-1]-dEticks[0]
    calculator = cvv2.Convolver(
        grid, hkl_start, hkl_end, 
        expansion_ratio=0.1, N_subpixels=5, res_func=None, res_range=resolution_range, transpose_res_matrix=False
    )
Slice_00L.convolution = convolution


class TestCase(ut.TestCase):
    
    def test_get_dispersions_along_slice(self):
        spinw_qs, omega, swInt \
            = spinw.get_dispersions_along_slice_using_spinw(ml_slice_func, slice=Slice_00L, Nq_disp=500)
        if plot:
            from matplotlib import pyplot as plt
            plt.plot(spinw_qs, omega[0])
            plt.show()
        # uncomment to save the expected results
        # np.save(os.path.join(expected_res_dir, "get_dispersions_along_slice.npy"), (spinw_qs, omega[0], swInt[0]))
        # checking
        expected = np.load(os.path.join(expected_res_dir, "get_dispersions_along_slice.npy"))
        self.assert_(np.allclose(expected, (spinw_qs, omega[0], swInt[0])))
        return
    
    def test_get_thin_slice(self):
        qmg, Emg, slice_img = spinw.get_thin_slice_using_spinw(ml_slice_func, slice=Slice_00L)
        if plot:
            from matplotlib import pyplot as plt
            plt.pcolormesh(qmg, Emg, slice_img.T)
            plt.colorbar()
            plt.show()
        # uncomment to save the expected results
        #     np.save(os.path.join(expected_res_dir, "get_thin_slice.npy"), (qmg, Emg, slice_img.T))
        # checking
        expected = np.load(os.path.join(expected_res_dir, "get_thin_slice.npy"))
        self.assert_(np.allclose(expected, (qmg, Emg, slice_img.T)))
        return

    def test_get_slice(self):
        qmg, Emg, slice_img = spinw.get_slice_using_spinw(ml_slice_func, slice=Slice_00L)
        if plot:
            from matplotlib import pyplot as plt
            plt.pcolormesh(qmg, Emg, slice_img.T)
            plt.colorbar()
            plt.show()
        # uncomment to save the expected results
        #     np.save(os.path.join(expected_res_dir, "get_slice.npy"), (qmg, Emg, slice_img.T))
        # checking
        expected = np.load(os.path.join(expected_res_dir, "get_slice.npy"))
        self.assert_(np.allclose(expected, (qmg, Emg, slice_img.T)))
        return


if __name__=='__main__': ut.main()
