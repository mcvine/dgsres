#!/usr/bin/env python
#
# Garrett Granroth <granrothge@ornl.gov>
import unittest
from dgsres.singlextal.mdh_tools import gen_lattice_vectors,slice_from_MDH
import numpy as np
import os

class tests(unittest.TestCase):
    def test_latt_vec_ortho(self):
        vs = gen_lattice_vectors(6,7,8,90,90,90)
        self.assertTrue(np.all(np.abs(vs[0]-np.array([6,0,0]))/6<1e-6))
        self.assertTrue(np.all(np.abs(vs[1]-np.array([0,7,0]))/7<1e-6))
        self.assertTrue(np.all(np.abs(vs[2]-np.array([0,0,8]))/7<1e-6))
        
    def test_latt_vec_hex(self):
        vs = gen_lattice_vectors(6,6,10,90,90,120)
        self.assertTrue(np.all(np.abs(vs[0]-np.array([6,0,0]))/6<1e-6))
        self.assertTrue(np.all(np.abs(vs[1]-np.array([-6/2,6/2*np.sqrt(3),0]))/7<1e-6))
        self.assertTrue(np.all(np.abs(vs[2]-np.array([0,0,10]))/10<1e-6))
        
    def test_mono_beta(self):
        vs = gen_lattice_vectors(5.81955, 10.08404, 6.89595, 90, 106.221, 90)
        self.assertTrue(np.all(np.abs(vs[0]-np.array([5.81955, 0, 0]))/5.81955<1e-6))
        self.assertTrue(np.all(np.abs(vs[1]-np.array([0, 10.08404, 0]))/10.08404<1e-6))
        self.assertTrue(np.all(np.abs(vs[2]-np.array([-1.92634, 0, 6.62143]))/6.89595<1e-6))
    
    def test_tri(self):
        vs = gen_lattice_vectors(8.478, 9.572, 10.487, 73.46, 68.591, 64.762)
        self.assertTrue(np.all(np.abs(vs[0]-np.array([8.478, 0, 0]))/8.478<1e-4))
        self.assertTrue(np.all(np.abs(vs[1]-np.array([4.081, 8.658, 0]))/9.572<1e-4))
        self.assertTrue(np.all(np.abs(vs[2]-np.array([3.828, 1.496, 9.648]))/10.487<1e-4))

    def test_slice_from_MDH(self):
        MDH_path = os.path.join('SEQUOIA_data','slice_0p5K0E_28meV_4K.nxs')
        sl = slice_from_MDH(MDH_path,'test')
        self.assertTrue(sl.name=='test')
        self.assertTrue(np.all(sl.hkl_projection-np.array([0,1,0])<1e-8))
        self.assertTrue(np.all(sl.hkl0-np.array([0.5,0,0])<1e-8))
        self.assertTrue(abs(sl.grid.qaxis.min+3.975)<1e-6)
        self.assertTrue(abs(sl.grid.qaxis.max-3.975)<1e-6)
        self.assertTrue(abs(sl.grid.qaxis.step-0.05)<1e-6)
        self.assertTrue(abs(sl.grid.Eaxis.min+2.925)<1e-6)
        self.assertTrue(abs(sl.grid.Eaxis.max-26.925)<5e-6)
        self.assertTrue(abs(sl.grid.Eaxis.step-.15)<1e-6)

if __name__ == '__main__': unittest.main()