#!/usr/bin/env python
#
# Garrett Granroth <granrothge@ornl.gov>
import unittest
import os
import dgsres.singlextal.mdh_tools as mdht
from mcvine.workflow import singlextal as sx
from dgsres.instruments import sequoia
from dgsres.singlextal import workflow
from dgsres.singlextal import fit_ellipsoid
from dgsres.singlextal.sim_config import config_cls
import numpy as np
import shutil
here = os.path.dirname(os.path.abspath('__file__'))
project_root = os.path.abspath(os.path.join(here, '..', '..'))


class tests(unittest.TestCase):
    def test_config(self):
      
        MDH_path = os.path.join(project_root, 'tests', 'data', 'SEQUOIA_data',
                                'slice_0p5K0E_28meV_4K.nxs')
        beam_path = os.path.join(project_root, 'tests', 'data', 'SEQUOIA_data',
                                 'beam', '')
        Ei = 28.0
        workdir = './test_workflow'
        try:
            os.listdir(workdir)
        except FileNotFoundError:
            os.mkdir(workdir)
        os.chdir(workdir)
        config = config_cls(Ei=Ei, beampth=os.path.abspath(beam_path))
        config.thisdir = os.path.abspath(os.path.curdir)
        config.psi_scan = mdht.angles_from_MDH(MDH_path)
        config.instrument = sequoia
        sl = mdht.slice_from_MDH(MDH_path, 'test_slice')
        self.assertTrue(np.all(np.abs(sl.hkl0 - np.array([0.5, 0., 0.])) < 1e-6))
        self.assertTrue(np.all(np.abs(sl.hkl_projection - np.array([0.0, 1.0, 0.0])) < 1e-6))
        self.assertTrue(np.abs(config.Ei - 28.0) < 1e-6)
        self.assertTrue(np.abs(config.psi_scan.min - 0.0) < 1e-6)
        self.assertTrue(np.abs(config.psi_scan.max - 20.0) < 1e-6)
        self.assertTrue(np.abs(config.psi_scan.step - 2.0) < 1e-6)
        os.chdir('..')

    def test_sim_n_fit(self):
        MDH_path = os.path.join(project_root, 'tests', 'data', 'SEQUOIA_data',
                                'slice_0p5K0E_28meV_4K.nxs')
        beam_path = os.path.join(project_root, 'tests', 'data', 'SEQUOIA_data',
                                 'beam', '')
        Ei = 28.0
        workdir = './test_workflow'
        try:
            os.listdir(workdir)
        except FileNotFoundError:
            os.mkdir(workdir)
        os.chdir(workdir)
        config = config_cls(Ei=Ei, beampth=os.path.abspath(beam_path))
        config.thisdir = os.path.abspath(os.path.curdir)
        yml_path = os.path.join(config.thisdir, 'test.yml')
        config.psi_scan = mdht.angles_from_MDH(MDH_path)
        config.instrument = sequoia
        sl = mdht.slice_from_MDH(MDH_path, 'test_slice')
        mdht.sample_from_MDH(MDH_path, yml_file=yml_path)
        config.sample_yaml = yml_path
        self.assertTrue(os.path.exists(config.sample_yaml))
        sl.res_2d_grid.qaxis = sx.axis(min=-0.2, max=0.2, step=0.005)
        sl.res_2d_grid.Eaxis = sx.axis(min=-2.5, max=2.5, step=0.1)
        sl.fitting = mdht.fitting()
        sl.grid.Eaxis.min = 21.5
        sl.grid.Eaxis.max = 26.5
        sl.grid.Eaxis.step = 1.0
        sl.grid.qaxis.min = -4
        sl.grid.qaxis.max = -3
        sl.grid.qaxis.step = 0.2
        config.slices.append(sl)
        output, failed = workflow.simulate_all_in_one(config)

        workflow.print_sim_all_in_one_output(output, failed)

        workflow.fit_all_in_one(config)
        os.chdir('..')
        shutil.rmtree(workdir)
