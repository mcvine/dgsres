#!/usr/bin/env python
#
# Garrett Granroth <granrothge@ornl.gov>
import os, numpy as np
import unittest
from dgsres.singlextal.PSF_Affine_Model import gaus

from dgsres.singlextal.workflow import fit
import imp
import time
here = os.path.abspath(os.path.dirname(__file__))

st = time.time()
class TestCase(unittest.TestCase):
    def test_NO_cache(self):
        config = imp.load_source('config', os.path.join(here,'SEQUOIA_data','Res_config_file_29.py'))
        workdir = os.path.join(here,'SEQUOIA_data','')
        os.chdir(workdir)
        print("load source {}".format(time.time()-st))
        sl=config.slices[0]
        q = 0.5 
        E = 12.0
        #datadir = os.path.join(config.thisdir,config.simdir(q,E,config.slices[0]))
        fitter = fit(q,E,sl,config,use_cache=True)
        return
    def test_cache(self):
        config = imp.load_source('config', os.path.join(here,'SEQUOIA_data','Res_config_file_29.py'))
        workdir = os.path.join(here,'SEQUOIA_data','')
        os.chdir(workdir)
        print("load source {}".format(time.time()-st))
        sl=config.slices[0]
        q = 0.5 
        E = 12.0
        #datadir = os.path.join(config.thisdir,config.simdir(q,E,config.slices[0]))
        fitter = fit(q,E,sl,config,use_cache=True)
        return

if __name__ == '__main__': 
    unittest.main()

# End of file