#!/usr/bin/env python
#
# Garrett Granroth <granrothge@ornl.gov>
import os, numpy as np, lmfit
from dgsres.singlextal.PSF_Affine_Model import gaus
from dgsres.singlextal.fit2ee import EEModel
from dgsres.singlextal.fit_2d_psf import Fit 
from dgsres.singlextal.workflow import fit
import imp
import time
here = os.path.abspath(os.path.dirname(__file__))

st = time.time()
def test():
    config = imp.load_source('config', os.path.join(here,'SEQUOIA_data','Res_config_file_29.py'))
    workdir = os.path.join(here,'SEQUOIA_data','')
    os.chdir(workdir)
    print("load source {}".format(time.time()-st))
    sl=config.slices[0]
    q = 0.5 
    E = 12.0
    #datadir = os.path.join(config.thisdir,config.simdir(q,E,config.slices[0]))
    fitter = fit(q,E,sl,config,use_cache=True)
def main():
    test()
    return


if __name__ == '__main__': main()

# End of file