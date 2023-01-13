#!/usr/bin/env python
#
# Garrett Granroth <granrothge@ornl.gov>
import os, numpy as np, lmfit
from dgsres.singlextal.PSF_Affine_Model import gaus
from dgsres.singlextal.fit2ee import EEModel
from dgsres.singlextal.fit_2d_psf import Fit 
import imp
import time
here = os.path.abspath(os.path.dirname(__file__))
st = time.time()
def test():
    config = imp.load_source('config', os.path.join(here,'SEQUOIA_data','Res_config_file_29.py'))
    print("load source {}".format(time.time()-st))
    sl=config.slices[0]
    datadir = os.path.join(config.thisdir,config.simdir(0.5,12.0,config.slices[0]))
    qaxis = sl.res_2d_grid.qaxis
    Eaxis = sl.res_2d_grid.Eaxis
    Fitobj = Fit(datadir,
                qaxis=(qaxis.min, qaxis.max, qaxis.step),
                Eaxis=(Eaxis.min, Eaxis.max, Eaxis.step),
    Ei=config.Ei, E=12.0)
    Fitobj.load_mcvine_psf_qE(adjust_energy_center=True)
    print("load psf {}".format(time.time()-st))
    fitting_params = dict([(k,v) for k,v in sl.fitting.__dict__.items() if not k.startswith('_')])
    print("populate fitting parameters {}".format(time.time()-st))
def main():
    test()
    return


if __name__ == '__main__': main()

# End of file