#!/usr/bin/env python
#
# Garrett Granroth <granrothge@ornl.gov>
import os, numpy as np, lmfit
from dgsres.singlextal.PSF_Affine_Model import gaus
from dgsres.singlextal.fit2ee import EEModel
from dgsres.singlextal.fit_2d_psf import Fit 
import imp
here = os.path.abspath(os.path.dirname(__file__))
def test():
    config = imp.load_source('config', os.path.join(here,'SEQUOIA_data','Res_config_file_29.py'))
    datadir = os.path.join(config.thisdir,config.simdir(0.5,12.0,config.slices[0]))
    qaxis = config.slices[0].res_2d_grid.qaxis
    Eaxis = config.slices[0].res_2d_grid.Eaxis
    Fitobj = Fit(datadir,
                qaxis=(qaxis.min, qaxis.max, qaxis.step),
                Eaxis=(Eaxis.min, Eaxis.max, Eaxis.step),
    Ei=config.Ei, E=12.0)
    Fitobj.load_mcvine_psf_qE(adjust_energy_center=True)
def main():
    test()
    return


if __name__ == '__main__': main()

# End of file