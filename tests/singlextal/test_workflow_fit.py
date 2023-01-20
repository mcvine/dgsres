#!/usr/bin/env python
#
# Garrett Granroth <granrothge@ornl.gov>
import os, numpy as np
import pytest
from dgsres.singlextal.PSF_Affine_Model import gaus

from dgsres.singlextal.workflow import fit
import imp
import time
import cloudpickle as pkl
interactive=False

ts = time.time()

@pytest.fixture
def init_dec():
    here = os.path.abspath(os.path.dirname(__file__))
    config = imp.load_source('config', os.path.join(here,'SEQUOIA_data','Res_config_file_29.py'))
    workdir = os.path.join(here,'SEQUOIA_data','')
    with open(os.path.join(here,'SEQUOIA_data','cnf-1p5k0-q_0.500-E_12.000-fitter.pkl'), 'rb') as fh:
        fit_chk = pkl.load(fh)
    os.chdir(workdir)
    sl = config.slices[0]
    q = 0.5
    E = 12.0
    print("time = {} s".format(time.time()-ts))
    return [q,E,sl,config,fit_chk]



def test_NO_cache(init_dec):
    fp = init_dec[:-1]
    fit_chk = init_dec[-1]
    fitter = fit(*fp,use_cache=False)
    tnum = np.sum((fitter.fit_result.best_fit-fit_chk.fit_result.best_fit)**2/fitter.fit_result.best_fit**2)
    tdnm = np.sum(1/fitter.fit_result.best_fit**2)
    print (np.sqrt(tnum/tdnm))
    assert np.sqrt(tnum/tdnm) < 1e-8
    return


def test_cache(init_dec):
    fp = init_dec[:-1]
    fit_chk = init_dec[-1]
    fitter = fit(*fp,use_cache=True)
    tnum = np.sum((fitter.fit_result.best_fit-fit_chk.fit_result.best_fit)**2/fitter.fit_result.best_fit**2)
    tdnm = np.sum(1/fitter.fit_result.best_fit**2)
    print (np.sqrt(tnum/tdnm))
    assert np.sqrt(tnum/tdnm) < 1e-8
    os.remove('1p5k0-q_0.500-E_12.000-fitter.pkl')
    return


def test_NO_multiprocessor(init_dec):
    fp = init_dec[:-1]
    fit_chk = init_dec[-1]
    fitter = fit(*fp,use_cache=False,multiprocess=False)
    tnum = np.sum((fitter.fit_result.best_fit-fit_chk.fit_result.best_fit)**2/fitter.fit_result.best_fit**2)
    tdnm = np.sum(1/fitter.fit_result.best_fit**2)
    print (np.sqrt(tnum/tdnm))
    assert np.sqrt(tnum/tdnm) < 1e-8
    return



# End of file
