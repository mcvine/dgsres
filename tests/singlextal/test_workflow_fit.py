#!/usr/bin/env python
#
# Garrett Granroth <granrothge@ornl.gov>
import os, numpy as np
import unittest
from dgsres.singlextal.PSF_Affine_Model import gaus

from dgsres.singlextal.workflow import fit
import imp
import time
import cloudpickle as pkl
interactive=False
here = os.path.abspath(os.path.dirname(__file__))
ts = time.time()

def init_dec(func):
    config = imp.load_source('config', os.path.join(here,'SEQUOIA_data','Res_config_file_29.py'))
    workdir = os.path.join(here,'SEQUOIA_data','')
    with open(os.path.join(here,'SEQUOIA_data','cnf-1p5k0-q_0.500-E_12.000-fitter.pkl'), 'rb') as fh:
        fit_chk = pkl.load(fh)
    os.chdir(workdir)
    sl = config.slices[0]
    q = 0.5
    E = 12.0
    func(q,E,sl,config,fit_chk)
    print("time = {} s".format(time.time()-ts))


@init_dec
def test_NO_cache(q,E,sl,config,fit_chk):
    #datadir = os.path.join(config.thisdir,config.simdir(q,E,config.slices[0]))
    fitter = fit(q,E,sl,config,use_cache=False)
    tnum = np.sum((fitter.fit_result.best_fit-fit_chk.fit_result.best_fit)**2/fitter.fit_result.best_fit**2)
    tdnm = np.sum(1/fitter.fit_result.best_fit**2)
    print (np.sqrt(tnum/tdnm))
    assert np.sqrt(tnum/tdnm) < 1e-8
    return

@init_dec
def test_cache(q,E,sl,config,fit_chk):
    fitter = fit(q,E,sl,config,use_cache=True)
    tnum = np.sum((fitter.fit_result.best_fit-fit_chk.fit_result.best_fit)**2/fitter.fit_result.best_fit**2)
    tdnm = np.sum(1/fitter.fit_result.best_fit**2)
    print (np.sqrt(tnum/tdnm))
    assert np.sqrt(tnum/tdnm) < 1e-8
    os.remove('1p5k0-q_0.500-E_12.000-fitter.pkl')
    return

@init_dec
def test_NO_multiprocessor(q,E,sl,config,fit_chk):
    fitter = fit(q,E,sl,config,use_cache=False,multiprocess=False)
    tnum = np.sum((fitter.fit_result.best_fit-fit_chk.fit_result.best_fit)**2/fitter.fit_result.best_fit**2)
    tdnm = np.sum(1/fitter.fit_result.best_fit**2)
    print (np.sqrt(tnum/tdnm))
    assert np.sqrt(tnum/tdnm) < 1e-8
    return

def main():
    test_NO_cache()
    test_cache()
    test_NO_multiprocessor()

#if __name__ == '__main__': main()

# End of file
