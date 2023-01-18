#!/usr/bin/env python
#
# Jiao Lin <jiao.lin@gmail.com>

import os, numpy as np, pytest
import mcvine.cli
from dgsres.singlextal import fit2icg
from dgsres import icg
here = os.path.abspath(os.path.dirname(__file__))

def test():
    cncs_geom = icg.Geom(l1=6.413, l2=36.2-6.413, l3=3.5)
    Ei = 12.
    E0 = 5.
    datadir = os.path.join(here, '..', 'data', 'PbTe-CNCS-HH4-Ei_12-E_5-q_0')
    fitter = fit2icg.Fit2ICG(datadir, cncs_geom, qaxis=(-0.3, 0.3, 0.01), Eaxis=(-2, 1, 0.05), Ei=Ei, E=E0)
    fitter.load_mcvine_psf_qE()
    fitted_params, res_x, res_y, yfit, E_profile = fitter.fit_E_profile(a=0.35, b=0.13, R=.8, sigma=2, t0=15)
    print(fitted_params)
    popt, res_x, res_y, yfit, q_profile = fitter.fit_q_profile(sigma=0.1)
    fitted_parameters, qgrid, Egrid, res_z, zfit, qE_profile = fitter.fit_qE_profile(dq_over_dE=0.02)
    print(fitted_parameters)
    return

def main():
    test()
    return


if __name__ == '__main__': main()

# End of file 
