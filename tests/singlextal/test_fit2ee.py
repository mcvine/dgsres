#!/usr/bin/env python
#
# Jiao Lin <jiao.lin@gmail.com>

import os, numpy as np
import mcvine.cli
from dgsres.singlextal import fit2ee
here = os.path.abspath(os.path.dirname(__file__))

@pytest.mark.skipif(True, reason='temporarily disabled')
def test():
    Ei = 12.
    E0 = 5.
    datadir = os.path.join(here, '..', 'data', 'PbTe-CNCS-HH4-Ei_12-E_5-q_0')
    fitter = fit2ee.Fit2EE(datadir, qaxis=(-0.3, 0.3, 0.01), Eaxis=(-2, 1, 0.05), Ei=Ei, E=E0)
    fitter.load_mcvine_psf_qE()
    fitted_params, res_x, res_y, yfit, E_profile = fitter.fit_E_profile(
        sigma_left=0.2, sigma_right=0.2, ef_width=0.1, weight_left=0.5, E_offset=0.)
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
