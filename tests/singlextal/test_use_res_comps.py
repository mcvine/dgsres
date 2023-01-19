#!/usr/bin/env python
#
# Jiao Lin <jiao.lin@gmail.com>

import mcvine.cli
import dgsres
from dgsres.singlextal import use_res_comps
import numpy as np, os
here = os.path.abspath(os.path.dirname(__file__))


def test_setup():
    sampleyml = os.path.join(here, "Si.yml")
    beam = "/SNS/users/lj7/simulations/ARCS/beam/100meV-n1e10"
    # for running on jenkins. data is downloaded by jenkins/getbeam.sh
    # beam = os.path.expanduser("~/beam/ARCS/100meV")
    beam = os.path.join(here, '..', 'data', 'beam', 'ARCS', '100meV')
    E = 40.
    hkl = [-16/3.,-8/3.,8/3.]
    hkl_projection = np.array([-1.,1.,-1.])/3
    psi_axis = dgsres.axis(min=-5, max=90., step=0.5)
    instrument = use_res_comps.instrument(
        name = 'ARCS',
        detsys_radius = "3.*meter",
        L_m2s = "13.6*meter",
        offset_sample2beam = "-0.15*meter" # offset from sample to saved beam
        )
    pixel = use_res_comps.pixel(
        radius = "0.5*inch",
        height = "meter/128",
        pressure = "10*atm",
        )
    outdir = os.path.join(here, '_tmp.out')
    if os.path.exists(outdir):
        import shutil
        shutil.rmtree(outdir)
    use_res_comps.setup(
        outdir, sampleyml, beam, E, hkl, hkl_projection,
        psi_axis, instrument, pixel)
    return


def main():
    test_setup()
    return


if __name__ == '__main__': main()

# End of file 
