#!/usr/bin/env python
#
# Jiao Lin <jiao.lin@gmail.com>

import mcvine.cli
import dgsres
from dgsres.singlextal import use_res_comps
import numpy as np, os
here = os.path.abspath(os.path.dirname(__file__))


def test_setup():
    def pixel_orientation_func(theta, phi):
        from mcni.neutron_coordinates_transformers.mcstasRotations import toMatrix
        if np.rad2deg(np.abs(theta))>34:
            m = np.dot(toMatrix(0, np.rad2deg(phi), 0), toMatrix(0, 0, 90.))
        else:
            m = np.dot(toMatrix(np.rad2deg(theta)-90., 0, 0), toMatrix(0, np.rad2deg(phi), 0))
        return m

    sampleyml = os.path.join(here, "Si.yml")
    # still use ARCS beam for now, just for testing purpose
    # beam = os.path.expanduser("~/beam/ARCS/100meV")
    beam = os.path.join(here, '..', 'data', 'beam', 'ARCS', '100meV')
    E = 40.
    hkl = [-16/3.,-8/3.,8/3.]
    hkl_projection = np.array([-1.,1.,-1.])/3
    psi_axis = dgsres.axis(min=-5, max=90., step=0.5)
    # this is not right -- using some ARCS parameters, just for testing purpose
    # the main point is we use "sphere" for detesys_shape
    instrument = use_res_comps.instrument(
        name = 'CHESS',
        detsys_radius = "3*meter", detsys_shape = 'sphere',
        L_m2s = "13.6*meter",
        offset_sample2beam = "-0.15*meter", # offset from sample to saved beam
        pixel_orientation_func = pixel_orientation_func
        )
    pixel = use_res_comps.pixel(
        radius = "0.5*inch",
        height = "1.5*meter/128",
        pressure = "6*atm",
        )
    outdir = os.path.join(here, '_tmp.out-chess')
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
