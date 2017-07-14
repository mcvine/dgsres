#!/usr/bin/env python
#
# Jiao Lin <jiao.lin@gmail.com>

import mcvine.cli
from mcvine_workflow.singlextal.resolution import use_covmat

import numpy as np

def test():
    sampleyml = "Si.yml"
    Ei = 100
    class dynamics:
        hkl0 = [-16/3.,-8/3.,8/3.]
        hkl_dir = np.array([-1.,1.,-1.])/3
        E = 40.
        dq = 0
    class scan:
        min, max, step = -5, 90., 0.5
    instrument = use_covmat.instrument(
        name = 'ARCS',
        detsys_radius = "3.*meter",
        L_m2s = "13.6*meter",
        L_m2fc = "11.61*meter",
        offset_sample2beam = "-0.15*meter" # offset from sample to saved beam
        )
    pixel = use_covmat.pixel(
        radius = "0.5*inch",
        height = "meter/128",
        pressure = "10*atm",
        )
    tofwidths = use_covmat.tofwidths(P=10,M=8)
    beamdivs = use_covmat.beamdivs(theta=0.01, phi=0.01)
    samplethickness = 0.01
    use_covmat.compute(
        sampleyml, Ei, dynamics, scan,
        instrument, pixel,
        tofwidths, beamdivs, samplethickness,
        plot=True)
    return


def main():
    test()
    return


if __name__ == '__main__': main()

# End of file 
