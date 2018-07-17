#!/usr/bin/env python
#
# Jiao Lin <jiao.lin@gmail.com>

import mcvine.cli
import dgsres
from dgsres.powder import use_ConstantQEKernel
import numpy as np, os
here = os.path.abspath(os.path.dirname(__file__))


def test():
    beam = "/SNS/users/lj7/simulations/ARCS/beam/100meV-n1e10"
    # for running on jenkins. data is downloaded by jenkins/getbeam.sh
    beam = os.path.expanduser("~/beam/ARCS/100meV")
    E = 50.
    Q = 5
    workdir = os.path.join(here, "work.useConstantQE")
    sim = use_ConstantQEKernel.Sim(
        instrument = 'ARCS',
        workdir = workdir,
        beamdir = beam,
        Ei = 100,
        Q = Q,
        dQ_axis=(-2, 2, 0.02),
        dE_axis=(-30, 30, 1.),
        ncount = 1e6,
        nodes = 5,
        )
    sim.run(E=E)
    return


def main():
    test()
    return


if __name__ == '__main__': main()

# End of file 
