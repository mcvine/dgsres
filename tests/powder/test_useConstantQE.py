#!/usr/bin/env python
#
# Jiao Lin <jiao.lin@gmail.com>

import mcvine.cli
import dgsres
from dgsres.powder import use_ConstantQEKernel

import numpy as np, os

def test():
    beam = "/SNS/users/lj7/simulations/ARCS/beam/100meV-n1e10"
    E = 80.
    Q = 7
    workdir = "work.useConstantQE"
    sim = use_ConstantQEKernel.Sim(
        workdir = workdir,
        beamdir = "/SNS/users/lj7/simulations/ARCS/beam/300meV-n3e9",
        Ei = 300,
        Q = Q,
        dQ_axis=(-2, 2, 0.02),
        dE_axis=(-60, 40, 1.),
        ncount = 1e7,
        nodes = 20,
        )
    sim.run(E=E)
    return


def main():
    test()
    return


if __name__ == '__main__': main()

# End of file 
