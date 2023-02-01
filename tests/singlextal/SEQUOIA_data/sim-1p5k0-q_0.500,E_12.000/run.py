#!/usr/bin/env python
import cloudpickle as pkl
import mcvine.cli
from numpy import array
from dgsres.singlextal import use_res_comps as urc
# parameters
beam_neutrons_path = '/home/19g/simulations/SEQUOIA/NiPS3/res_files_29/beam/out/neutrons'
instrument = pkl.load(open('/home/19g/simulations/SEQUOIA/NiPS3/res_files_29/sim-1p5k0-q_0.500,E_12.000/instrumentk6ixy5sm.pkl', 'rb'))
samplexmlpath = '/home/19g/simulations/SEQUOIA/NiPS3/res_files_29/sim-1p5k0-q_0.500,E_12.000/sample/sampleassembly.xml'
psi = 0.6586104396148563
hkl2Q = array([[ 0.90043725, -0.66221495, -0.12259223],
       [-0.36971606, -0.5014957 , -0.00659223],
       [ 0.13799321, -0.08944738, -0.93455866]])
pixel = pkl.load(open('/home/19g/simulations/SEQUOIA/NiPS3/res_files_29/sim-1p5k0-q_0.500,E_12.000/pixel41bqx17t.pkl', 'rb'))
t_m2p = 0.011549456114227113
Q = array([ 1.16579784, -1.24407028, -0.18718446])
E = 12.0
hkl_projection = array([0, 1, 0])
# mc parameters
mc_p_path = './mc_params.yml'
import yaml, os
if os.path.exists(mc_p_path):
    mc_params = yaml.safe_load(open(mc_p_path))
else:
    mc_params = dict(Nbuffer=10000, Nrounds_beam=1)
Nbuffer = 100000
Nrounds_beam = 1
# run
urc.run(
    beam_neutrons_path, instrument, samplexmlpath, psi, hkl2Q, pixel, t_m2p,
    Q, E, hkl_projection, **mc_params)
