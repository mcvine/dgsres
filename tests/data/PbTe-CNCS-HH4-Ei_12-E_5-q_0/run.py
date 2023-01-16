#!/usr/bin/env python
import mcvine.cli
from numpy import array
from dgsres.singlextal import use_res_comps as urc
# parameters
beam_neutrons_path = '/SNS/users/lj7/simulations/CNCS/AgBiSe2/PbTe-try1/beam/out/neutrons'
instrument = urc.instrument('CNCS', '3.5*meter', '36.264*meter', '-0.15*meter')
samplexmlpath = '/SNS/users/lj7/simulations/CNCS/AgBiSe2/PbTe-try1/test/sample/sampleassembly.xml'
psi = -0.3696248467562983
hkl2Q = array([[  2.48422992e-01,   6.41204719e-01,  -6.87646330e-01],
       [  2.48422992e-01,   6.41204719e-01,   6.87646330e-01],
       [  9.06800410e-01,  -3.51323165e-01,   5.03375347e-17]])
pp = array([ -2.28235251e+00,   2.65346321e+00,  -3.80187844e-16])
pixel = urc.pixel('0.5*inch', '2.*meter/128', '6*atm', position=(pp[1], pp[2], pp[0]))
t_m2p = 0.026950797917531365
Q = array([  3.62720164e+00,  -1.40529266e+00,   2.01350139e-16])
E = 5
hkl_projection = array([1, 1, 0])
# mc parameters
mc_p_path = './mc_params.yml'
import yaml, os
if os.path.exists(mc_p_path):
    mc_params = yaml.load(open(mc_p_path))
else:
    mc_params = dict(Nbuffer=10000, Nrounds_beam=1)
Nbuffer = 100000
Nrounds_beam = 1
# run
urc.run(
    beam_neutrons_path, instrument, samplexmlpath, psi, hkl2Q, pixel, t_m2p,
    Q, E, hkl_projection, **mc_params)
