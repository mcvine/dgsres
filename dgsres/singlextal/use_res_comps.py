# -*- Python -*-
#
# Jiao Lin <jiao.lin@gmail.com>
#

"""
Tools to run mcvine simulation that computes the resolution function.
The idea is to first calculate the position of the pixel that corresponds
to a point hklE and also create a "resolution sample" that scatters
neutrons to that pixel.
Then a mcvine simulation is done using precalculated beam,  the sample,
and the pixel.
The simulated neutrons are gathered into events with hklE coordinates,
which can be hisogrammed to resolution function.

The procedure is to 
* Setup the simulation
* Run the simulation

MC parameters:
* MC parameters are used to control the Monte Carlo aspect of the simulation
* They can be independently set by an yml file called mc_params.yml. Parameters:
  - Nbuffer
  - Nrounds_beam
* The file needs to be present in the simulation directory ("outdir" of setup method)

Limitations:
* assume cylindrical arrangement of det tubes for the det system

Inputs:

* sample assembly
* beam
* dynamics
  - from these compute Q in instrument coordinate system (z vertical up)
    and also hkl2Q matrix
* instrument geometry info
* pixel
"""

import os, sys, numpy as np, cloudpickle as pkl, tempfile as tpf
import mcvine
from mcni.utils import conversion as Conv


# input data structures
from . import instrument, pixel


def setup(outdir, sampleyml, beam, E, hkl, hkl_projection, psi_axis, instrument, pixel, log=None):
    """setup the simulation directory and scripts and input files

    - outdir: output dir
    - sampleyml: path. should contain
      * name
      * chemical_formula
      * lattice
      * orientation
      * shape
    - beam: mcvine beam simulation path
    - E: energy transfer
    - hkl: momentum transfer
    - hkl_projection: hkl projection axis for slice
    - psi_axis: psi scan
    - instrument: instrument object with geometry data such as L1 and L2
    - pixel: pixel object with pixe, height, pressure. and position
    - log: logger object
    """
    if log is None:
        import sys
        log = sys.stdout
    # load beam
    from ._beam import computeEi_and_t0
    Ei, t0 = computeEi_and_t0(beam, instrument)
    log.write( "Ei=%s, t0=%s\n" % (Ei, t0) )
    # load sample
    from mcvine.workflow.sample import loadSampleYml
    sample = loadSampleYml(sampleyml)
    # the sample kernel need information of E and hkl
    Q, hkl2Qmat, psi = calcQ(sampleyml, Ei, E, hkl, psi_axis, Npsisegments=10)
    log.write( "Computed:\n" )
    log.write( "* psi=%s degree\n" % (psi/np.pi*180,) )
    log.write( "* Q=%s\n" % (Q,) )
    log.write( "* hkl2Qmat=%s\n" % (hkl2Qmat,) )
    kfv, Ef = computeKf(Ei, E, Q, log)
    log.write( "* Ef=%s\n" % (Ef,))
    pixel_position, pixel_orientation = computePixelPositionOrientation(kfv, instrument, log)
    # at this point the coordinates have convention of z vertical up
    # ** coordinate system for calculated position: z is vertical **
    # this pixel_position is in the instrument coordinate system.
    # we need, however, the pixel_position in the instrument
    # coordinate system rorated psi angle.
    from numpy import sin, cos
    x,y,z=pixel_position
    x2,y2,z2 = x*cos(psi)+y*sin(psi), -x*sin(psi)+y*cos(psi), z
    # convert to z along beam
    pixel_position3 = y2, z2, x2
    # compute nominal tof from mod to sample
    vi = Conv.e2v(Ei)
    L_m2s = mcvine.units.parse(instrument.L_m2s)/mcvine.units.meter
    t_m2s = L_m2s/vi + t0*1e-6
    # nominal tof from mod to pixel
    vf = Conv.e2v(Ef)
    t_s2p = np.linalg.norm(pixel_position)/vf
    t_m2p = t_m2s + t_s2p
    log.write( "t_m2s=%s, t_s2p=%s, t_m2p=%s\n" % (t_m2s, t_s2p, t_m2p))
    # tof passing through the pixel
    r = mcvine.units.parse(pixel.radius)/mcvine.units.meter
    h = mcvine.units.parse(pixel.height)/mcvine.units.meter
    dtof = np.sqrt(4*r*r+h*h)*1.1/vf
    # decorate the sample kernel
    kernel = sample.excitations[0]
    kernel.target_position = "%s*meter,%s*meter,%s*meter" % pixel_position3
    kernel.target_radius = "%s*meter" % (np.sqrt(r*r+h*h/4)*1.1,)
    kernel.tof_at_target = "%s*microsecond" % (t_m2p*1e6)
    kernel.dtof = "%s*microsecond" % (dtof*1e6,)
    # create sample assembly
    from mcvine.workflow.singlextal.scaffolding import createSampleAssembly
    sampledir = os.path.abspath(os.path.join(outdir, 'sample'))
    createSampleAssembly(sampledir, sample, add_elastic_line=False)
    # save instrument object
    ofstream = tpf.NamedTemporaryFile(
        dir=outdir, prefix='instrument', suffix='.pkl', delete=False)
    pkl.dump(instrument, ofstream)
    instr_fn = os.path.abspath(ofstream.name)
    ofstream.close()
    # save pixel object
    pixel.position = pixel_position
    pixel.orientation = pixel_orientation
    ofstream = tpf.NamedTemporaryFile(
        dir=outdir, prefix='pixel', suffix='.pkl', delete=False)
    pkl.dump(pixel, ofstream)
    pixel_fn = os.path.abspath(ofstream.name)
    ofstream.close()
    # create sim script
    params = dict(
        beam_neutrons_path = os.path.join(beam, 'out', 'neutrons'),
        instr_fn = instr_fn, pixel_fn = pixel_fn,
        samplexmlpath = os.path.join(sampledir, "sampleassembly.xml"),
        psi = psi,
        hkl2Q = hkl2Qmat,
        t_m2p = t_m2p,
        Q = Q,
        E = E,
        hkl_projection = hkl_projection,
        )
    script = sim_script_template % params
    open(os.path.join(outdir, 'run.py'), 'wt').write(script)
    return
sim_script_template = """#!/usr/bin/env python
import cloudpickle as pkl
import mcvine.cli
from numpy import array
from dgsres.singlextal import use_res_comps as urc
# parameters
beam_neutrons_path = %(beam_neutrons_path)r
instrument = pkl.load(open(%(instr_fn)r, 'rb'))
samplexmlpath = %(samplexmlpath)r
psi = %(psi)r
hkl2Q = %(hkl2Q)r
pixel = pkl.load(open(%(pixel_fn)r, 'rb'))
t_m2p = %(t_m2p)r
Q = %(Q)r
E = %(E)r
hkl_projection = %(hkl_projection)r
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
""" 

def computeKf(Ei, E, Q, log):
    ki = Conv.e2k(Ei);
    log.write( "* ki=%s\n" % (ki,) )
    kiv = np.array([ki, 0, 0])
    kfv = kiv - Q
    log.write( "* vectors ki=%s, kf=%s\n" % (kiv, kfv) )
    Ef = Ei - E
    # ** Verify the momentum and energy transfers **
    log.write( "These two numbers should be very close:\n")
    log.write( "  %s\n" % (Ei-Conv.k2e(np.linalg.norm(kfv)),) )
    log.write( "  %s\n" % (Ei-Ef,) )
    assert np.isclose(Ef, Conv.k2e(np.linalg.norm(kfv)))
    log.write( "  Ei=%s, Ef=%s\n" % (Ei,Ef) )
    return kfv, Ef

def computePixelPositionOrientation(kfv, instrument, log):
    # ** compute nominal TOF at detector pixel **
    # where is detector pixel?
    # kfv is in instrument scientist coordinate system
    R = mcvine.units.parse(instrument.detsys_radius)/mcvine.units.meter
    kf = np.linalg.norm(kfv)
    if instrument.detsys_shape == 'cylinder':
        t_sample2pixel = R/(kfv[0]**2 + kfv[1]**2)**.5
        rotmat = np.eye(3)
    elif instrument.detsys_shape == 'sphere':
        t_sample2pixel = R/kf
        cos_theta = kfv[2]/kf; theta = np.arccos(cos_theta)
        phi = np.arctan2(kfv[1], kfv[0])
        rotmat = instrument.pixel_orientation_func(theta, phi)
    else:
        raise NotImplementedError("detector system shape {}".format(instrument.detsys_shape))
    pixel_pos = kfv*t_sample2pixel
    from mcni.neutron_coordinates_transformers.mcstasRotations import toAngles
    pixel_ori = toAngles(rotmat, unit='degree')
    log.write( "* pixel positon=%s orientation=%s\n" % (pixel_pos, pixel_ori))
    return pixel_pos, pixel_ori

def calcQ(sampleyml, Ei, E, hkl, psi_axis, Npsisegments=10):
    from mcvine.workflow.singlextal.io import loadXtalOriFromSampleYml
    xtalori = loadXtalOriFromSampleYml(sampleyml)
    psimin, psimax, dpsi = psi_axis.min, psi_axis.max, psi_axis.step
    from mcvine.workflow.singlextal.solve_psi import solve
    results = solve(
        xtalori, Ei, hkl, E, psi_min=psimin, psi_max=psimax,
        Nsegments = Npsisegments)
    assert len(results)
    from mcvine.workflow.singlextal.coords_transform import hkl2Q
    r0 = results[0]
    psi_in_degrees = r0
    xtalori.psi = psi = r0*np.pi/180.
    Q = hkl2Q(hkl, xtalori) # here the convension is z vertical
    hkl2Qmat = xtalori.hkl2cartesian_mat()
    return Q, hkl2Qmat, psi


def run(beam_neutrons_path, instrument, samplexmlpath, psi, hkl2Q, pixel, t_m2p,
        Q, E, hkl_projection, Nbuffer=100000, Nrounds_beam=1):
    """Run mcvine simulation with using the resolution sample and the resolution pixel,
    and save the results as numpy arrays

    - beam_neutrons_path: path to the "neutrons" file of a beam simulation
    - instrument: instrument object with geometry data such as L1 and L2
    - samplexmlpath: path to the sampleassembly xml file
    - psi: sample rotation angle
    - hkl2Q: matrix to convert hkl to Q: Q = hkl dot hkl2Q
    - pixel: pixel object with pixe, height, pressure. and position
    - t_m2p: exepcted tof from moderator to pixel
    - Q: expected Q
    - E: expected E
    - hkl_projection: the simulated data is projected to this axis and E axis for easy inspection
    - Nbuffer: neutron buffer size
    - Nrounds_beam: number of rounds replaying neutrons in the beam
    """
    from mcni.components.NeutronFromStorage import NeutronFromStorage
    source = NeutronFromStorage('source', path=beam_neutrons_path)
    from mccomponents.sample import samplecomponent
    sample = samplecomponent( 'sample', samplexmlpath)
    # dummy component to save the state of scattered neutrons
    from mcni.components.Dummy import Dummy
    sample_location = Dummy('sample_location')
    # det pixel
    from mccomponents.components.DGSSXResPixel import DGSSXResPixel
    pressure = mcvine.units.parse(pixel.pressure)/mcvine.units.parse("kg/m/s/s")
    r = mcvine.units.parse(pixel.radius)/mcvine.units.meter
    h = mcvine.units.parse(pixel.height)/mcvine.units.meter
    pixel_comp = DGSSXResPixel(
        "pixel",
        pressure=pressure, tof=t_m2p,
        radius=r, height=h)
    # build instrument simulation chain
    import mcni
    sim_chain = mcni.instrument( [source, sample, sample_location, pixel_comp] )
    # put components into place
    geometer = mcni.geometer()
    z_beam = mcvine.units.parse(instrument.offset_sample2beam)/mcvine.units.meter
    # z along beam
    geometer.register( source, (0,0,z_beam), (0,0,0) )
    geometer.register( sample, (0,0,0), (0,psi*180/np.pi,0) )
    geometer.register( sample_location, (0,0,0), (0,0,0) )
    geometer.register( pixel_comp, pixel.position, pixel.orientation )
    #
    Q2hkl = np.linalg.inv(hkl2Q)
    # lengh of hkl_projection squared
    hkl_proj_len2 = np.dot(hkl_projection, hkl_projection)
    # neutron buffer
    from mcni.neutron_storage.idf_usenumpy import count
    N0 = count(beam_neutrons_path) * Nrounds_beam
    dxs_all = None; dEs_all = None; probs_all=None; dhkls_all=None
    start = 0
    for i in range(int(np.ceil(N0/Nbuffer))+1):
    # for i in range(10):
        end = start + Nbuffer
        end = min(end, N0)
        if end<=start:
            continue
        sys.stdout.write("%s-%s: " % (start, end-1)); sys.stdout.flush()
        neutrons = mcni.neutron_buffer(end-start)
        # simulate
        tracer = NeutronTracer()
        mcni.simulate( sim_chain, geometer, neutrons, tracer=tracer)
        #
        before_dummy_start, after_dummy_start, \
            before_incident, after_incident, \
            before_scattered, after_scattered, \
            at_sample_location, at_sample_location2, \
            before_detected, after_detected, \
            before_dummy_end, after_dummy_end = tracer._store
        incident = after_incident
        # has to be `at_sample_location` because both sample_location and
        # beam have no relative rotation. 
        # should not used after_scattered.
        # it is in the sample's coordinate system, which is rotated by angle psi.
        # should not use before_detected. It could be rotated in spherical case
        scattered = at_sample_location
        detected = after_detected
        del (before_dummy_start, after_dummy_start,
             before_incident, after_incident,
             before_scattered, after_scattered,
             before_detected, after_detected,
             before_dummy_end, after_dummy_end)
        is_scattered = incident.v != scattered.v
        is_scattered = np.logical_or(
            is_scattered[:,0],
            np.logical_or(is_scattered[:,1], is_scattered[:,2])
            )
        good = np.logical_and(is_scattered, detected.p>np.finfo(float).eps)

        vi = incident.v[good]
        vf = scattered.v[good]
        probs = p = detected.p[good]

        from mcni.utils import conversion
        Ei = conversion.VS2E * (vi*vi).sum(axis=-1)
        Ef = conversion.VS2E * (vf*vf).sum(axis=-1)
        Es = Ei - Ef

        vQ = vi-vf
        Qs = vQ * conversion.V2K
        Qs = np.array([Qs[:, 2], Qs[:, 0], Qs[:, 1]]).T
        # print Qs
        # print Es
        dQs = Qs - Q
        dEs = Es - E
        dhkls = np.dot(dQs, Q2hkl)
        # print dhkls
        # print dEs
        # print p
        dxs = np.dot( dhkls, np.array(hkl_projection) )/hkl_proj_len2
        if dxs_all is None:
            dxs_all = dxs
            dEs_all = dEs
            probs_all = probs
            dhkls_all = dhkls
        else:
            dxs_all = np.concatenate((dxs_all, dxs))
            dEs_all = np.concatenate((dEs_all, dEs))
            probs_all = np.concatenate((probs_all, probs))
            dhkls_all = np.concatenate((dhkls_all, dhkls))
        print()
        start = end
        continue
    # reverse x and E
    # the negative sign here makes the result the PSF
    dxs_all *= -1
    dEs_all *= -1
    dhkls_all *= -1
    # save results
    np.save("dhkls.npy", dhkls_all)
    np.save("dxs.npy", dxs_all)
    np.save("dEs.npy", dEs_all)
    np.save("probs.npy", probs_all)
    h, xedges, yedges = np.histogram2d(
        dxs_all, dEs_all, bins=100, weights=probs_all)
    import histogram as H, histogram.hdf as hh
    xaxis = H.axis('x', boundaries=xedges)
    Eaxis = H.axis('E', boundaries=yedges)
    res = H.histogram('res', (xaxis, Eaxis), data=h)
    hh.dump(res, 'res.h5')
    sys.stdout.write("Done.\n"); sys.stdout.flush()
    return


def main():
    # test_process()
    run()
    return


class _N:
    def __init__(self, r, v, p):
        self.r, self.v, self.p = r,v,p
        return

# custom tracer
from mcni.pyre_support.AbstractNeutronTracer import AbstractNeutronTracer
class NeutronTracer(AbstractNeutronTracer):

    def __init__(self, name='neutron-tracer'):
        super(NeutronTracer, self).__init__(name)
        self._store = []
        return

    
    def __call__(self, neutrons, context=None):
        if context:
            context.identify(self)

        from mcni.neutron_storage import neutrons_as_npyarr
        a = neutrons_as_npyarr(neutrons)
        a.shape = -1, 10
        r = a[:, :3]
        v = a[:, 3:6]
        p = a[:, -1]
        self._store.append(_N(r, v, p))
        sys.stdout.write(".")
        sys.stdout.flush()
        return


    def onBefore(self, context):
        return

    def onProcessed(self, context):
        return

# End of file 
