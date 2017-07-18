# -*- Python -*-
#
# Jiao Lin <jiao.lin@gmail.com>
#

"""
Limitations:
* assume cylindrical arrangement of det tubes for the det system

Inputs

* sample assembly
  - shape
  - material (??.xyz)
  - pixel position
  - pixel radius
  - tof at pixel
  - dtof
* beam_path
  - a quick simulation to generate event nexus file and let mantid compute t0 and Ei?
    - t0
    - Ei
* dynamics
  - E
  - hkl
  - sample.yml
  - from these compute Q in instrument coordinate system (z vertical up)
    and also hkl2Q matrix
* L_m2s
* pixel_pos (this was calculated using R=3, cylinder)
* pixel_pressure, radius, height
* distance from saved neutrons to sample position
* Nbuffer
"""

import os, sys, numpy as np, mcvine
from mcni.utils import conversion as Conv


# input data structures
from . import instrument, pixel


def setup(outdir, sampleyml, beam, E, hkl, hkl_projection, psi_axis, instrument, pixel, log=None):
    if log is None:
        import sys
        log = sys.stdout
    # load beam
    from ._beam import computeEi_and_t0
    Ei, t0 = computeEi_and_t0(beam, instrument)
    log.write( "Ei=%s, t0=%s\n" % (Ei, t0) )
    # load sample
    from mcvine_workflow.sample import loadSampleYml
    sample = loadSampleYml(sampleyml)
    # the sample kernel need information of E and hkl
    Q, hkl2Qmat, psi = calcQ(sampleyml, Ei, E, hkl, psi_axis, Npsisegments=10)
    log.write( "Computed:\n" )
    log.write( "* psi=%s degree\n" % (psi/np.pi*180,) )
    log.write( "* Q=%s\n" % (Q,) )
    log.write( "* hkl2Qmat=%s\n" % (hkl2Qmat,) )
    kfv, Ef = computeKf(Ei, E, Q, log)
    log.write( "* Ef=%s\n" % (Ef,))
    pixel_position = computePixelPosition(kfv, instrument, log)
    log.write( "* pixel_position=%s\n" % (pixel_position,) )
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
    from mcvine_workflow.singlextal.scaffolding import createSampleAssembly
    sampledir = os.path.abspath(os.path.join(outdir, 'sample'))
    createSampleAssembly(sampledir, sample, add_elastic_line=False)
    # create sim script
    params = dict(
        beam_neutrons_path = os.path.join(beam, 'out', 'neutrons'),
        instr_name =  instrument.name,
        instr_dr = instrument.detsys_radius,
        instr_L_m2s = instrument.L_m2s,
        instr_s2b = instrument.offset_sample2beam,
        p_r = pixel.radius, p_h=pixel.height, p_p = pixel.pressure,
        pixel_position = pixel_position,
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
import mcvine.cli
from numpy import array
from dgsres.singlextal import use_res_comps as urc
beam_neutrons_path = %(beam_neutrons_path)r
instrument = urc.instrument(%(instr_name)r, %(instr_dr)r, %(instr_L_m2s)r, %(instr_s2b)r)
samplexmlpath = %(samplexmlpath)r
psi = %(psi)r
hkl2Q = %(hkl2Q)r
pp = %(pixel_position)r
pixel = urc.pixel(%(p_r)r, %(p_h)r, %(p_p)r, position=(pp[1], pp[2], pp[0]))
t_m2p = %(t_m2p)r
Q = %(Q)r
E = %(E)r
hkl_projection = %(hkl_projection)r
urc.run(
    beam_neutrons_path, instrument, samplexmlpath, psi, hkl2Q, pixel, t_m2p,
    Q, E, hkl_projection, Nbuffer=100000)
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

def computePixelPosition(kfv, instrument, log):
    # ** compute nominal TOF at detector pixel **
    # where is detector pixel?
    # cylinder radius = 3meter. 
    R = mcvine.units.parse(instrument.detsys_radius)/mcvine.units.meter
    t_sample2pixel = R/(kfv[0]**2 + kfv[1]**2)**.5
    pixel_pos = kfv*t_sample2pixel
    log.write( "* pixel positon=%s\n" % (pixel_pos,))
    return pixel_pos

def calcQ(sampleyml, Ei, E, hkl, psi_axis, Npsisegments=10):
    from mcvine_workflow.singlextal.io import loadXtalOriFromSampleYml
    xtalori = loadXtalOriFromSampleYml(sampleyml)
    psimin, psimax, dpsi = psi_axis.min, psi_axis.max, psi_axis.step
    from mcvine_workflow.singlextal.solve_psi import solve
    results = solve(
        xtalori, Ei, hkl, E, psi_min=psimin, psi_max=psimax,
        Nsegments = Npsisegments)
    assert len(results)
    from mcvine_workflow.singlextal.coords_transform import hkl2Q
    r0 = results[0]
    psi_in_degrees = r0
    xtalori.psi = psi = r0*np.pi/180.
    Q = hkl2Q(hkl, xtalori) # here the convension is z vertical
    hkl2Qmat = xtalori.hkl2cartesian_mat()
    return Q, hkl2Qmat, psi


def run(beam_neutrons_path, instrument, samplexmlpath, psi, hkl2Q, pixel, t_m2p,
        Q, E, hkl_projection, Nbuffer=100000):
    from mcni.components.NeutronFromStorage import NeutronFromStorage
    source = NeutronFromStorage('source', path=beam_neutrons_path)
    from mccomponents.sample import samplecomponent
    sample = samplecomponent( 'sample', samplexmlpath)
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
    sim_chain = mcni.instrument( [source, sample, pixel_comp] )
    # put components into place
    geometer = mcni.geometer()
    z_beam = mcvine.units.parse(instrument.offset_sample2beam)/mcvine.units.meter
    geometer.register( source, (0,0,z_beam), (0,0,0) )
    geometer.register( sample, (0,0,0), (0,psi*180/np.pi,0) )
    #  8.69538059e-01,   2.87121987e+00,  -1.66786532e-16
    # (2.8712198737660715, -1.6678653212531173e-16, 0.86953805925373207)
    geometer.register( pixel_comp, pixel.position, (0,0,0) )
    # 
    Q2hkl = np.linalg.inv(hkl2Q)
    # lengh of hkl_projection squared
    hkl_proj_len2 = np.dot(hkl_projection, hkl_projection)
    # neutron buffer
    from mcni.neutron_storage.idf_usenumpy import count
    N0 = count(beam_neutrons_path)
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
            before_detected, after_detected, \
            before_dummy_end, after_dummy_end = tracer._store
        incident = after_incident
        # has to be this. pixel and beam has no relative rotation. 
        # should not used after_scattered.
        # it is in the sample's coordinate system,
        # which is rotated by angle psi.
        scattered = before_detected 
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
        print
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
