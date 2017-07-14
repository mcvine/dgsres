# -*- Python -*-
#
# Jiao Lin <jiao.lin@gmail.com>
#

import os, numpy as np, mcvine


def computeEi_and_t0(beampath, instrument):
    """use Ei saved in props.json. use saved neutrons information to compute t0
    """
    Ei = getEi(beampath)
    tof = getTof(beampath)
    from mcni.utils import conversion as Conv
    vi = Conv.e2v(Ei)
    L_m2s = mcvine.units.parse(instrument.L_m2s)/mcvine.units.meter
    t0 = tof - L_m2s/vi
    # print Ei, t0*1e6
    return Ei, t0*1e6


def getEi(beampath):
    import json, os
    beam_outdir = os.path.join(beampath, 'out')
    props_path = os.path.join(beam_outdir, 'props.json')
    s = open(props_path).read()
    s = s.replace("'", '"') # ' -> "
    props = json.loads(s)
    Ei, unit = props['average energy'].split()
    return float(Ei)


def getTof(beampath):
    import json, os
    beam_outdir = os.path.join(beampath, 'out')
    tof_hist_path = os.path.join(beam_outdir, 'itof.h5')
    import histogram.hdf as hh
    tofhist = hh.load(tof_hist_path)
    tofcenter = (tofhist.tof*tofhist.I).sum()/tofhist.I.sum()
    return tofcenter


def computeEi_and_t0_usingMantid_ARCS(beampath, instrument='ARCS'):
    """use mantid to compute Ei and t0
    """
    # only implementation right now
    assert instrument=='ARCS'
    # create dummy nxs
    import hashlib
    key = hashlib.sha224(beampath).hexdigest()
    outdir = 'computeEi_and_t0-%s' % key
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    dummynxs = os.path.join(outdir, 'dummy.nxs')
    if not os.path.exists(dummynxs):
        create_dummy_nxs(dummynxs, beampath)
    results = os.path.join(outdir, 'results.txt')
    if not os.path.exists(results):
        # call mantid
        from mantid import simpleapi as mtdsa
        ws = mtdsa.Load(dummynxs, LoadMonitors=True)
        Ei, firstMonitorPeak, FirstMonitorIndex, t0 = mtdsa.GetEi(ws[1])
        open(results, 'wt').write('%s,%s' % (Ei,t0))
    else:
        Ei,t0 = eval(open(results).read().strip())
    return Ei, t0

def create_dummy_nxs(out, beam):
    import shutil, sys
    from mcvine.instruments.ARCS.nxs.raw import nxs_template, populateEiData
    shutil.copyfile(nxs_template, out)
    import time; time.sleep(0.5)
    import h5py
    f = h5py.File(out, 'a')
    entry = f['entry']
    populateEiData(entry, os.path.join(beam, 'out'))
    return out

def test():
    beam = "/SNS/users/lj7/simulations/ARCS/beam/100meV-n1e10"
    Ei,t0 = computeEi_and_t0(beam)
    assert np.isclose(Ei, 100.482385711)
    assert np.isclose(t0, 18.7623408857)
    return

if __name__ == '__main__': test()

# End of file 
