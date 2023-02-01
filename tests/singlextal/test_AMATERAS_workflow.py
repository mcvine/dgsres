#!/usr/bin/env python

import os, numpy as np, tempfile, shutil, pytest
thisdir = os.path.abspath(os.path.dirname(__file__))

def test():
    workdir = tempfile.mkdtemp(dir=thisdir)
    print(workdir)
    # link beam dir
    beamdir = os.path.join(workdir, 'beam')
    # beamdir_orig = os.path.expanduser('~/beam/AMATERAS/2.63meV')
    beamdir_orig = os.path.join(thisdir, '..', 'data', 'beam', 'AMATERAS', '2.63meV')
    os.symlink(beamdir_orig, beamdir)
    # copy sample yaml file
    shutil.copyfile(
        os.path.join(thisdir, 'amateras_workflow_data/sample.yaml'),
        os.path.join(workdir, 'sample.yaml'),
    )
    import imp
    # copy config file. it should be alongside with beam dir
    config_src = os.path.join(
        thisdir, 'amateras_workflow_data/resolution_workflow_config.py'
    )
    config_dest = os.path.join(workdir, 'confgi.py')
    shutil.copyfile(config_src, config_dest)
    # load config
    config = imp.load_source('config', config_dest)
    # run
    from dgsres.singlextal import workflow
    os.chdir(workdir)
    outputs, failed = workflow.simulate_all_in_one(config)
    f0 = failed[0]
    for item, tb in f0.items():
        print(item)
        print(tb)
    assert len(failed[0]) == 12 and len(outputs[0]) == 8
    nofits = workflow.fit_all_in_one(config)
    # print(outputs)
    # print(failed)
    # print(len(failed[0]), len(outputs[0]), len(nofits[0]))
    assert len(nofits[0]) == 12
    assert sorted(nofits[0]) == sorted(list(failed[0].keys()))
    return

if __name__ == '__main__': test()
