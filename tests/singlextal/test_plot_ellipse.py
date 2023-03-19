#!/usr/bin/env python
#
# Jiao Lin <jiao.lin@gmail.com>

import os, glob, sys, contextlib

from matplotlib import pyplot as plt
import numpy as np

import dgsres
from dgsres.singlextal import plot as resplot, workflow, violini

thisdir = os.path.abspath(os.path.dirname(__file__))

def test_plot_ellipse():
    configdir = os.path.join(thisdir, '../../notebooks/singlextal/Mn3Si2Te6-SEQ')
    mcvinesim = os.path.join(thisdir, '../data/SEQUOIA_data/Mn3Si2Te6/mcvine-res-sim')
    workdir = os.path.abspath("work.plotellipse")
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    copy_config_files(configdir, workdir)
    with chdir(workdir):
        import imp
        config = imp.load_source('config', 'convolution_config.py')
        sl = config.rwc.slices[0]
        print(sl)
        class McvineResData(resplot.McvineResolutionData):
            def path(self, q, E):
                return os.path.join(self.parent_dir, config.rwc.simdir(q, E, sl))
        mrd = McvineResData(mcvinesim, dirname_template=None)
        q1, E1 = 3.1, 19.
        hkl1 = sl.hkl0+sl.hkl_projection*q1
        mcvine_pc1 = mrd.loadPointCloud(q1, E1)
        mcvine_cm1 = mrd.computeCovMat(q1, E1, dhkl_ranges=[(-0.2,0.2)]*3)
        print(mcvine_cm1)
        violini_model = workflow.create_violini_model(config.rwc, 0.01)
        violini_cm1 = violini_model.computeCovMat(hkl1, E1)
        print(violini_cm1)
        violini_pc1 = violini_model.computePointCloud(hkl1, E1)
        plot_qE_with_violini(mcvine_pc1, violini_cm1, 'res_qE_with_violini.png')

def plot_qE_with_violini(mcvine_pc1, violini_cm1, outfile):
    plt.figure(figsize=(10,4))
    directions = [
        ('h', [1,0,0]),
        ('k', [0,1,0]),
        ('l', [0,0,1]), 
    ]
    for i, direction in enumerate(directions):
        name, vector = direction
        axis1 = name, -0.5, 0.5, 0.002
        axis2 = 'E', -5, 5., 0.02
        otheraxes = [ (n, -0.015, 0.015) for n,v in directions if n !=name]

        plt.subplot(1,3,i+1)
        hg1, Eg1, I1 = mcvine_pc1.getThinSlice(axis1, axis2, *otheraxes)
        plt.pcolormesh(hg1, Eg1, I1.T)
        plt.ylim(-1.5,1.5)
        plt.xlim(-.125, .125)
        plt.xlabel(name)
        plt.ylabel('E (meV)')
        plt.clim(0, np.max(I1)/2)
        resplot.plot_qE_ellipse(violini_cm1, vector, 'r') #, label='Violini')
        # plt.legend()
    plt.tight_layout()
    plt.savefig(outfile)
    print(f"made plot {outfile}")
    plt.close()
    return

def copy_config_files(src, dest):
    files = 'sample.yaml convolution_config.py resolution_workflow_config.py '.split()
    import shutil
    for f in files:
        shutil.copyfile(os.path.join(src, f), os.path.join(dest, f))
    return

@contextlib.contextmanager
def chdir(path):
    """Sets the cwd within the context
    """
    origin = os.path.abspath(os.curdir)
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)

def main():
    test_plot_ellipse()

if __name__ == '__main__': main()

# End of file
