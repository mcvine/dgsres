#!/usr/bin/env bash

echo "Install"

export PATH=$HOME/mc/bin:$PATH
conda remove -n test-dgsres --all
conda create -n test-dgsres
source activate test-dgsres

conda config --add channels conda-forge
conda config --add channels diffpy
conda config --add channels mantid  # need mantid-framework
# conda install mcvine-core
conda install -c mcvine/label/unstable mcvine
