#!/usr/bin/env bash

echo "Testing dgsres"

export PATH=$HOME/mc/bin:$PATH
conda remove -n test-dgsres --all
conda create -n test-dgsres
source activate test-dgsres

set -e
conda config --add channels mantid  # need mantid-framework
conda install pytest
conda install -c mcvine/label/unstable mcvine
# conda install mantid-framework
conda list mcvine

export AWS_S3_PROFILE_NAME=ndav_mcvine
py.test -s
