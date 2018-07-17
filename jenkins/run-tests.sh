#!/usr/bin/env bash

echo "Testing dgsres"

export PATH=$HOME/mc/bin:$PATH
source activate test-dgsres

# check mcvine
conda list mcvine

# check mcvine
python -c "import mcvine"
mcvine

# check mantid
python -c "import mantid"

set -e
conda install pytest
py.test -s
