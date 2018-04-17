#!/usr/bin/env bash

echo "Testing dgsres"

export PATH=$HOME/mc/bin:$PATH
source activate test-dgsres

set -e

# check mcvine
conda list mcvine

# check mantid
python -c "import mantid"

conda install pytest
py.test -s
