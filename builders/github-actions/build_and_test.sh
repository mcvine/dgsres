#!/bin/bash

set -x
set -e

# checking installation
export LC_ALL=C.UTF-8
export LANG=C.UTF-8
which py.test
python -c "import mantid"
python -c "import mcvine"
mcvine

python setup.py install

py.test
python tests/singlextal/test_workflow_fit.py
