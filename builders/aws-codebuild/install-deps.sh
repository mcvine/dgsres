#!/bin/bash

set -x
set -e
export PATH=$HOME/mc/bin:$PATH
conda create -n test tbb=2020.3 mcvine=1.4 tqdm pytest awscli cloudpickle lmfit pylatex scikit-image=0.15 python=$PYTHON_VERSION
source activate test
conda list mcvine
conda list mcvine-core
THIS_SCRIPT_DIR=`dirname $0`
python $THIS_SCRIPT_DIR/init_mantid_user_config.py
python -c "from mantid import simpleapi"
