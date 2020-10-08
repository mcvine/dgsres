#!/bin/bash

set -x
set -e
export PATH=$HOME/mc/bin:$PATH
source activate test

# get beam data
THIS_SCRIPT_DIR=`dirname $0`
$THIS_SCRIPT_DIR/getdata.sh

# checking installation
export LC_ALL=C.UTF-8
export LANG=C.UTF-8
which py.test
python -c "import mcvine"
mcvine

python setup.py install
py.test 
