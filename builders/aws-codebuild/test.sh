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

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
echo "localhost slots=8" > $(dirname $(dirname $(which python)))/etc/openmpi-default-hostfile

py.test
