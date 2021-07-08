#!/bin/bash

set -e

echo "Get beam"

export PATH=$HOME/mc/bin:$PATH
source activate test

THIS_SCRIPT_DIR=`dirname $0`
python $THIS_SCRIPT_DIR/setup-aws-testconfig.py

# export AWS_S3_PROFILE_NAME=ndav_mcvine
aws s3 sync s3://sts-mcvine/beam/ $HOME/beam/ --profile ${AWS_S3_PROFILE_NAME}
aws s3 sync s3://sts-mcvine/dgsres/ $HOME/test-data/mcvine-dgsres --profile ${AWS_S3_PROFILE_NAME}
rm -f tests/data
ln -s $HOME/test-data/mcvine-dgsres tests/data
