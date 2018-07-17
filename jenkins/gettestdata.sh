#!/bin/bash

echo "Get test data"

export PATH=$HOME/mc/bin:$PATH
source activate testenv

export AWS_S3_PROFILE_NAME=ndav_mcvine
aws s3 sync s3://ndav-mcvine/dgsres/ $HOME/test-data/mcvine-dgsres --profile ${AWS_S3_PROFILE_NAME}
rm -f tests/data
ln -s $HOME/test-data/mcvine-dgsres tests/data
