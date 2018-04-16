#!/bin/bash

set -e

echo "Get test data"

export PATH=$HOME/mc/bin:$PATH
source activate testenv

export AWS_S3_PROFILE_NAME=ndav_mcvine
aws s3 sync s3://ndav-mcvine/dgsres/ test/data/ --profile ${AWS_S3_PROFILE_NAME}
