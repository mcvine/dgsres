#!/usr/bin/env python

import os, sys

id = os.environ['AWS_S3_ACCESS_ID']
secret = os.environ['AWS_S3_ACCESS_SECRET']
profile_name = os.environ['AWS_S3_PROFILE_NAME']

CONFIG="""[profile %(profile_name)s]
output = text
region = us-east-1
"""

CRED_template="""[%(profile_name)s]
aws_access_key_id = %(id)s
aws_secret_access_key = %(secret)s
"""

dir = os.path.expanduser("~/.aws")
if not os.path.exists(dir):
    os.makedirs(dir)

config_path = os.path.join(dir, "config")
if os.path.exists(config_path):
    raise IOError("%s already exists" % config_path)
with open(config_path, 'wt') as stream:
    stream.write(CONFIG % locals())


cred_path = os.path.join(dir, "credentials")
if os.path.exists(cred_path):
    raise IOError("%s already exists" % cred_path)
with open(cred_path, 'wt') as stream:
    stream.write(CRED_template % locals())
