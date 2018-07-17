#!/usr/bin/env python

import os
from setuptools import setup, find_packages

here = os.path.dirname(__file__)
version_ns = {}
with open(os.path.join(here, 'dgsres', '_version.py')) as f:
    exec(f.read(), {}, version_ns)

# define distribution
setup(
    name = "dgsres",
    version = version_ns['__version__'],
    packages = find_packages(".", exclude=['tests', 'notebooks', 'jenkins']),
    package_dir = {'': "."},
    data_files = [],
    test_suite = 'tests',
    install_requires = [
    ],
    dependency_links = [
    ],
    author = "MCViNE team",
    description = "Resoluation calculators for neutron DGS instruments",
    license = 'BSD',
    keywords = "instrument, neutron",
    url = "https://github.com/mcvine/dgsres",
    # download_url = '',
)

# End of file
