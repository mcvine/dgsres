<!-- [![Build Status](http://35.168.96.122:8080/buildStatus/icon?job=dgsres)](http://35.168.96.122:8080/job/dgsres/) -->
[![Build Status](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoibUZKS0FlY080d1I3ZVFyM2ZMVm8reFJGemNYcTFBOE9mVG1rWkVVcnVYOStmaHNIOGUvb1piL2FTT2NWbGJVRXlOVkFKK2N3OXQ2ZzlGOXZTL0MrZ3pRPSIsIml2UGFyYW1ldGVyU3BlYyI6Iit6YVNCQXZHL3lUUVUzdXkiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=master)](https://console.aws.amazon.com/codesuite/codebuild/668650830132/projects/mcvine-dgsres-py2)
[![DOI](https://zenodo.org/badge/97215709.svg)](https://zenodo.org/badge/latestdoi/97215709)

# dgsres
Neutron Direct Geometry Spectrometer Resolution Calculators

This repository is a collection of resolution calculators for direct-geometry inelastic neutron spectrometers.

Resolution for direct geometry spectrometers is asymmetric and non-stationary. The following is the mcvine simulation and model fitting results for an ARCS experiment:
![image](https://user-images.githubusercontent.com/1796155/59158473-a8947e00-8a88-11e9-9e4c-0158ee5e0443.png)


## Dependencies

* tqdm
* latex (texlive in ubuntu)
* pylatex
* lmfit
* cloudpickle
* scikit-image 0.15
* mcpl (optional. install use pip)
* dill (obsolete. now use cloudpickle)

## Installation
`$ python setup.py install`


## Related projects
* [Online application](http://rez.mcvine.ornl.gov) and its [source code](https://github.com/sns-chops/resolution/tree/master/dashui)
* [Experimental data and PyChop modeling](https://github.com/sns-chops/resolution)
