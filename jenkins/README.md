# scripts to run at jenkins slaves to build this repo

## Install prerequisites
Modify and run https://github.com/mcvine/systemtests-instruments/blob/master/jenkins/install.sh
to isntall required packages on any jenkins slave for this project.
This should only need to be done once.

## Build
Then in a Jenkins build project configuration, the "Build" panel, the "Execute shell" sub-panel
the "Command" is:

  `$ ./jenkins/build.sh`

The build uses jenkins env vars:
* JOB_NAME
* BUILD_NUMBER