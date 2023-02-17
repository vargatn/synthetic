#!/bin/bash
# this script must be run in interactive mode
# bash -i install.sh

conda create -y -n synth python=3.9 anaconda

conda activate synth

conda env list

conda install -y -c conda-forge ngmix=1.3.9 # this is a really complex package, if this installs, probably everything else works
conda install -y -c conda-forge galsim
conda install -y -c conda-forge fitsio
conda install -y -c conda-forge esutil
conda install -y -c conda-forge healpy

pip install git+https://github.com/esheldon/images.git \
            git+https://github.com/esheldon/meds \
            git+https://github.com/esheldon/psfex.git \

python ./tutorial/env_checkup.py

