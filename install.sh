#!/bin/bash
# this script must be run in interactive mode
# bash -i install.sh

#conda create -y -n synth python=3.9
conda env create -f base_env.yml

conda activate synth

conda env list

conda install pip

conda install -y -c conda-forge ngmix=1.3.9
conda install -y -c conda-forge galsim
conda install -y -c conda-forge fitsio
conda install -y -c conda-forge esutil
conda install -y -c conda-forge healpy
conda install -y -c conda-forge scikit-learn

pip install git+https://github.com/esheldon/images.git \
            git+https://github.com/esheldon/meds \
            git+https://github.com/esheldon/psfex.git \

python setup.py install

python ./tutorial/env_checkup.py

