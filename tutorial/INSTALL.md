# Manual installation guide

The scope of this activity is to set up a conda environment in which
the synthetic package and its dependencies are installed properly.

Please consider that some of the dependencies are really fragile, and follow the instructions carefully

## anaconda setup

Be sure to have anaconda setup on your machine

    conda 23.1.0
    python-3.9

Now create a new environment

    conda create -n synth python=3.9 anaconda

This environment can be activated via 

    conda activate synth

    conda install -c conda-forge ngmix=1.3.9
    conda install -c conda-forge fitsio
    conda install -c conda-forge galsim
    conda install -c conda-forge healpy
    conda install -c conda-forge esutil

    pip install git+https://github.com/esheldon/images.git
    pip install git+https://github.com/esheldon/meds
    pip install git+https://github.com/esheldon/psfex.git
