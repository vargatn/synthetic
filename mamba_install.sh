
conda create -y --name synth python=3.9
conda activate synth

conda env list

conda install mamba -y -c conda-forge
mamba env update --file base_env.yml --prune

mamba install -y pip

mamba install -y -c conda-forge ngmix=1.3.9
mamba install -y -c conda-forge galsim
mamba install -y -c conda-forge fitsio
mamba install -y -c conda-forge esutil
mamba install -y -c conda-forge healpy
mamba install -y -c conda-forge scikit-learn
mamba install -y -c conda-forge astromatic-source-extractor
mamba install pytables -y -c conda-forge

pip install git+https://github.com/esheldon/images.git \
            git+https://github.com/esheldon/meds \
            git+https://github.com/esheldon/psfex.git \

python setup.py install

python ./tutorial/env_checkup.py
