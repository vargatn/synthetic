for the newest conda version run 
    
    conda init bash

to recover the bash command prompt, then create new environment named `inkind`

    conda create -n inkind python=3

then activate the environment by typing 

    conda activate inkind

Packages which should be installed

    conda install -c conda-forge numpy scipy numba pandas matplotlib scikit-learn astropy  fitsio esutil

note the first import of matplotlib will take a longer time for compiling
Now install galsim

    conda install -c conda-forge galsim

install `meds` 

    git clone git@github.com:esheldon/meds.git
    cd meds
    python setup.py install

and `psfex`,  (this is just a python wrapper)

    git clone https://github.com/esheldon/psfex
    cd psfex
    python setup.py install

intall ngmix 1.3.9 from the repo

    git clone git@github.com:esheldon/ngmix.git
    cd ngmix
    git checkout v1.3.9
    python setup.py install

Install sextractor and psfex

     conda install -c conda-forge astromatic-source-extractor
    conda install -c conda-forge psfex

