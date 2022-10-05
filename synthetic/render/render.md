Galsim backend
==============
to render a catalog line of sight into an observed image canvas


Preparing to run
-----------------

Start by installing galsim https://github.com/GalSim-developers/GalSim

It's fine to do this install in a separate conda environment.
Let's make that as 

    conda create -n galsim python=3

which you should activate by `conda activate galsim`

Intall the actual galsim package by

    conda install -c conda-forge galsim

Install meds
https://github.com/esheldon/meds
needs also esutil to function 
conda install -c conda-forge esutil

TODO expand:
Following this install ngmix and ngmixer

Default Use case

Load catalog for Field (This should assume a not yet cutout catalog, must include a step which filters for objects which are in the image)

add a step which injects objects only to the overlap area between the stamp and the canvas

Render field

Load catalog for cluster

Render cluster

calculate ICL model for cluster

Overlay ICL model

Save image for cluster

Save location of bounding boxes for objects

Run a simple galsim shear estimator

Run a sextractor on image??? let's hope this works out

Install sextractor as  conda install -c conda-forge astromatic-source-extractor 

test the correct path by using "which sex" command, should point to the conda environment



