# Getting started with the SYNTHETIC package


## Test the installation

Right now the best way to run this package is via a pre set conda environment. 

At the moment this is available in the LMU USM server, and is planned to be added for the NERSC LSST DESC server.

On the USM you should use the following command to activate the environment 
    
    conda activate /home/moon/vargatn/anaconda3/envs/galsim

and use the [checkup script](./env_checkup.py)

    python env_checkup.py

If this runs without errors, then you are set to go!

In case of a manual install, read more about the required dependencies [here](../DEPENDENCIES.md)

## Part A: The generative model

:red_circle: **A1** [Prepare catalog](A1_prepare_catalogs.ipynb)
* status:  TBA

:red_circle: **A2** [Train cluster model](A2_train_cluster_model.ipynb)
* status:  TBA

:yellow_circle: **A3**  [Work with galaxy distributions](A3_work_with_galaxy_distributions.ipynb)
* status:  scope added, needs further figures and descriptions

:red_circle: **A4**  [Draw cluster catalog](A4_draw_cluster_catalog.ipynb)
* status:  TBA

## Part B: Rendering images

:green_circle: **B1**  [Render images from catalog](B1_render_image.ipynb)
* Learn to turn a galaxy catalog into an image
* Turn fits images into color composites

:green_circle: **B2**  [Inject galaxy cluster](B2_inject_image.ipynb)
* Learn to inject a galaxy cluster at the catalog level 
* Apply shear to background sources, based on lens geometry

:green_circle: **B3**  [Add intra-cluster light (ICL) to images](B3_add_ICL.ipynb)
* Learn to render intra cluster light as a synthetic image
* apply intra cluster light to the rendered images

## Part C: Processing with metacalibration

:green_circle: **C1**  [Running sextractor](C1_running_sextractor.ipynb)
* Learn to use the synthetic package to automatically wrap sextractor

:green_circle: **C2**  [create MEDS](C2_creat_MEDS.ipynb)
* Learn to create MEDS files
* Prepare for running metacalibration

:yellow_circle: **C3**  [metacalibration on a grid](C3_metacal_on_a_grid.ipynb)
* run a minimal version of metacalibration 
* recover constant shear from a grid of postage stamp galaxies

:red_circle: **C4**  [metacalibration of cluster injection](C4_metacal_on_cluster_injections.ipynb)
* status:  TBA
