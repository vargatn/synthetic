# Getting started with the SYNTHETIC package

I. Introduction

A. Overview of synthetic package

B. Purpose of the getting started tutorial

Common Use Cases

A. Explanation of common use cases for synthetic data

B. Code examples for each use case


## Test the installation

Right now the best way to run this package is via a pre set conda environment. 

At the moment this is available in the LMU USM server, and is planned to be added for the NERSC LSST DESC server.

On the USM you should use the following command to activate the environment 
    
    conda activate /home/moon/vargatn/anaconda3/envs/galsim

and use the [checkup script](./env_checkup.py)

    python env_checkup.py

If this runs without errors, then you are set to go!

*TBA* hands on dependency install instructions..

## Part A: The generative model

**A1** [Prepare catalog](A1_prepare_catalogs.ipynb)

**A2** [Train cluster model](A2_train_cluster_model.ipynb)

**A3** [Work with galaxy distributions](A3_work_with_galaxy_distributions.ipynb)
status: :large_orange_diamond: scope added

**A4** [Draw cluster catalog](A4_draw_cluster_catalog.ipynb)

## Part B: Rendering images

B1 [Render images from catalog](B1_render_image.ipynb)

B2 [Inject images to catalog](B2_inject_image.ipynb)

B3 [Add intra-cluster light (ICL) to images](B3_add_icl.ipynb)


## Processing with metacalibration


