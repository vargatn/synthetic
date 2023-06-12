# Getting started with the SYNTHETIC package



There are two ways to run this package

### Following a manual install

Instructions for this are found [here](../README.md)

You can read more about the required dependencies [here](../DEPENDENCIES.md)

###  Pre set conda environment.  
At the moment this is available in the LMU USM server, where you should use the following command to activate the environment 
    
    conda activate /home/moon/vargatn/anaconda3/envs/galsim


In either case, you can use the [checkup script](./env_checkup.py)

    python env_checkup.py

If this runs without errors, then you are set to go!

**Tutorial Status:**

1) Part A :green_circle: functional and data hosted externally
2) Part C :green_circle: functional and data packaged
3) Part B :green_circle: functional and data packaged


## Part A: The generative model

This section is run first at a place where one can access the DC2 dataset along with the `desc_stack_weekly`

**note** this is quite time intensive, is primarily provided for reference on how the data is prepared 

:large_blue_circle: **A0** [Prepare catalog on NERSC](A0_prepare_catalogs_on_NERSC.ipynb)
* status:  Must be validated again, should seek people who can do it
* data is hosted on NERSC

Then the rest can be ran on your local machine / server as well. In case you want to process the DC2 data locally,
then overview notebook A0 to understand how to prepare for the rest of the examples. Otherwise you can go ahead with the pre-packaged data products for each tutorial

:green_circle: **A1** [Train cluster model](A1_train_cluster_model.ipynb)
* Load and collate galaxy catalogs 
* Set up cluster line-of-sight emulation script
* INTENSIVE Draw samples from the proposal galaxy catalog and calculate survival scores for rejection sampling.
* data hosted separately, see [DATA ACCESS](DATA.md)

:green_circle: **A2**  [Work with galaxy distributions](A2_work_with_galaxy_distributions.ipynb)
* Define KDE model for cluster member galaxy features
* Visualize model against simulations
* data hosted separately, currently on google drive

:green_circle: **A3**  [Draw cluster catalog](A3_draw_cluster_catalog.ipynb)
* Draw random realizations from the KDE galaxy cluster model
* data hosted separately, currently on google drive

## Part B: Rendering images

:green_circle: **B1**  [Render images from catalog](B1_render_image.ipynb)
* Learn to turn a galaxy catalog into an image
* Turn fits images into color composites
* data included in repo
* Validated by G. Queirolo

:green_circle: **B2**  [Inject galaxy cluster](B2_inject_image.ipynb)
* Learn to inject a galaxy cluster at the catalog level 
* Apply shear to background sources, based on lens geometry
* data included in repo
* Validated by G. Queirolo

:green_circle: **B3**  [Add intra-cluster light (ICL) to images](B3_add_ICL.ipynb)
* Learn to render intra cluster light as a synthetic image
* apply intra cluster light to the rendered images
* data included in repo
* Validated by G. Queirolo

## Part C: Processing with metacalibration

:green_circle: **C1**  [Running sextractor](C1_running_sextractor.ipynb)
* Learn to use the synthetic package to automatically wrap sextractor
* data included in repo
* Validated by G. Queirolo

:green_circle: **C2**  [create MEDS](C2_create_MEDS.ipynb)
* Learn to create MEDS files
* Prepare for running metacalibration
* data included in repo
* Validated by G. Queirolo

:green_circle: **C3**  [metacalibration on a grid](C3_metacal_on_a_grid.ipynb)
* run a minimal version of metacalibration 
* recover constant shear from a grid of postage stamp galaxies
* No external data needed, mock galaxy grid created within the example

:green_circle: **C4**  [metacalibration of cluster injection](C4_metacal_on_cluster_injections.ipynb)
* Use the metacalibration algorithm in the cluster lensing scenario
* Measure the tangential shear field induced by a galaxy cluster
* data included in repo

