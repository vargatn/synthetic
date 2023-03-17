# How to access data for the part A of the tutorial

The training data extracted from the LSST DC2 simulations is too big to be included in the git repository

At the moment the key files are accessed from the following google drive links:

[clust_dc2-sim-LOS_v1.h5](https://drive.google.com/file/d/1Zzm7PJm96xvtxfb1xqkip4w9SjVPV84J/view?usp=share_link)

[cosmoDC2_v1.1.4_refpixels.h5](https://drive.google.com/file/d/1v6k1V4kaDFfLE-qgWQbQ54x2D3BYCA5J/view?usp=share_link)

[dc2-alpha_concentric_sample-v01_test-03.tar.gz](https://drive.google.com/file/d/1Rwe-ZSkIvUMbTgGDiP3HLP0nC7fiQIhv/view?usp=sharing)

Alternatively the same files can be produced via notebooks A0, A1 and A2

If you would prefer a different way of accessing the data files, please contact us at t.varga@physik.lmu.de

## Setting up the file structure

The data files for this example calculation are pre packaged, and should be downloaded from a link provided upon request.

    1 dc2-alpha_concentric_sample-v01_test-03.tar.gz
    2 cosmoDC2_v1.1.4_refpixels.h5
    3 clust_dc2-sim-LOS_v1.h5
    
These should be downloaded and placed in a file structure such that

    /root/
    |----/resamples/ 
    |----/dc2-alpha_concentric_sample-v01_test-03.tar.gz
    |----/dc2_cluster_sim_cutouts/cosmoDC2_v1.1.4_refpixels.h5
    |----/dc2_cluster_sim_cutouts/clust_dc2-sim-LOS_v1.h5
    
from within the root folder, extract the .tar.gz file using the command

    tar xzf dc2-alpha_concentric_sample-v01_test-03.tar.gz -C  resamples --strip-components 1    
    
This should yield a file structure as below
 
     /root/
    |----/resamples/ 
    |--------------/dc2-alpha_concentric_sample-v01_test-03_run0_1846435878_rbin0.p
    |--------------/dc2-alpha_concentric_sample-v01_test-03_run0_1846435878_rbin0_samples.fits
    |--------------/dc2-alpha_concentric_sample-v01_test-03_run0_1846435878_rbin0_scores.fits
            .
            .
            .
    |--------------/dc2-alpha_concentric_sample-v01_test-03_run3_664487101_rbin3.p
    |--------------/dc2-alpha_concentric_sample-v01_test-03_run3_664487101_rbin3_samples.fits
    |--------------/dc2-alpha_concentric_sample-v01_test-03_run3_664487101_rbin3_scores.fits            
    |----/dc2-alpha_concentric_sample-v01_test-03.tar.gz
    |----/dc2_cluster_sim_cutouts/cosmoDC2_v1.1.4_refpixels.h5
    |----/dc2_cluster_sim_cutouts/clust_dc2-sim-LOS_v1.h5
