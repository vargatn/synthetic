"""
This script is intended to test the correct setup of the python environment

On the USM this script can be run with the following conda environment

    conda activate /home/moon/vargatn/anaconda3/envs/galsim

"""

def test_basic():
    print('testing basic packages')
    import math
    import os
    import time
    import multiprocessing as mp
    import copy
    import glob
    import healpy
    import pickle
    import subprocess
    import scipy

    import numpy as np
    print('\tnumpy version:', np.__version__)
    import pandas as pd
    print('\tpandas version:', pd.__version__)

    try:
        from collections.abc import Iterable
    except ImportError:
        from collections import Iterable

    import fitsio as fio
    print('\tfitsio version:', fio.__version__)

    print('\tall imported successfully')

def test_sklearn():
    print('testing basic packages')
    import sklearn
    print('\tsklearn version:', sklearn.__version__)

    import sklearn.neighbors as neighbors
    import sklearn.decomposition as decomp
    print('\tall imported successfully')

def test_render():
    print('testing galsim & rendering packages')
    import galsim
    print('\tgalsim version:', galsim.__version__)
    import scipy.special
    import scipy.optimize as optimize

    import astropy.cosmology as cosmology
    import astropy.units as u
    import images


    import ngmix
    print('\tngmix version:', ngmix.__version__)
    from ngmix.gmix import GMixBDF
    print('galaxy writer object imported from ngmix')

    print('\tall imported successfully')


def test_metacal():
    print('testing metacal & ngmix packages')
    import ngmix
    print('\tngmix version:', ngmix.__version__)

    from ngmix.medsreaders import NGMixMEDS

    import astropy
    print('\tastropy version:', astropy.__version__)

    from astropy.table import Table, vstack, hstack
    print('\tall imported successfully')



packages = (
    test_basic,
)

if __name__ == "__main__":
    print('\nOn the USM this script can be run with the following conda environment:')
    print('\tconda activate /home/moon/vargatn/anaconda3/envs/galsim\n')
    print("starting environment checkup")
    test_basic()
    test_sklearn()
    test_render()
    test_metacal()