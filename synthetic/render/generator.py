"""
Mockup generator, create a plausible line-of-sight catalog on the fly for testing purposes
"""

import numpy as np
import pandas as pd


class MockupSurvey(object):
    _names = None
    def __init__(self, maglim=24, seed=None, pos_limits=("rect", ((), ()))):
        """
        Mockup generator, create a plausible line-of-sight catalog on the fly for testing purposes

        quantities to mockup:

            ra, dec, flux_[griz], mag_[griz], z, size, bulge_disc_fraction

        pos_limits:

            ("rect", ((-10, 10), (-10, 10))) type, x_extent, y_extent
            ("circ", ((10.), (0., 0))) type, radius, xy_extent
        """
        self.maglim = maglim
        self.rng = np.random.RandomState(seed)





    def draw(self, num=1e4):
        num = float(num)


