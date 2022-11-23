import copy

import numpy as np
import pandas as pd
import sklearn.neighbors as neighbors
import sklearn.decomposition as decomp


def weighted_mean(values, weights):
    average = np.average(values, axis=0, weights=weights)
    return average


def weighted_std(values, weights):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, axis=0, weights=weights)
    variance = np.average((values - average) ** 2, axis=0, weights=weights)
    return np.sqrt(variance)


class KDEContainer(object):
    """Represent dataset as a KDE in it's eigne-frame (PCA)"""
    _default_subset_sizes = (2000, 5000, 10000)
    _kernel = "gaussian"
    # _kernel = "tophat"
    _atol = 1e-6
    _rtol = 1e-6
    _breadth_first = False
    _jacobian_matrix = None
    _jacobian_det = None

    def __init__(self, raw_data, weights=None, transform_params=None, seed=None):
        """
        Represent dataset as a KDE in it's eigne-frame (PCA)
        Parameters
        ----------
        raw_data: pd.DataFrame
            input data to represent
        weights: np.array
            weight for each element
        transform_params: bool
            deprecated
        seed: int
            np random seed
        """
        if seed is not None:
            self.rng = np.random.RandomState(seed)
        else:
            self.rng = np.random.RandomState()

        self.data = raw_data
        self.columns = raw_data.columns
        if weights is None:
            self.weights = np.ones(len(raw_data), dtype=float)
        else:
            self.weights = weights.astype(float)
        self.ndim = self.data.shape[1]

    def set_seed(self, seed):
        self.rng.seed(seed)

    @staticmethod
    def _weight_multiplicator(arr, weights):
        """
        turns element wise-weight into multiplicitiy of elements in list
        Parameters
        ----------
        arr: np.array
            array with weights
        weights: np.array
            weights
        Returns
        -------

        """
        multiplier = np.round(weights)
        newarr = []
        for i in np.arange(len(arr)):
            for j in np.arange(multiplier[i]):
                newarr.append(arr[i])
        newarr = np.vstack(newarr)
        return newarr

    def shuffle(self):
        self.sample(n=None, frac=1.)

    def sample(self, n=None, frac=1.):
        inds = np.arange(len(self.data))
        print(self.data.shape)
        print(inds.shape)
        tab = pd.DataFrame()
        tab["IND"] = inds
        if n and n > len(tab):
            n = None
        inds = tab.sample(n=n, frac=frac)["IND"].values

        self.data = self.data.iloc[inds].copy().reset_index(drop=True)
        self.weights = self.weights.iloc[inds].copy().reset_index(drop=True)

    def fit_pca(self):
        """Standardize -> PCA -> Standardize"""
        # _data = self.data
        self.mean1 = weighted_mean(self.data, self.weights)
        _data = self.data - self.mean1
        #
        # Add here a PCA weights pre-burner,
        # draw a subset of 100k rows, then multiplicate them according to weights
        # fit the PCA on theose new rows
        subset = self.select_subset(_data, self.weights, nsample=100000)
        # subset = self._weight_multiplicator(subset.values, self.weights)
        self.pca = decomp.PCA()
        self.pca.fit(subset)

        _data = self.pca.transform(_data)
        self.std2 = weighted_std(_data, self.weights)

        rotation_matrix = self.pca.components_
        scale_matrix = np.diag(1. / self.std2)

        # this is the forward transformation from raw data to processed data
        self._jacobian_matrix = np.dot(scale_matrix, rotation_matrix)
        # self._jacobian_matrix = rotation_matrix

        # this is the inverse transformation from processed data to raw data
        self._jacobian_matrix_inv = np.linalg.inv(self._jacobian_matrix)
        # in the KDE we need the Jacobi determinat of the inverse transformation
        self._jacobian_det = np.linalg.det(self._jacobian_matrix_inv)
        # self._jacobian_det = 1.

        self.pca_params = {
            "mean1": self.mean1.copy(),
            "std2": self.std2.copy(),
            "pca": copy.deepcopy(self.pca),
        }

    def pca_transform(self, data):
        # _data = data
        _data = data - self.mean1
        _data = self.pca.transform(_data)
        _data /= self.std2
        return _data

    def pca_inverse_transform(self, data):
        # _data = data
        _data = data * self.std2
        _data = self.pca.inverse_transform(_data)
        _data = _data + self.mean1
        res = pd.DataFrame(_data, columns=self.columns)
        return res

    def standardize_data(self):
        self.fit_pca()
        self._data = self.pca_transform(self.data)
        # self._data = self.data

    def select_subset(self, data, weights, nsample=10000):
        # if nsample > len(data):
        # nsample = len(data)
        indexes = np.arange(len(data))
        ww = weights / weights.sum()
        inds = self.rng.choice(indexes, size=int(nsample), p=ww, replace=True)
        subset = data.iloc[inds]
        return subset

    def construct_kde(self, bandwidth):
        """"""
        self.bandwidth = bandwidth
        self.kde = neighbors.KernelDensity(bandwidth=self.bandwidth, kernel=self._kernel,
                                           atol=self._atol, rtol=self._rtol, breadth_first=self._breadth_first)
        self.kde.fit(self._data, sample_weight=self.weights)

    def random_draw(self, num, rmin=None, rmax=None, rcol="LOGR"):
        """draws random samples from KDE maximum radius"""
        _res = self.kde.sample(n_samples=int(num), random_state=self.rng)
        self.res = self.pca_inverse_transform(_res)
        if (rmin is not None) or (rmax is not None):
            if rmin is None:
                rmin = self.data[rcol].min()

            if rmax is None:
                rmax = self.data[rcol].max()

            # these are the indexes to replace, not the ones to keep...
            inds = (self.res[rcol] > rmax) | (self.res[rcol] < rmin)
            while inds.sum():
                vals = self.kde.sample(n_samples=int(inds.sum()), random_state=self.rng)
                _res[inds, :] = vals
                self.res = self.pca_inverse_transform(_res)
                inds = (self.res[rcol] > rmax) | (self.res[rcol] < rmin)

        self.res = self.pca_inverse_transform(_res)
        return self.res

    def score_samples(self, arr):
        """Assuming that arr is in the data format"""

        arr = self.pca_transform(arr)
        res = self.kde.score_samples(arr)
        return res, self._jacobian_det

    def drop_kde(self):
        self.pca = None
        self.kde = None

    def drop_col(self, colname):
        self.data = self.data.drop(columns=colname)
        self.columns = self.data.columns
        self.ndim = len(self.columns)