import copy
import multiprocessing as mp

import numpy as np
import pandas as pd

from ..tools import partition


class KFoldCV(object):
    def __init__(self, cont, nfold=5, nprocess=100, seed=None):
        """
        K-fold crossvalidator for KDE containers
        Parameters
        ----------
        cont: KDEContainer
            KDE container for the data which will be cross validated
        nfold: int
            number of folds
        nprocess: int
            number of CPUs to use for scoring
        seed: int
            random seed for numpy random
        """
        self.cont = cont
        self.cont.standardize_data()
        self.nfold = nfold
        self.nprocess = nprocess
        self.seed = seed
        self.rng = np.random.RandomState(seed)

        self.labels = self.rng.randint(0, self.nfold, len(self.cont.data))

    def loop_bandwidths(self, bandwidths):
        """
        Loop over each bandwidth and evaluate K-fold scores
        Parameters
        -----------
        bandwidths: list
            bandwidths to evaluate
        Returns
        -------
            list of scores for each bandwidth (each element is a list of K-scores for each fold)
        """
        scores = []
        for bw in bandwidths:
            print(bw)
            mscores = self._loop_cv(bw)
            scores.append(mscores)
        return scores

    def _loop_cv(self, bandwidth):

        mscores = []
        for ifold in np.arange(self.nfold):
            print("\t", ifold)
            train, test = self._split(ifold)

            iarr = np.arange(len(test.data))
            iparts = partition(iarr, self.nprocess)

            infodicts = []
            for i in np.arange(len(iparts)):
                info = {
                    "train": train,
                    "test_data": test.data.iloc[iparts[0]],
                    "test_weights": test.weights.iloc[iparts[0]],
                    "bandwidth": bandwidth,
                }
                infodicts.append(info)

            result = run_cv_scores(infodicts)
            mscore = np.sum(result["weights"] * (result["scores"] + np.log(result["jac"]))) / np.sum(result["weights"])
            mscores.append(mscore)

        return mscores

    def _split(self, label):
        """splits data into train and test k-fold"""

        _ind = self.labels != label
        train = self._shrink(_ind)

        _ind = self.labels == label
        test = self._shrink(_ind)
        return train, test

    def _shrink(self, index):
        _cont = copy.copy(self.cont)
        _cont.data = _cont.data.iloc[index].copy()
        _cont.weights = _cont.weights.iloc[index].copy()

        _cont._data = _cont.pca_transform(_cont.data)
        return _cont


def run_cv_scores(infodicts):
    pool = mp.Pool(processes=len(infodicts))
    try:
        pp = pool.map_async(_score_cv_samples, infodicts)
        # the results here should be a list of score values
        result = pp.get(86400)  # apparently this counters a bug in the exception passing in python.subprocess...
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        pool.join()
    else:
        pool.close()
        pool.join()

    return pd.concat(result)


def _score_cv_samples(info):
    train = info["train"]

    train.construct_kde(bandwidth=info["bandwidth"])

    scores, jac = train.score_samples(info["test_data"])

    result = pd.DataFrame()
    result["scores"] = scores
    result["jac"] = np.ones(len(result)) * jac
    result["weights"] = info["test_weights"].values

    return