"""
Emulator.py
"""

import fitsio as fio
import numpy as np
import pandas as pd
import copy
import glob

from .kde import KDEContainer

try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable

BADVAL = -9999

ENDIANS = {
    "little": "<",
    "big": ">",
}

import matplotlib as mpl
try:
    import matplotlib.pyplot as plt
except:
    mpl.use("Agg")
    import matplotlib.pyplot as plt

import multiprocessing as mp

from ..tools import partition


def get_angle(num, rng):
    """draw random angles in radians"""
    angle = rng.uniform(0, np.pi, size=num)
    return angle


class BaseContainer(object):
    def __init__(self):
        """
        Base class for Feature Space Container

        Main goal is to containe

            * alldata
            * features extracted / consturcted from alldata
            * weights for each entry
        """
        self.alldata = None
        self.features = None
        self.weights = None

    def construct_features(self, columns, limits=None, logs=None, **kwargs):
        """
        Create features from raw data, this includes complicated operations as well

        Examples
        --------
            "columns": [
                ("GABS", ("ellipticity_1_true", "ellipticity_2_true", "SQSUM")),
                ("SIZE", "size_true"),
                ("MAG_I", "mag_i"),
                ("COLOR_G_R", ("mag_g", "mag_r", "-")),
                ("COLOR_R_I", ("mag_r", "mag_i", "-")),
                ("COLOR_I_Z", ("mag_i", "mag_z", "-")),
                ("STELLAR_MASS", "stellar_mass"),
                ("HALO_MASS", "halo_mass")
            ],
            "logs": [False, True, False, False, False, False, True, True],
            "limits": [(0., 1.), (-1, 5), (17, 25), (-1, 3), (-1, 3), (-1, 3), (10**3, 10**13), (10**9, 10**16)],

        Thus we can specify that a desired feature is the squared sum of two columns in the alldata

        Supported operations are + - * / SQSUM


        Parameters
        ----------
        columns: tuple
            instructions about the feature construction, see example above
        limits: tuple
            value limits after the feature is constructed
        logs: list
            bool list if the feature should be taken the log of or not
        kwargs
        """
        self.columns = columns
        self.limits = limits
        self.logs = logs
        self.features = pd.DataFrame()

        self.inds = np.ones(len(self.alldata), dtype=bool)
        for i, col in enumerate(columns):
            if isinstance(col[1], str):
                res = self.alldata[col[1]]
            else:
                if len(col[1]) == 3:
                    if isinstance(col[1][0], str):
                        col1 = self.alldata[col[1][0]]
                    elif isinstance(col[1][0], (list, tuple)):
                        col1 = self.alldata[col[1][0][0]][:, col[1][0][1]]
                    else:
                        col1 = col[1][0]

                    if isinstance(col[1][1], str):
                        col2 = self.alldata[col[1][1]]
                    elif isinstance(col[1][1], (list, tuple)):
                        col2 = self.alldata[col[1][1][0]][:, col[1][1][1]]
                    else:
                        col2 = col[1][1]

                    if col[1][2] == "-":
                        res = col1 - col2
                    elif col[1][2] == "+":
                        res = col1 + col2
                    elif col[1][2] == "*":
                        res = col1 * col2
                    elif col[1][2] == "/":
                        res = col1 / col2
                    elif col[1][2] == "SQSUM":
                        res = np.sqrt(col1**2. + col2**2.)
                    else:
                        raise KeyError("only + - * / are supported at the moment")

                elif len(col[1]) == 2:
                    res = self.alldata[col[1][0]][:, col[1][1]]
                else:
                    raise KeyError

            self.features[col[0]] = res.astype("float64")
            #
            if limits is not None:
                self.inds &= (self.features[col[0]] > limits[i][0]) & (self.features[col[0]] < limits[i][1])

        self.features = self.features[self.inds]

        try:
            self.weights = self.alldata["WEIGHT"][self.inds]
        except:
            self.weights = pd.Series(data=np.ones(len(self.features)), name="WEIGHT")

        for i, col in enumerate(columns):
            if logs is not None and logs[i]:
              self.features[col[0]] = np.log10(self.features[col[0]])

    def to_kde(self, **kwargs):
        """Convert to KDE container"""
        res = KDEContainer(self.features, weights=self.weights)
        return res


class FeatureSpaceContainer(BaseContainer):
    def __init__(self, info):
        """
        Container for galaxy features as in number profiles and target -- galaxy pair samples

        This has a radial galaxy number profile in bins described by rcens, redges, rareas

        A SurveyData and TargetData object and a number profile in the radial bins, as well as a galaxy sample profile

        The alldata of the BaseContainer parent class is created from the samples

        Parameters
        ----------
        info: indexer.IndexedDataContainer
            the result fo the indexing operation
        """

        BaseContainer.__init__(self)

        self.rcens = info.rcens
        self.redges = info.redges
        self.rareas = info.rareas

        self.survey = info.survey
        self.target = info.target

        self.numprof = info.numprof
        self.samples = info.samples

        valid_elements = np.nonzero([(len(tmp) > 0) for tmp in self.samples])[0].astype(int)
        if len(valid_elements) != len(self.samples):
            self.alldata = pd.concat(np.array(self.samples)[valid_elements]).reset_index(drop=True)
        else:
            self.alldata = pd.concat(self.samples).reset_index(drop=True)

        self.nobj = self.target.nrow

    def surfdens(self, icol=0, scaler=1):
        """Calculate the surface density of galaxies"""
        if self.logs[icol]:
            arr = 10**self.features.values[:, icol]
        else:
            arr = self.features.values[:, icol]
        vals = np.histogram(arr, bins=self.redges, weights=self.weights)[0] / self.nobj / self.rareas * scaler
        return vals

    def downsample(self, nmax=10000, r_key="LOGR", nbins=40, rng=None, **kwargs):
        """
        Radially balanced downsampling

        Parameters
        ----------
        nmax: int
            max number of objects to keep
        r_key: str
            key for the radial position column in the dataframe
        nbins: number of bins
        rng: np.random.RandomState
        kwargs

        Returns
        -------
        KDEContainer
            a downsampled version of the dataset in KDE container format
        """

        if rng is None:
            rng = np.random.RandomState()

        rarr = self.features[r_key]
        # rbins = np.sort(rng.uniform(low=rarr.min(), high=rarr.max(), size=nbins+1))
        rbins = np.linspace(rarr.min(), rarr.max(), nbins+1)

        tmp_features = []
        tmp_weights = []
        for i, tmp in enumerate(rbins[:-1]):
            selinds = (self.features[r_key] > rbins[i]) & (self.features[r_key] < rbins[i + 1])
            vals = self.features.loc[selinds]
            ww = self.weights.loc[selinds]

            if len(vals) < nmax:
                tmp_features.append(vals)
                tmp_weights.append(ww)
            else:
                inds = np.arange(len(vals))
                pp = ww / ww.sum()
                chindex = rng.choice(inds, size=nmax, replace=False, p=pp)

                newvals = vals.iloc[chindex]
                newww = ww.iloc[chindex] * len(ww) / nmax

                tmp_features.append(newvals)
                tmp_weights.append(newww)

        features = pd.concat(tmp_features)
        weights = pd.concat(tmp_weights)

        res = KDEContainer(features, weights=weights)
        # res = DualContainer(features.columns, **kwargs)
        # res.set_data(features, weights=weights)
        return res


class DeepFeatureContainer(BaseContainer):
    def __init__(self, data):
        """Feature container for DeepField datasets, This has no radial element, so it's much simpler"""
        BaseContainer.__init__(self)
        self.alldata = data
        self.weights = pd.Series(data=np.ones(len(self.alldata)), name="WEIGHT")

    @classmethod
    def from_file(cls, fname, flagsel=True):
        """
        Reads from file

        Parameters
        ----------
        fname: str
            file name to read from FITS or PANDAS with data key
        flagsel: bool
            if true, uses the 'flags" == 0 check to select entries
        """
        if ".fit" in fname:
            _deep = fio.read(fname)
        else:
            _deep = pd.read_hdf(fname, key="data").to_records()

        if flagsel:
            inds = _deep["flags"] == 0
            deep = _deep[inds]
        else:
            deep = _deep
        return cls(deep)


##########################################################################

def construct_wide_container(dataloader, settings, nbins=100, nmax=5000, seed=None, drop=None, **kwargs):
    fsc = FeatureSpaceContainer(dataloader)
    fsc.construct_features(**settings)
    # cont = fsc.to_dual(r_normalize=r_normalize)
    cont_small = fsc.downsample(nbins=nbins, nmax=nmax, kwargs=kwargs)
    cont_small.set_seed(seed)
    cont_small.shuffle()
    if drop is not None:
        cont_small.drop_col(drop)
    # cont_small.standardize_data()
    settings = copy.copy(settings)
    settings.update({"container": cont_small})
    return settings


def construct_deep_container(data, settings, seed=None, frac=1., drop=None):
    fsc = DeepFeatureContainer(data)
    fsc.construct_features(**settings)
    cont = fsc.to_kde()
    if drop is not None:
        cont.drop_col(drop)
    cont.set_seed(seed)
    cont.sample(frac=frac)
    # cont.standardize_data()
    settings = copy.copy(settings)
    settings.update({"container": cont})
    return settings

##########################################################################

def make_classifier_infodicts(wide_cr_clust, wide_r_ref, wide_cr_rands,
                              deep_c, deep_smc, columns,
                              nsamples=1e5, nchunks=1, bandwidth=0.1,
                              rmin=None, rmax=None, rcol="LOGR"):

    deep_smc_emu = deep_smc["container"]
    deep_smc_emu.standardize_data()
    deep_smc_emu.construct_kde(bandwidth)

    wide_r_emu = wide_r_ref["container"]
    wide_r_emu.standardize_data()
    wide_r_emu.construct_kde(bandwidth)

    samples_smc = deep_smc_emu.random_draw(nsamples)
    samples_r = wide_r_emu.random_draw(nsamples, rmin=rmin, rmax=rmax)
    samples = pd.merge(samples_smc, samples_r, left_index=True, right_index=True)
    sample_inds = partition(list(samples.index), nchunks)

    deep_smc_emu.drop_kde()
    wide_r_emu.drop_kde()
    infodicts = []
    for i in np.arange(nchunks):
        info = {
            "columns": columns,
            "bandwidth": bandwidth,
            "wide_cr_clust": wide_cr_clust,
            "wide_cr_rands": wide_cr_rands,
            "deep_c": deep_c,
            "wide_r_ref": wide_r_ref,
            "sample": samples.loc[sample_inds[i]],
            "rmin": rmin,
            "rmax": rmax,
        }
        infodicts.append(info)
    return infodicts, samples


def calc_scores2(info):
    scores = pd.DataFrame()
    try:
        columns = info["columns"]
        bandwidth = info["bandwidth"]
        sample = info["sample"]

        scores = pd.DataFrame()

        dc_emu = info["deep_c"]["container"]
        dc_emu.standardize_data()
        dc_emu.construct_kde(bandwidth=bandwidth)
        _score, _jacobian = dc_emu.score_samples(sample[columns["cols_dc"]])
        scores["dc"] = _score
        scores["dc_jac"] = _jacobian
        # scores["dc_jac"] = 1.

        wr_emu = info["wide_r_ref"]["container"]
        wr_emu.standardize_data()
        wr_emu.construct_kde(bandwidth=bandwidth)
        _score, _jacobian = wr_emu.score_samples(sample[columns["cols_wr"]])
        scores["wr"] = _score
        scores["wr_jac"] = _jacobian
        # scores["wr_jac"] = 1.

        wcr_emu = info["wide_cr_clust"]["container"]
        wcr_emu.standardize_data()
        wcr_emu.construct_kde(bandwidth=bandwidth)
        _score, _jacobian = wcr_emu.score_samples(sample[columns["cols_wcr"]])
        scores["wcr_clust"] = _score
        scores["wcr_clust_jac"] = _jacobian


        wcr_emu = info["wide_cr_rands"]["container"]
        wcr_emu.standardize_data()
        wcr_emu.construct_kde(bandwidth=bandwidth)
        _score, _jacobian = wcr_emu.score_samples(sample[columns["cols_wcr"]])
        scores["wcr_rands"] = _score
        scores["wcr_rands_jac"] = _jacobian
        # scores["wcr_jac"] = 1.

    except KeyboardInterrupt:
        pass

    return scores


def run_scores2(infodicts):
    pool = mp.Pool(processes=len(infodicts))
    try:
        pp = pool.map_async(calc_scores2, infodicts)
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


def make_naive_infodicts(wide_cr, wide_r, deep_c, deep_smc, columns,
                         nsamples=1e5, nchunks=1, bandwidth=0.1,
                         rmin=None, rmax=None, rcol="LOGR"):
    deep_smc_emu = deep_smc["container"]
    deep_smc_emu.standardize_data()
    deep_smc_emu.construct_kde(bandwidth)

    wide_r_emu = wide_r["container"]
    wide_r_emu.standardize_data()
    wide_r_emu.construct_kde(bandwidth)

    samples_smc = deep_smc_emu.random_draw(nsamples)
    samples_r = wide_r_emu.random_draw(nsamples, rmin=rmin, rmax=rmax)
    samples = pd.merge(samples_smc, samples_r, left_index=True, right_index=True)
    sample_inds = partition(list(samples.index), nchunks)

    deep_smc_emu.drop_kde()
    wide_r_emu.drop_kde()
    infodicts = []
    for i in np.arange(nchunks):
        info = {
            "columns": columns,
            "bandwidth": bandwidth,
            "wide_cr": wide_cr,
            "deep_c": deep_c,
            "wide_r": wide_r,
            "sample": samples.loc[sample_inds[i]],
            "rmin": rmin,
            "rmax": rmax,
        }
        infodicts.append(info)
    return infodicts, samples


def calc_scores(info):
    try:
        columns = info["columns"]
        bandwidth = info["bandwidth"]
        sample = info["sample"]

        scores = pd.DataFrame()

        dc_emu = info["deep_c"]["container"]
        dc_emu.standardize_data()
        dc_emu.construct_kde(bandwidth=bandwidth)
        _score, _jacobian = dc_emu.score_samples(sample[columns["cols_dc"]])
        scores["dc"] = _score
        scores["dc_jac"] = _jacobian
        # scores["dc_jac"] = 1.

        wr_emu = info["wide_r"]["container"]
        wr_emu.standardize_data()
        wr_emu.construct_kde(bandwidth=bandwidth)
        _score, _jacobian = wr_emu.score_samples(sample[columns["cols_wr"]])
        scores["wr"] = _score
        scores["wr_jac"] = _jacobian
        # scores["wr_jac"] = 1.

        wcr_emu = info["wide_cr"]["container"]
        wcr_emu.standardize_data()
        wcr_emu.construct_kde(bandwidth=bandwidth)
        _score, _jacobian = wcr_emu.score_samples(sample[columns["cols_wcr"]])
        scores["wcr"] = _score
        scores["wcr_jac"] = _jacobian
        # scores["wcr_jac"] = 1.

    except KeyboardInterrupt:
        pass

    return scores


def run_scores(infodicts):
    pool = mp.Pool(processes=len(infodicts))
    try:
        pp = pool.map_async(calc_scores, infodicts)
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


##############################################################

def read_concentric(score_path_expr, m_factor=20, seed=6):
    fname_scores = np.sort(glob.glob(score_path_expr))
    fname_samples = []
    for _fname in fname_scores:
        fname_samples.append(_fname.replace("scores", "samples"))

    samples = []
    for _fname in fname_samples:
        _tab = fio.read(_fname)
        _tab = pd.DataFrame.from_records(_tab)
        samples.append(_tab)
    samples = pd.concat(samples)

    scores = []
    for _fname in fname_scores:
        _tab = fio.read(_fname)
        _tab = pd.DataFrame.from_records(_tab)
        scores.append(_tab)
    scores = pd.concat(scores)

    dc_score = np.exp(scores["dc"]) * np.abs(scores["dc_jac"])
    wr_score = np.exp(scores["wr"]) * np.abs(scores["wr_jac"])
    wcr_score = np.exp(scores["wcr"]) * np.abs(scores["wcr_jac"])

    rng = np.random.RandomState(seed)
    uniform = rng.uniform(0, 1, len(samples))
    p_proposal = m_factor * dc_score * wr_score
    p_ref = wcr_score

    inds = uniform < p_ref / p_proposal
    resamples = samples[inds]
    return resamples
