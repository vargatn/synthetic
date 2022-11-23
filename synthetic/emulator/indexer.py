"""
Module to handle survey data Processing and IO
# TODO refactor all settings into somewhere else
# TODO this should be runnable for DES and LSST
# TODO logging should be set up for this package
# TODO Testing and mock data should be implemented for this data
"""
from __future__ import print_function


import glob
import numpy as np
import pandas as pd
import fitsio as fio
import healpy as hp
import copy
import multiprocessing as mp
import pickle

from .utils import to_pandas, radial_bins, partition
from .paths import setup_logger, logfile_info, config


BADVAL = -9999.

DEFAULT_FLAGS = [
    ("MOF_CM_FLAGS", "==", 0),
    ("MOF_CM_T", "in", (0., 100)),
    ("MOF_CM_MAG_CORRECTED_I", "in", (14, 24)),
    (("MOF_CM_MAG_CORRECTED_G", "-", "MOF_CM_MAG_CORRECTED_R"), "in", (-4, 4)),
    (("MOF_CM_MAG_CORRECTED_R", "-", "MOF_CM_MAG_CORRECTED_I"), "in", (-4, 4)),
    (("MOF_CM_MAG_CORRECTED_I", "-", "MOF_CM_MAG_CORRECTED_Z"), "in", (-4, 4)),
]


# logger = setup_logger("INDEXER", level=config["logging_level"], logfile_info=logfile_info)


def get_theta_edges(nbins, theta_min, theta_max, eps):
    """
    Creates logarithmically space angular bins which include +- EPS linear range around zero
    The binning scheme looks the following::
        theta_edges = [ -eps, eps, theta_min, ... , theta_max]
    hence there are in total :code:`nbins + 2` bins.
    Parameters
    ----------
    nbins: int
        number of radial bins
    theta_min: float
        start of log10 spaced bins
    theta_max: float
        end of log10 spaced bins
    eps: float
        linear padding around zero
    Returns
    -------
    theta_edges: np.array
        radial edges
    rcens: np.array
        cemters of radial rings (starting at theta_min)
    redges: np.array
        edges of radial rings (starting at theta_min)
    rareas: np.array
        2D areas of radial rings (starting at theta_min)
    """
    rcens, redges, rareas = radial_bins(theta_min, theta_max, nbins)
    theta_edges = np.concatenate((np.array([-eps, eps, ]), redges))
    # logger.debug("theta_edges " + str(theta_edges))
    # logger.debug("rcens " + str(rcens))
    # logger.debug("redges " + str(redges))
    # logger.debug("rareas " + str(redges))
    return theta_edges, rcens, redges, rareas


def shuffle(tab, rng):
    """
    Returns a shuffled version of the passed DataFrame
    Uses :py:meth:`subsample` in the backend
    Parameters
    ----------
    tab: pd.DataFrame
        input table
    rng: np.random.RandomState
        random number generator, if None uses np.random directly
    Returns
    -------
    pd.DataFrame
        shuffled table
    """
    # logger.debug("shuffling table in place")
    # print("shuffling table in place")
    return subsample(tab, len(tab), rng, replace=False)


def get_ndraw(nsample, nchunk):
    """
    Calculates length of each approximately equal partitions of list
    Parameters
    ----------
    nsample: int
        number of entries in original list
    nchunk: int
        number of desired chunks
    Returns
    -------
    np.array
        length of each partition
    """
    division = float(nsample) / float(nchunk)
    arr = np.array([int(round(division * (i+1))) - int(round(division * (i)))
                for i in range(nchunk) ])
    return arr


def subsample(tab, nrows=1000, rng=None, replace=False):
    """
    Choose rows randomly from pandas DataFrame
    Parameters
    ----------
    tab: pd.DataFrame
        input table
    nrows: int
        number of rows to choose, automatically capped at table length
    rng: np.random.RandomState
        random number generator, if None uses np.random directly
    replace: bool
        draw with replacement or not
    Returns
    -------
    pd.DataFrame
        random row subset of input table
    """

    if rng is None:
        rng = np.random.RandomState()

    nrows=np.min((len(tab), int(round(nrows))))
    allinds = np.arange(len(tab))
    inds = allinds[rng.choice(allinds, nrows, replace=replace)]
    # logger.debug("subsampling " + str(nrows) + " objects out of " + str(len(tab)))
    return tab.iloc[inds], inds


class TargetData(object):
    def __init__(self, fname, mode=None):
        """
        Wrapper for unified handling of clusters and random point tables
        Exposes richness, redshift, ra, dec columns
        Supports selecting subsets of based on richness or other parameters
        Parameters
        ----------
        fname: str
            File name for fits table to use
        mode: str
            "clust" or "rands", if None figures out automatically
        """

        self.fname = fname
        # logger.debug(self.fname)
        self.mode = mode
        # logger.debug(self.mode)

        _data = fio.read(self.fname)
        self.alldata = to_pandas(_data)
        self.data = self.alldata
        # logger.debug("data shape:" + str(self.alldata.shape))

        self.inds = None
        self.pars = None
        self.limits = None
        self.assign_values()
        # logger.info("initiated TargetDate in mode " + str(self.mode) + " from " + str(self.fname))

    def assign_values(self):
        """Tries to guess 'mode' and exposes richness and redshift columns"""
        if self.mode is not None:
            if self.mode == "clust":
                self.richness = self.data.LAMBDA_CHISQ
                self.redshift = self.data.Z_LAMBDA
            elif self.mode == "rands":
                self.richness = self.data.AVG_LAMBDAOUT
                self.redshift = self.data.ZTRUE
        else:
            try:
                self.richness = self.data.LAMBDA_CHISQ
                self.redshift = self.data.Z_LAMBDA
                self.mode = "clust"
            except:
                self.richness = self.data.AVG_LAMBDAOUT
                self.redshift = self.data.ZTRUE
                self.mode = "rands"

        # logger.debug("z: " + str(np.array(self.redshift)))
        # logger.debug("lambda: " + str(np.array(self.richness)))
        self.ra = self.data.RA
        self.dec = self.data.DEC
        self.nrow = len(self.data)
        # logger.info("Number of targets: " + str(self.nrow))

    def reset_data(self):
        """Resets data to original table"""
        self.data, self.inds = self.alldata, None
        self.assign_values()
        # logger.info("resetting TargetData with filename " + str(self.fname))

    def draw_subset(self, nrows, rng=None):
        """draw random to subset of rows"""
        self.data, self.inds = subsample(self.data, nrows, rng=rng)
        self.assign_values()
        # logger.info("drawing " + str(nrows) + " subset from  TargetData with filename " + str(self.fname))

    def select_inds(self, inds, bool=True):
        """
        Selects subset based on index
        Parameters
        ----------
        inds: np.array
            indexing array
        bool: bool
            whether indexing array is bool or integer
        """

        if bool:
            self.inds = np.nonzero(inds)
        else:
            self.inds = inds

        self.data = self.alldata.iloc[self.inds]
        self.nrow = len(self.data)
        # logger.info("selected inds (" + str(len(self.data)) + " subset) from  TargetData with filename " + str(self.fname))

    def select_range(self, pars, limits):
        """
        Selects single parameter bin from underlying data table
        In addition to columns, "redshift" and "richness" are also valid keys, these automatically refer to the
        appropriate column
        Parameters
        ----------
        pars: str or list
            Column name or list of Column names
        limits: list
            value limits for each column
        """

        self.reset_data()

        self.pars = pars
        self.limits = limits
        # logger.info("selecting subset from  TargetData with filename " + str(self.fname))
        # logger.info("pars:" + str(self.pars))
        # logger.info("limits:" + str(self.limits))
        bool_inds = np.ones(len(self.data), dtype=bool)
        for par, lim in zip(pars, limits):
            if par == "redshift":
                _ind = (self.redshift > lim[0]) & (self.redshift < lim[1])
                bool_inds[np.invert(_ind)] = False
            elif par == "richness":
                _ind = (self.richness > lim[0]) & (self.richness < lim[1])
                bool_inds[np.invert(_ind)] = False
            else:
                _ind = (self.alldata[par] > lim[0]) & (self.alldata[par] < lim[1])
                bool_inds[np.invert(_ind)] = False

        self.inds = np.nonzero(bool_inds)
        self.data = self.alldata.iloc[self.inds]
        self.nrow = len(self.data)
        self.assign_values()

    def to_dict(self):
        """Extracts metadata of self int a dictionary for lean storage"""
        info = {
            "fname": self.fname,
            "inds": self.inds,
            "pars": self.pars,
            "limits": self.limits,
            "mode": self.mode,
            "nrow": self.nrow

        }
        return info

    @classmethod
    def from_dict(cls, info):
        """recreate full object from dictionary"""
        res = cls(info["fname"], info["mode"])

        if (info["pars"] is not None):
            res.select_range(info["pars"], info["limits"])
        if (info["inds"] is not None):
            res.select_inds(info["inds"], bool=False)

        return res

    @classmethod
    def from_config(cls, mode, config):
        """
        Automatically reads from config
        Parameters
        ----------
        mode: str
            clust or rands
        config: dict
            Config dictionary
        """

        if mode == "clust":
            fname = config["catalogs"]["targets"]["clust"]
        elif mode == "rands":
            fname = config["catalogs"]["targets"]["rands"]
        else:
            raise KeyError("Currently only clust and rands mode is supported")
        # logger.info("constructing TargetData from config file")
        return cls(fname, mode)


class SurveyData(object):
    def __init__(self, fnames, nside=16):
        self.fnames = fnames
        self.nchunks = len(fnames)
        self.nside = nside

    def get_data(self, ind):
        fname = self.fnames[ind]
        print(fname)
        self.tab = pd.read_hdf(fname)
        self.itab = ind

    def drop_data(self):
        """Resets SurveyData table to None"""
        self.tab = None
        self.pixels = None
        # logger.info("resetting SurveyData")

    def lean_copy(self):
        """Returns a low-memory version of the SurveyData"""
        return SurveyData(self.fnames, self.nside)

    def to_dict(self):
        infodict = {
            "fnames": self.fnames,
            "nside": self.nside,
        }
        return infodict

    @classmethod
    def from_dict(cls, infodict):
        return cls(**infodict)


def convert_on_disk(fnames, nprocess, nside=16):
    nchunks = len(fnames)
    if nprocess > nchunks:
        nprocess = nchunks

    infodicts = []
    for fname in fnames:
        info = {
            "fname": fname,
            "nside": nside,
        }
        infodicts.append(info)
    info_chunks = partition(infodicts, nchunks)

    pool = mp.Pool(processes=nprocess)
    try:
        pp = pool.map_async(_convert_chunk_run, info_chunks)
        pp.get(86400)  # apparently this counters a bug in the exception passing in python.subprocess...
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        pool.join()
    else:
        pool.close()
        pool.join()

def _convert_chunk_run(chunks):
    try:
        for infodict in chunks:
            _converter(infodict)
    except KeyboardInterrupt:
        pass

def _converter(info):
    fname = info["fname"]
    nside = info["nside"]
    print("converting",fname)

    data = fio.read(fname)
    data = to_pandas(data)
    data["IPIX"] = hp.ang2pix(nside, data.RA, data.DEC, lonlat=True)
    data.nside = nside
    data.to_hdf(fname.replace("fits", "h5"), key="data")


def get_flags(tab, cutlist):
    flags = np.ones(len(tab), dtype=bool)
    for sel in cutlist:
        # calculate input column
        if isinstance(sel[0], str):
            col = tab[sel[0]]
        else:
            if sel[0][1] == "+":
                col = tab[sel[0][0]] + tab[sel[0][2]]
            elif sel[0][1] == "-":
                col = tab[sel[0][0]] - tab[sel[0][2]]
            elif sel[0][1] == "*":
                col = tab[sel[0][0]] * tab[sel[0][2]]
            elif sel[0][1] == "/":
                col = tab[sel[0][0]] / tab[sel[0][2]]
            else:
                raise KeyError

        # applying cut
        if sel[1] == "==":
            vals = (col == sel[2])
        elif sel[1] == ">":
            vals = (col > sel[2])
        elif sel[1] == "<":
            vals = (col < sel[2])
        elif sel[1] == "in":
            vals = (col >= sel[2][0]) & (col < sel[2][1])
        else:
            raise KeyError

        flags &= vals

    return flags


class MultiIndexer(object):
    def __init__(self, survey, target, fname_root, search_radius=360.,
                 nbins=50, theta_min=0.1, theta_max=100, eps=1e-3):

        self.survey = survey
        self.target = target
        self.fname_root = fname_root

        self.search_radius = search_radius

        self.nbins = nbins
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.eps = eps
        self.theta_edges, self.rcens, self.redges, self.rareas = get_theta_edges(nbins, theta_min, theta_max, eps)

    def _get_infodicts(self):
        nchunks = self.survey.nchunks
        infodicts = []
        for i in np.arange(nchunks):
            info = {
                "ind": i,
                "target": self.target.to_dict(),
                "survey": self.survey.to_dict(),
                "search_radius": self.search_radius,
                "nbins": self.nbins,
                "theta_edges": self.theta_edges,
                "rcens": self.rcens,
                "redges": self.redges,
                "rareas": self.rareas,
                "fname": self.fname_root + "_" + str(i) + ".p"
            }
            infodicts.append(info)
        return infodicts

    def run(self, nprocess=1):

        infodicts = self._get_infodicts()
        if nprocess > len(infodicts):
            nprocess = len(infodicts)

        info_chunks = partition(infodicts, nprocess)
        print("starting calculation in", len(info_chunks), "processes")

        pool = mp.Pool(processes=nprocess)
        try:
            pp = pool.map_async(_indexer_chunk_run, info_chunks)
            pp.get(86400)  # apparently this counters a bug in the exception passing in python.subprocess...
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            pool.join()
        else:
            pool.close()
            pool.join()


def _indexer_chunk_run(chunks):
    try:
        for infodict in chunks:
            print("running", infodict["fname"])
            maker = SurveyIndexer(**infodict.copy())
            maker.run(infodict["fname"])
    except KeyboardInterrupt:
        pass


class IndividualSurveyIndexer(object):
    def __init__(self, survey, target, theta_edges, rcens, redges, rareas, search_radius=360., nbins=50, ind=0,
                 flags=DEFAULT_FLAGS, **kwargs):

        if isinstance(survey, dict):
            self.survey = SurveyData.from_dict(survey)
        else:
            self.survey = survey

        if isinstance(target, dict):
            self.target = TargetData.from_dict(target)
        else:
            self.target = target
        self.nrow = self.target.nrow

        self.search_radius = search_radius

        self.nbins = nbins
        self.theta_edges = theta_edges
        self.rcens = rcens
        self.redges = redges
        self.rareas = rareas
        self.ind = ind

        self.flags = flags

        self._get_data()


    def _get_data(self):
        self.survey.get_data(self.ind)
        if self.flags is not None:
            flags = get_flags(self.survey.tab, self.flags)
            self.survey.tab = self.survey.tab[flags]
        else:
            self.survey.tab = self.survey.tab

    def run(self, fname="test"):
        print(fname)
        # pass
        # self.index()
        # result = self.draw_samples()
        # pickle.dump(result, open(fname, "wb"))

    def index(self):
        print("starting survey indexing")
        self.numprof = np.zeros(self.nbins + 2)
        self.numprofiles = np.zeros((self.target.nrow, self.nbins + 2))
        self.container = [[] for tmp in np.arange(self.nbins + 2)]
        print("indexing samples")
        for i in np.arange(self.target.nrow):
            print(str(i) + "\t / \t" + str(self.target.nrow))
            # logger.debug(str(i) + "/" + str(self.target.nrow))

            trow = self.target.data.iloc[i]
            tvec = hp.ang2vec(trow.RA, trow.DEC, lonlat=True)

            _radius = self.search_radius / 60. / 180. * np.pi
            dpixes = hp.query_disc(self.survey.nside, tvec, radius=_radius)

            gals = []
            for dpix in dpixes:
                cmd = "IPIX == " + str(dpix)
                gals.append(self.survey.tab.query(cmd))
            gals = pd.concat(gals)

            darr = np.sqrt((trow.RA - gals.RA) ** 2. + (trow.DEC - gals.DEC) ** 2.) * 60. # converting to arcmin
            gals["DIST"] = darr

            tmp = np.histogram(darr, bins=self.theta_edges)[0]
            self.numprof += tmp
            self.numprofiles[i] = tmp

            for j in np.arange(self.nbins + 2):
                cmd = str(self.theta_edges[j]) + " < DIST < " + str(self.theta_edges[j + 1])
                rsub = gals.query(cmd)
                self.container[j].append(rsub.index.values)

        self.indexes, self.counts = [], []
        for j in np.arange(self.nbins):
            _uniqs, _counts = np.unique(np.concatenate(self.container[j]), return_counts=True)
            self.indexes.append(_uniqs)
            self.counts.append(_counts)

        result = IndexedDataContainer(self.survey.lean_copy(), self.target.to_dict(),
                                      self.numprof, self.indexes, self.counts,
                                      self.theta_edges, self.rcens, self.redges, self.rareas)
        # logger.info("finished survey indexing")
        print("finished survey indexing")
        return result

    # def draw_samples(self, nsample=10000, rng=None):
    #     # logger.info("starting drawing random subsample with nsample=" + str(nsample))
    #     print("starting drawing random subsample with nsample=" + str(nsample))
    #
    #     if rng is None:
    #         rng = np.random.RandomState()
    #
    #     num_to_draw = np.min((self.numprof, np.ones(self.nbins + 2) * nsample), axis=0).astype(int)
    #     limit_draw = num_to_draw == nsample
    #
    #     self.sample_nrows = np.zeros(self.nbins + 2)
    #     samples = [[] for tmp in np.arange(self.nbins + 2)]
    #     self.samples = samples
    #     print("drawing samples")
    #     for i in np.arange(self.target.nrow):
    #         print(str(i) + "\t / \t" + str(self.target.nrow))
    #         # logger.debug(str(i) + "/" + str(self.target.nrow))
    #
    #         trow = self.target.data.iloc[i]
    #         tvec = hp.ang2vec(trow.RA, trow.DEC, lonlat=True)
    #
    #         _radius = self.search_radius / 60. / 180. * np.pi
    #         dpixes = hp.query_disc(self.survey.nside, tvec, radius=_radius)
    #
    #         gals = []
    #         for dpix in dpixes:
    #             cmd = "IPIX == " + str(dpix)
    #             gals.append(self.survey.tab.query(cmd))
    #         gals = pd.concat(gals)
    #
    #         darr = np.sqrt((trow.RA - gals.RA) ** 2. + (trow.DEC - gals.DEC) ** 2.) * 60.  # converting to arcmin
    #         gals["DIST"] = darr
    #
    #         digits = np.digitize(darr, self.theta_edges) - 1.
    #         gals["DIGIT"] = digits
    #
    #         for d in np.arange(self.nbins + 2):
    #             cmd = "DIGIT == " + str(d)
    #             rows = gals.query(cmd)
    #
    #             # tries to draw a subsample from around each cluster in each radial range
    #             if limit_draw[d]:
    #                 _ndraw = get_ndraw(nsample - self.sample_nrows[d], self.target.nrow - i)[0]
    #                 ndraw = np.min((_ndraw, len(rows)))
    #                 self.sample_nrows[d] += ndraw
    #                 if ndraw > 0:
    #                     ii = rng.choice(np.arange(len(rows)), ndraw, replace=False)
    #                     samples[d].append(rows.iloc[ii])
    #             else:
    #                 # print("  ",len(rows))
    #                 self.sample_nrows[d] += len(rows)
    #                 samples[d].append(rows)
    #
    #     for d in np.arange(self.nbins + 2):
    #         self.samples[d] = pd.concat(samples[d], ignore_index=True)
    #
    #     result = IndexedDataContainer(self.survey.lean_copy(), self.target.to_dict(),
    #                                   self.numprof, self.indexes, self.counts,
    #                                   self.theta_edges, self.rcens, self.redges, self.rareas,
    #                                   self.samples, self.sample_nrows)
    #     # logger.info("finished random draws")
    #     print("finished random draws")
    #
    #     return result


class SurveyIndexer(object):
    def __init__(self, survey, target, theta_edges, rcens, redges, rareas, search_radius=360., nbins=50, ind=0,
                 flags=DEFAULT_FLAGS, **kwargs):

        if isinstance(survey, dict):
            self.survey = SurveyData.from_dict(survey)
        else:
            self.survey = survey

        if isinstance(target, dict):
            self.target = TargetData.from_dict(target)
        else:
            self.target = target

        self.search_radius = search_radius

        self.nbins = nbins
        self.theta_edges = theta_edges
        self.rcens = rcens
        self.redges = redges
        self.rareas = rareas
        self.ind = ind

        self.flags = flags

        self._get_data()

    def _get_data(self):
        self.survey.get_data(self.ind)
        if self.flags is not None:
            flags = get_flags(self.survey.tab, self.flags)
            self.survey.tab = self.survey.tab[flags]
        else:
            self.survey.tab = self.survey.tab

    def run(self, fname="test"):
        print(fname)
        # pass
        self.index()
        result = self.draw_samples()
        pickle.dump(result, open(fname, "wb"))

    def index(self):

        print("starting survey indexing")
        self.numprof = np.zeros(self.nbins + 2)
        self.numprofiles = np.zeros((self.target.nrow, self.nbins + 2))
        self.container = [[] for tmp in np.arange(self.nbins + 2)]
        print("indexing samples")
        for i in np.arange(self.target.nrow):
            print(str(i) + "\t / \t" + str(self.target.nrow))
            # logger.debug(str(i) + "/" + str(self.target.nrow))

            trow = self.target.data.iloc[i]
            tvec = hp.ang2vec(trow.RA, trow.DEC, lonlat=True)

            _radius = self.search_radius / 60. / 180. * np.pi
            dpixes = hp.query_disc(self.survey.nside, tvec, radius=_radius)

            gals = []
            for dpix in dpixes:
                cmd = "IPIX == " + str(dpix)
                gals.append(self.survey.tab.query(cmd))
            gals = pd.concat(gals)

            darr = np.sqrt((trow.RA - gals.RA) ** 2. + (trow.DEC - gals.DEC) ** 2.) * 60. # converting to arcmin
            gals["DIST"] = darr

            tmp = np.histogram(darr, bins=self.theta_edges)[0]
            self.numprof += tmp
            self.numprofiles[i] = tmp

            for j in np.arange(self.nbins + 2):
                cmd = str(self.theta_edges[j]) + " < DIST < " + str(self.theta_edges[j + 1])
                rsub = gals.query(cmd)
                self.container[j].append(rsub.index.values)

        self.indexes, self.counts = [], []
        for j in np.arange(self.nbins):
            _uniqs, _counts = np.unique(np.concatenate(self.container[j]), return_counts=True)
            self.indexes.append(_uniqs)
            self.counts.append(_counts)

        result = IndexedDataContainer(self.survey.lean_copy(), self.target.to_dict(),
                                      self.numprof, self.indexes, self.counts,
                                      self.theta_edges, self.rcens, self.redges, self.rareas)
        # logger.info("finished survey indexing")
        print("finished survey indexing")
        return result

    def draw_samples(self, nsample=10000, rng=None):
        # logger.info("starting drawing random subsample with nsample=" + str(nsample))
        print("starting drawing random subsample with nsample=" + str(nsample))

        if rng is None:
            rng = np.random.RandomState()

        num_to_draw = np.min((self.numprof, np.ones(self.nbins + 2) * nsample), axis=0).astype(int)
        limit_draw = num_to_draw == nsample

        self.sample_nrows = np.zeros(self.nbins + 2)
        samples = [[] for tmp in np.arange(self.nbins + 2)]
        self.samples = samples
        print("drawing samples")
        for i in np.arange(self.target.nrow):
            print(str(i) + "\t / \t" + str(self.target.nrow))
            # logger.debug(str(i) + "/" + str(self.target.nrow))

            trow = self.target.data.iloc[i]
            tvec = hp.ang2vec(trow.RA, trow.DEC, lonlat=True)

            _radius = self.search_radius / 60. / 180. * np.pi
            dpixes = hp.query_disc(self.survey.nside, tvec, radius=_radius)

            gals = []
            for dpix in dpixes:
                cmd = "IPIX == " + str(dpix)
                gals.append(self.survey.tab.query(cmd))
            gals = pd.concat(gals)

            darr = np.sqrt((trow.RA - gals.RA) ** 2. + (trow.DEC - gals.DEC) ** 2.) * 60.  # converting to arcmin
            gals["DIST"] = darr

            digits = np.digitize(darr, self.theta_edges) - 1.
            gals["DIGIT"] = digits

            for d in np.arange(self.nbins + 2):
                cmd = "DIGIT == " + str(d)
                rows = gals.query(cmd)

                # tries to draw a subsample from around each cluster in each radial range
                if limit_draw[d]:
                    _ndraw = get_ndraw(nsample - self.sample_nrows[d], self.target.nrow - i)[0]
                    ndraw = np.min((_ndraw, len(rows)))
                    self.sample_nrows[d] += ndraw
                    if ndraw > 0:
                        ii = rng.choice(np.arange(len(rows)), ndraw, replace=False)
                        samples[d].append(rows.iloc[ii])
                else:
                    # print("  ",len(rows))
                    self.sample_nrows[d] += len(rows)
                    samples[d].append(rows)

        for d in np.arange(self.nbins + 2):
            self.samples[d] = pd.concat(samples[d], ignore_index=True)

        result = IndexedDataContainer(self.survey.lean_copy(), self.target.to_dict(),
                                      self.numprof, self.indexes, self.counts,
                                      self.theta_edges, self.rcens, self.redges, self.rareas,
                                      self.samples, self.sample_nrows)
        # logger.info("finished random draws")
        print("finished random draws")

        return result


class IndexedDataContainer(object):
    def __init__(self, survey=None, target=None, numprof=None, indexes=None, counts=None, theta_edges=None, rcens=None,
                 redges=None, rareas=None, samples=None, samples_nrows=None):
        """
        Container for Indexed Survey Data
        It serves only as a data wrapper which can be pickled easily. The bulk of the survey data or target data
        should not be contained, and can be dropped and recovered when necessary.
        All parameters are class variables with the same name.
        Parameters
        ----------
        survey: :py:meth:`SurveyData` instance
            Container for the survey data
        target: :py:meth:`TargetData` instance
            Container for the target data
        numprof: np.array
            number profile of objects around the targets
        indexes: list
            index of unique galaxies at each radial bin around targets
        counts: list of list
            multiplicity of unique galaxies at each radial bin around targets
        theta_edges: np.array
                radial edges
        rcens: np.array
                centers of radial rings (starting at theta_min)
        redges: np.array
            edges of radial rings (starting at theta_min)
        rareas: np.array
            2D areas of radial rings (starting at theta_min)
        samples: list of pd.DataFrame
            table of random galaxy draws from each radial bin (capped in size at :code:`nsamples`)
        samples_nrows: np.array
            number of galaxies drawn from each radial bin
        """
        self.survey = survey
        self.target = target
        self.numprof = numprof
        self.indexes = indexes
        self.counts = counts
        self.theta_edges = theta_edges
        self.rcens = rcens
        self.redges = redges
        self.rareas = rareas
        self.samples = samples
        self.samples_nrows = samples_nrows

    def expand_data(self):
        """Recover all data from disk"""
        self.target = TargetData.from_dict(self.target)
        # self.survey.read_all()

    def drop_data(self):
        """Drops all data and keeps only necessary values"""
        self.survey = self.survey.drop_data()
        self.target = self.target.to_dict()


class MultiDataLoader(object):
    def __init__(self, fnames, force_target):
        self.fnames = fnames

        self._get_info(force_target)

    def _get_info(self, force_target):
        pp = pickle.load(open(self.fnames[0], "rb"))

        self.rcens = pp.rcens.copy()
        self.redges = pp.redges.copy()
        self.rareas = pp.rareas.copy()
        self.target_nrow = pp.target["nrow"]
        self.nrbins = len(pp.numprof)

        if force_target:
            self.target = force_target
        else:
            self.target = TargetData.from_dict(pp.target)
        self.survey = pp.survey

    def collate_samples(self):
        _samples = [[] for tmp in np.arange(self.nrbins)]
        self.numprof = np.zeros(self.nrbins)
        pps = []
        for i, fname in enumerate(self.fnames):
            print(fname)
            pp = pickle.load(open(fname, "rb"))
            pps.append(pp)

            self.numprof += pp.numprof
            for j in np.arange(self.nrbins):
                if len(pp.samples[j]):
                    _samples[j].append(pp.samples[j])

        self.samples = [[] for tmp in np.arange(self.nrbins)]
        for j in np.arange(self.nrbins):
            print(j)
            if len(_samples[j]):
                self.samples[j] = pd.concat(_samples[j]).reset_index(drop=True)
                self.samples[j]["WEIGHT"] = self.numprof[j] / len(self.samples[j])

    def to_cont(self):
        cont = IndexedDataContainer(self.survey, self.target, numprof=self.numprof,
                                    rcens=self.rcens, redges=self.redges, rareas=self.rareas, samples=self.samples)
        return cont



