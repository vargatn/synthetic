"""
General Utilities
"""

import numpy as np
import pandas as pd
import math


def radial_bins(rmin, rmax, nbin):
    """
    Calculates nominal edges and centers for logarithmic radial bins(base10 logarithm)
    Edges and areas are exact, "center" values are estimated as CIRCUMFERENCE weighted
    mean radius
    Parameters
    ----------
    rmin : float
        inner edge
    rmax : float
        outer edge
    nbin : int
        number of bins
    Returns
    -------
    np.array, np.array, np.array
        centers, edges, areas
    """
    edges = np.logspace(math.log10(rmin), math.log10(rmax), nbin + 1,
                        endpoint=True)
    cens = np.array([(edges[i + 1] ** 3. - edges[i] ** 3.) * 2. / 3. /
                     (edges[i + 1] ** 2. - edges[i] ** 2.)
                     for i, edge in enumerate(edges[:-1])])

    areas = np.array([np.pi * (edges[i + 1] ** 2. - edges[i] ** 2.)
                      for i, val in enumerate(edges[:-1])])

    return cens, edges, areas


def to_pandas(recarr):
    """
    Converts potentially nested record array (such as a FITS Table) into Pandas DataFrame
    FITS tables sometimes have multidimensional columns, which are not supported for DataFrames
    Pandas DataFrames however provide many nice features, such as SQL speed database matchings.
    The approach is to flatten out multidimensional column [[COL]] into [COL_1, COL_2, ..., COL_N]
    Examples
    --------
    Just pass the loaded FITS table::
        import fitsio as fio
        import xpipe.io.catalogs as catalogs
        raw_data = fio.read("data.fits")
        data = catalogs.to_pandas(raw_data)
    Parameters
    ----------
    recarr : numpy.array
        array to be converted to DataFrame
    Returns
    -------
    pandas.DataFrame
        array as DataFrame
    """

    newarr = flat_copy(recarr)
    res = pd.DataFrame.from_records(newarr.byteswap().newbyteorder(), columns=newarr.dtype.names)
    return res


def flat_type(recarr):
    """
    Assigns the dtypes to the flattened array
    Parameters
    ----------
    recarr : numpy.array
        array to be converted to DataFrame
    Returns
    -------
    list
        dtypes of flattened array
    """

    newtype = []
    for dt in recarr.dtype.descr:
        if len(dt) == 3:
            for i in np.arange(dt[2][0]):
                newtype.append((dt[0] + '_' + str(i), dt[1]))
        else:
            newtype.append(dt)
    return newtype


def flat_copy(recarr):
    """
    Copies the record array into a new recarray which has only 1-D columns
    Parameters
    ----------
    recarr : numpy.array
        array to be converted to DataFrame
    Returns
    -------
    numpy.array
        array with 1-D columns
    """

    newtype = flat_type(recarr)
    newarr = np.zeros(len(recarr), dtype=newtype)

    oldnames = recarr.dtype.names
    j = 0
    for i, dt in enumerate(recarr.dtype.descr):
        if len(dt) == 3:
            for c in np.arange(dt[2][0]):
                #                 print newtype[j]
                newarr[newtype[j][0]] = recarr[oldnames[i]][:, c]
                j += 1

        else:
            #             print newtype[j]
            newarr[newtype[j][0]] = recarr[oldnames[i]]
            j += 1
    return newarr


def partition(lst, n):
    """
    Divides a list into N roughly equal chunks
    Examples
    --------
    Define some test list, and look at the obtained chunks with different :code:`n` values::
        >>> lst = np.arange(20)
        >>> lst
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
        17, 18, 19])
        >>> partition(lst, n=5)
        [array([0, 1, 2, 3]),
         array([4, 5, 6, 7]),
         array([ 8,  9, 10, 11]),
         array([12, 13, 14, 15]),
         array([16, 17, 18, 19])]
        >>> partition(lst, n=6)
        [array([0, 1, 2]),
         array([3, 4, 5, 6]),
         array([7, 8, 9]),
         array([10, 11, 12]),
         array([13, 14, 15, 16]),
         array([17, 18, 19])]
    As we can see, even when :code:`n` is not a divisor of :code:`len(lst)`, it returns
    roughly balanced chunks
    Parameters
    ----------
    lst : list
        list to split up
    n : int
        chunks to make
    Returns
    -------
    list of lists
        list of chunks
    """
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))]
            for i in range(n) ]

def partition(lst, n):
    """
    Divides a list into N roughly equal chunks
    Examples
    --------
    Define some test list, and look at the obtained chunks with different :code:`n` values::
        >>> lst = np.arange(20)
        >>> lst
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
        17, 18, 19])
        >>> partition(lst, n=5)
        [array([0, 1, 2, 3]),
         array([4, 5, 6, 7]),
         array([ 8,  9, 10, 11]),
         array([12, 13, 14, 15]),
         array([16, 17, 18, 19])]
        >>> partition(lst, n=6)
        [array([0, 1, 2]),
         array([3, 4, 5, 6]),
         array([7, 8, 9]),
         array([10, 11, 12]),
         array([13, 14, 15, 16]),
         array([17, 18, 19])]
    As we can see, even when :code:`n` is not a divisor of :code:`len(lst)`, it returns
    roughly balanced chunks
    Parameters
    ----------
    lst : list
        list to split up
    n : int
        chunks to make
    Returns
    -------
    list of lists
        list of chunks
    """
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))]
            for i in range(n) ]