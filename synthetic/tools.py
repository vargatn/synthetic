
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