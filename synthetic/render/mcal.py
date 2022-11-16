import ngmix
from ngmix.medsreaders import NGMixMEDS
import numpy as np
import os
import time
from astropy.table import Table, vstack, hstack
import fitsio as fio
import multiprocessing as mp


def infomaker(maxnum, medsfile, outfile_root, nchunk=100):
    lst = partition(np.arange(maxnum), n=nchunk)

    infodicts = []
    for i, tmp in enumerate(lst):
        info = {
            "medsfile": medsfile,
            "outfile": outfile_root + "_{:03d}.fits".format(i),
            "outdir": None,
            "start": tmp.min(),
            "end": tmp.max() + 1,  # due to how lists work, we need to include the one higher value here
        }
        infodicts.append(info)
    return infodicts


def collater(infodicts, ):
    tab = []
    names = []
    for i, info in enumerate(infodicts):

        fname = info['outfile']
        tmp = fio.read(fname)

        if i == 0:
            for i, name in enumerate(tmp.dtype.names):
                if 'err' not in name:
                    names.append(name)

        tab.append(tmp[names])

    return np.hstack(tab)


class MetacalFitter(object):
    def __init__(self, medsfile, seed=None):
        self.medsfile = medsfile

        # for some reason this is bugged if you
        self.meds = NGMixMEDS(medsfile)
        self.cat = self.meds.get_cat()
        self.Nobjs = len(self.cat)

        self.set_seed(seed)

        return

    def set_seed(self, seed=None):
        if seed is None:
            # current time in microseconds
            seed = int(1e6*time.time())

        self.seed = seed

        return

    def get_obs_list(self, iobj):
        return self.meds.get_obslist(iobj)

    def get_jacobians(self, iobj):
        Njac = len(self.meds.get_jacobian_list(iobj))

        jacobians = [self.meds.get_ngmix_jacobian(iobj, icutout)
                     for icutout in range(Njac)]

        return jacobians

    def get_obj_info(self, iobj):
        '''
        Setup object property dictionary used to compile fit params later on
        '''

        obj = self.meds[iobj]

        # Mcal object properties
        obj_info = {}

        obj_info['meds_indx'] = iobj
        obj_info['id'] = obj['id']
        obj_info['ra'] = obj['ra']
        obj_info['dec'] = obj['dec']
        obj_info['X_IMAGE'] = obj['X_IMAGE']
        obj_info['Y_IMAGE'] = obj['Y_IMAGE']

        return obj_info

    def get_prior(self, pixel_scale):
        '''
        pix_scale: The pixel scale of the image, in arcsec / pixel
        '''

        # This bit is needed for ngmix v2.x.x
        # won't work for v1.x.x
        #rng = np.random.RandomState(self.seed)

        # prior on ellipticity.  The details don't matter, as long
        # as it regularizes the fit.  This one is from Bernstein & Armstrong 2014

        g_sigma = 0.3
        g_prior = ngmix.priors.GPriorBA(g_sigma)

        # 2-d gaussian prior on the center
        # row and column center (relative to the center of the jacobian, which would be zero)
        # and the sigma of the gaussians

        # units same as jacobian, probably arcsec
        row, col = 0.0, 0.0
        row_sigma, col_sigma = pixel_scale, pixel_scale # use pixel_scale as a guess
        cen_prior = ngmix.priors.CenPrior(row, col, row_sigma, col_sigma)

        # T prior.  This one is flat, but another uninformative you might
        # try is the two-sided error function (TwoSidedErf).
        # NOTE: T units are arcsec^2 but can be slightly negative, especially for
        # stars if the PSF is mis-estimated

        Tminval = -1.0 # arcsec squared
        Tmaxval = 1000
        T_prior = ngmix.priors.FlatPrior(Tminval, Tmaxval)

        # similar for flux.  Make sure the bounds make sense for
        # your images

        Fminval = -1.e1
        Fmaxval = 1.e5
        F_prior = ngmix.priors.FlatPrior(Fminval, Fmaxval)

        # now make a joint prior.  This one takes priors
        # for each parameter separately
        prior = ngmix.joint_prior.PriorSimpleSep(
            cen_prior,
            g_prior,
            T_prior,
            F_prior
            )

        return prior

    def add_mcal_responsivities(self, mcal_res, mcal_shear):
        '''
        Compute and add the mcal responsivity values to the output
        result dict from get_metacal_result()

        NOTE: These are only for the selection-independent component!
        '''

        # Define full responsivity matrix, take inner product with shear moments
        r11 = (mcal_res['1p']['g'][0] - mcal_res['1m']['g'][0]) / (2*mcal_shear)
        r12 = (mcal_res['2p']['g'][0] - mcal_res['2m']['g'][0]) / (2*mcal_shear)
        r21 = (mcal_res['1p']['g'][1] - mcal_res['1m']['g'][1]) / (2*mcal_shear)
        r22 = (mcal_res['2p']['g'][1] - mcal_res['2m']['g'][1]) / (2*mcal_shear)

        R = [ [r11, r12], [r21, r22] ]
        Rinv = np.linalg.inv(R)
        gMC = np.dot(Rinv,
                     mcal_res['noshear']['g']
                     )

        MC = {
            'r11':r11, 'r12':r12,
            'r21':r21, 'r22':r22,
            'g1_MC':gMC[0], 'g2_MC':gMC[1]
        }

        mcal_res['MC'] = MC

        return mcal_res

    def mcal_dict2tab(self, mcal, obj_info):
        '''
        mcal is the dict returned by ngmix.get_metacal_result()

        obj_info is an array with MEDS identification info like id, ra, dec
        not returned by the function
        '''

        # Annoying, but have to do this to make Table from scalars
        for key, val in obj_info.items():
            obj_info[key] = np.array([val])

        tab_names = ['noshear', '1p', '1m', '2p', '2m','MC']
        for name in tab_names:
            tab = mcal[name]

            for key, val in tab.items():
                tab[key] = np.array([val])

            mcal[name] = tab

        id_tab = Table(data=obj_info)

        tab_noshear = Table(mcal['noshear'])
        tab_1p = Table(mcal['1p'])
        tab_1m = Table(mcal['1m'])
        tab_2p = Table(mcal['2p'])
        tab_2m = Table(mcal['2m'])
        tab_MC = Table(mcal['MC'])

        join_tab = hstack([id_tab, hstack([tab_noshear,
                                           tab_1p,
                                           tab_1m,
                                           tab_2p,
                                           tab_2m,
                                           tab_MC
                                           ],
                                          table_names=tab_names)
                           ]
                          )

        return join_tab

    def fit_obj(self, iobj, pars=None, ntry=4, psf_model='gauss',
                gal_model='gauss', vb=False):
        '''
        Run metacal fit for a single object of given index

        pars: mcal running parameters
        '''

        obj_info = self.get_obj_info(iobj)

        # Fits need a list of ngmix.Observation objects
        obs_list = self.get_obs_list(iobj)

        # Get pixel scale from image jacobian
        jac_list = self.get_jacobians(iobj)
        pixel_scale = jac_list[0].get_scale()

        if pars is None:
            # standard mcal run parameters
            mcal_shear = 0.01
            lm_pars = {'maxfev':2000, 'xtol':5.0e-5, 'ftol':5.0e-5}
            max_pars = {'method':'lm', 'lm_pars':lm_pars, 'find_center':True}
            metacal_pars = {'step':mcal_shear}
        else:
            mcal_shear = metacal_pars['step']
            max_pars = pars['max_pars']
            metacal_pars = pars['metacal_pars']

        prior = self.get_prior(pixel_scale)

        Tguess = 4*pixel_scale**2

        # setup run bootstrapper
        mcb = ngmix.bootstrap.MaxMetacalBootstrapper(obs_list)

        # Run the actual metacalibration step on the observed source
        mcb.fit_metacal(psf_model, gal_model, max_pars, Tguess, prior=prior,
                        ntry=ntry, metacal_pars=metacal_pars)

        mcal_res = mcb.get_metacal_result() # this is a dict

        # Add selection-independent responsitivities
        mcal_res = self.add_mcal_responsivities(mcal_res, mcal_shear)

        if vb is True:
            r11 = mcal_res['MC']['r11']
            r22 = mcal_res['MC']['r22']
            print(f'i={iobj}: R11: {r11:.3}; R22: {r22:.3} ')

        mcal_tab = self.mcal_dict2tab(mcal_res, obj_info)

        return mcal_tab


def run_mcal(info):
    medsfile = info["medsfile"]
    outfile = info["outfile"]
    outdir = info["outdir"]
    start = info["start"]
    end = info["end"]
    clobber = True
    vb = False  # args.vb # if True, prints out values of R11/R22 for every galaxy

    if outdir is not None:
        outfile = os.path.join(outdir, outfile)

    fitter = MetacalFitter(medsfile)

    if start is None:
        start = 0
    if end is None:
        end = fitter.Nobjs

    # Can set mcal parameters here if you want something beyond the default
    # in fit_obj()
    pars = None

    Tstart = time.time()

    mcal_tab = []
    for iobj in range(start, end):
        mcal_tab.append(
            fitter.fit_obj(iobj, pars=pars, vb=vb)
        )
    mcal_tab = vstack(mcal_tab)

    Tend = time.time()

    T = Tend - Tstart
    print(f'Total fitting and stacking time: {T} seconds')

    if vb is True:
        print(f'Writing out mcal results to {outfile}...')
    mcal_tab.write(outfile, overwrite=clobber)

    if vb is True:
        print('Done!')


def call_chunks(infodicts):
    """Technically this is inside the multiprocessing call, so it should not change any self properties"""
    for info in infodicts:
        run_mcal(info)



def multi_mcal(infodicts, nprocess=1):
    """
    OpenMP style parallelization
    Separates tasks into chunks, and passes each chunk for an independent process
    for serial evaulation via :py:func:`call_chunks`
    Parameters
    ----------
    infodict : dict
        A single list element returned from :py:func:`create_infodict`
    nprocess : int
        Number of processes (cores) to use. Maximum number is always set by ``len(infodicts)``
    """
    # at most as many processes can be used as there are independent tasks...
    if nprocess > len(infodicts):
        nprocess = len(infodicts)

    print('starting xshear calculations in ' + str(nprocess) + ' processes')
    fparchunks = partition(infodicts, nprocess)
    pool = mp.Pool(processes=nprocess)

    try:
        pp = pool.map_async(call_chunks, fparchunks)
        res = pp.get()  # apparently this counters a bug in the exception passing in python.subprocess...

    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        pool.join()
    else:
        pool.close()
        pool.join()


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