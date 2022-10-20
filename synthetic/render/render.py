"""
Galsim rendering backend

the first stretch goal is to build a stand alone renderer which is encapsulated in a class. The

# TODO add input catalog specification here

# TODO add color composite image creator function



"""

import numpy as np
import galsim
import ngmix
import multiprocessing as mp
from ..tools import partition


def draw_info(info):
    bdf_pars = info["bdf_pars"]
    psf = info["psf"]
    pixel_scale = info["pixel_scale"]
    x_cen = info["x_cen"]
    y_cen = info["y_cen"]
    offset = info["offset"]

    galmaker = ngmix.gmix.GMixBDF(bdf_pars)
    gs_profile = galmaker.make_galsim_object()
    final_gal = galsim.Convolve([psf, gs_profile])
    stamp_size = final_gal.getGoodImageSize(pixel_scale)

    bound = galsim.BoundsI(x_cen - stamp_size / 2 + 1, x_cen + stamp_size / 2,
                           y_cen - stamp_size / 2 + 1, y_cen + stamp_size / 2)

    stamp = galsim.ImageF(bound, scale=pixel_scale, )
    final_gal.drawImage(stamp, offset=offset, )

    return stamp, bound, info["id"]

def call_chunks(infodicts):
    """Technically this is inside the multiprocessing call, so it should not change any self properties"""
    stamps = []
    bounds = []
    ids = []
    for info in infodicts:
        stamp, bound, i = draw_info(info)
        stamps.append(stamp)
        bounds.append(bounds)
        ids.append(i)

    return stamps, bounds, ids

#TODO make more elegant via galsim

class DrawField(object):
    def __init__(self, canvas_size, catalog, band="g", pixel_scale=0.264, sky_level=1.e2, psf_fwhm=0.9):
        """
        This assumes a gaussian PSF
        """
        self.canvas_size = canvas_size
        self.canvas_cen = np.array((canvas_size / 2., canvas_size / 2.))
        self.catalog = catalog
        self.band = band
        self.pixel_scale = pixel_scale
        self.sky_level = sky_level
        self.psf_fwhm = psf_fwhm

        self.stamps = []
        self.stamps_bounds = []
        self.positions = []
        self.offsets = []

    def render(self, nprocess=10):
        self.prepare()
        self.make_infodicts()
        self.multi_render(nprocess)
        self.collate_stamps()

    def make_canvas(self):
        self.xx = self.catalog['X'] - self.canvas_size / 2
        self.yy = self.catalog['Y'] - self.canvas_size / 2
        self.canvas = galsim.ImageF(self.canvas_size, self.canvas_size, scale=self.pixel_scale)
        #self.canvas.array[:, :] = 0  # this might be redundant

    def make_wcs(self):

        # If you wanted to make a non-trivial WCS system, could set theta to a non-zero number
        theta = 0.0 * galsim.degrees
        dudx = np.cos(theta) * self.pixel_scale
        dudy = -np.sin(theta) * self.pixel_scale
        dvdx = np.sin(theta) * self.pixel_scale
        dvdy = np.cos(theta) * self.pixel_scale
        image_center = self.canvas.true_center
        affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=self.canvas.true_center)
        sky_center = galsim.CelestialCoord(ra=90. * galsim.degrees, dec=0 * galsim.degrees)

        self.wcs = galsim.TanWCS(affine, sky_center, units=galsim.arcsec)
        self.canvas.wcs = self.wcs


    def make_psf(self):
        self.psf = galsim.Gaussian(fwhm=self.psf_fwhm)
        self.image_epsf = self.psf.drawImage(scale=self.pixel_scale)

    def make_bdf_pars(self):
        self.bdf_pars = np.zeros((len(self.catalog), 7))
        # There might be a way to add xoffset and yoffset here also as 0 and 1 element
        self.bdf_pars[:, 2] = self.catalog["G1"][:]
        self.bdf_pars[:, 3] = self.catalog["G2"][:]
        self.bdf_pars[:, 4] = self.catalog["TSIZE"][:]
        self.bdf_pars[:, 5] = self.catalog["FRACDEV"][:]
        self.bdf_pars[:, 6] = self.catalog["FLUX_" + self.band.upper()][:]

    def make_positions(self):

        self.xx = self.catalog['X'][:] - self.canvas_cen[0]
        self.yy = self.catalog['Y'][:] - self.canvas_cen[1]

        self.x_cen = np.floor(self.xx)
        self.y_cen = np.floor(self.yy)

        self.offsets = np.vstack((self.xx - self.x_cen, self.yy - self.y_cen)).T

    def prepare(self):
        self.make_canvas()
        self.make_psf()
        self.make_bdf_pars()
        self.make_positions()

    def make_infodicts(self):
        """Prepare instruction set dictionaries to be passed for multiprocessing calculations"""
        self.infodicts = []
        for i in np.arange(len(self.catalog)):
            info = {
                "id": i,
                "bdf_pars": self.bdf_pars[i],
                "psf": self.psf,
                "pixel_scale": self.pixel_scale,
                "x_cen": self.x_cen[i],
                "y_cen": self.y_cen[i],
                "offset": self.offsets[i]
            }
            self.infodicts.append(info)

    def multi_render(self, nprocess=1):
        """
        OpenMP style parallelization for xshear tasks
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
        if nprocess > len(self.infodicts):
            nprocess = len(self.infodicts)

        print('starting xshear calculations in ' + str(nprocess) + ' processes')
        fparchunks = partition(self.infodicts, nprocess)
        pool = mp.Pool(processes=nprocess)

        self.stamps, self.bounds, self.ids = [], [], []
        try:
            pp = pool.map_async(call_chunks, fparchunks)
            res = pp.get(172800)  # apparently this counters a bug in the exception passing in python.subprocess...

            for tmp in res:
                self.stamps += tmp[0]
                self.bounds += tmp[1]
                self.ids += tmp[2]

        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            pool.join()
        else:
            pool.close()
            pool.join()

    def collate_stamps(self):
        for i in np.arange(len(self.catalog)):
            stamp = self.stamps[i]
            bb = stamp.bounds & self.canvas.bounds
            self.canvas[bb] += stamp[bb]

    def add_icl(self, arr):
        self.canvas += arr






def scale_image(canvas):
    try:
        res = np.arcsinh(canvas) / canvas
    except:
        res = np.arcsinh(canvas.array) / canvas.array
    return






def color_composite():
    pass

