import galsim
import numpy as np

import astropy.constants as constants
import astropy.units as u


def sigma_crit_inv(zclust, z, cosmo):
    if z > zclust:
        prefac = (4. * np.pi * constants.G / (constants.c ** 2.)).to(u.pc / u.Msun)

        Ds = cosmo.angular_diameter_distance(z).to(u.pc)
        Dl = cosmo.angular_diameter_distance(zclust).to(u.pc)
        Dls = cosmo.angular_diameter_distance_z1z2(zclust, z).to(u.pc)

        val = prefac * Dl * Dls / Ds
        resval = np.max((0., val.value))
    else:
        resval = 0
    return resval


class Shear(object):
    def __init__(self, canvas, image_epsf, maskmap, pixel_scale):
        """ITs very important for the images and stamps and psf images to have a proper image format, includeing WCS"""
        self.canvas = canvas.copy()
        self.image_epsf = galsim.ImageF(image_epsf, wcs=galsim.PixelScale(pixel_scale))
        self.BADVAL = -9999
        self.maskmap = galsim.ImageF(maskmap, wcs=galsim.PixelScale(pixel_scale)).copy()
        self.pixel_scale = pixel_scale

    def extract_stamps(self, centers, imasks, sizes):

        self.stamps = []
        self.masks = []
        for i in np.arange(len(centers)):
            half_size = sizes[i] // 2 + 2
            cen = centers[i]
            bb = galsim.bounds.BoundsI(np.round(cen[0]) - half_size,
                                       np.round(cen[0]) + half_size - 1,
                                       np.round(cen[1]) - half_size,
                                       np.round(cen[1]) + half_size - 1)

            stamp = self.canvas[bb].array.copy()
            mask = self.maskmap[bb].array.copy()
            mm = (mask != imasks[i]) & (mask != 0)
            muniques = np.unique(mask)
            stamp[mm] = 0
            mask[mm] = 0

            self.masks.append(self.maskmap[bb].array.copy())
            self.stamps.append(galsim.ImageF(stamp, wcs=galsim.PixelScale(self.pixel_scale)))

    def estimate_shear(self, sky_var=0, shear_est="KSB"):
        """
        shear modes = REGAUSS’, ‘LINEAR’, ‘BJ’, or ‘KSB’
        which is raw or canvas, default to canvas
        """

        #         stamps = self.raw_stamps
        # if which == "canvas":
        stamps = self.stamps

        self.logs = []
        self.shears = []
        self.success = []
        self.fluxes = []
        self.shears_error = []
        for i, im in enumerate(stamps):
            try:
                res = galsim.hsm.EstimateShear(im, self.image_epsf, shear_est=shear_est, sky_var=sky_var)
                #                 print(res)
                self.logs.append(res)
                self.success.append(True)

                tmp = np.array([res.corrected_g1, res.corrected_g2])

                self.shears.append(tmp)
                self.fluxes.append(res.moments_amp * (self.pixel_scale) ** 2)
                self.shears_error.append(res.corrected_shape_err)
            #                             print("here")
            except:
                self.success.append(False)
                self.logs.append(None)
                self.fluxes.append(None)
                tmp = np.array([self.BADVAL, self.BADVAL])
                self.shears.append(tmp)
                self.shears_error.append(None)

        self.fluxes = np.array(self.fluxes)
        self.shears = np.array(self.shears)
        self.success = np.array(self.success)
        self.shears_error = np.array(self.shears_error)

