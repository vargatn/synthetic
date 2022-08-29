"""
Interface to DES ICL measurements

Based on Varga et al 2021


# TODO encapsulate the data product better
"""

import numpy as np

import pandas as pd
import fitsio as fio
import scipy.special
import scipy.optimize as optimize

import astropy.cosmology as cosmology
import astropy.units as u

def sigmoid(arr, pos, k,):
    return 1 / (1 + np.exp((arr - pos) * k))

def tomag(flux):
    return 30. -2.5 * np.log10(flux)

def toflux(mag):
    return 10**((mag - 30) / (-2.5))


class DrawICL(object):
    def __init__(self, mass, z, bcg, galpath, mstarpath, jk_profile_root,
                 pixel_scale=0.264, canvas_size=5000, H0=70, Om0=0.3):
        self.mass = mass
        self.z = z
        self.bcg = bcg
        self.pixel_scale = pixel_scale
        self.canvas_size = canvas_size
        self.galpath = galpath
        self.mstarpath = mstarpath
        self.jk_profile_root = jk_profile_root
        self.H0 = H0
        self.Om0 = Om0

    def load_data(self):
        self.procmag = ProcMag(self.galpath, self.mstarpath)
        self.iclprof = ICLProf(self.jk_profile_root, self.procmag, H0=self.H0, Om0=self.Om0)

    def prepare_icl(self):
        sval = np.sqrt(self.bcg["size"])

        xsize = self.canvas_size
        ysize = self.canvas_size
        xcen = xsize / 2.
        ycen = ysize / 2.

        _xx = np.array([np.arange(xsize), ] * ysize) - xcen
        _yy = np.array([np.arange(ysize), ] * xsize).T - ycen
        xx = (1 - self.bcg["g1"]) * _xx - self.bcg["g2"] * _yy
        yy = - self.bcg["g2"] * _xx + (1 + self.bcg["g1"]) * _yy
        self.rmatrix = np.sqrt(xx ** 2. + yy ** 2.)

        _rvals = self.rmatrix.flatten()
        _iclarr = self.iclprof.icl_pixes(_rvals, mass=self.mass, z=self.z)
        iclarr = _iclarr * (1 - sigmoid(_rvals * self.pixel_scale, k=1, pos=sval))

        self.iclvals = iclarr.reshape((ysize, xsize))
        self.iclvals[self.iclvals < 1.] = 1.

    def colorize_icl(self):
        self.magvals_i = tomag(self.iclvals)
        self.magvals_g = self.bcg["color_gr"] + self.bcg["color_ri"] + self.magvals_i
        self.magvals_r = self.bcg["color_ri"] + self.magvals_i
        self.magvals_z = self.magvals_i - self.bcg["color_iz"]

        self.flux_g = toflux(self.magvals_g).T
        self.flux_r = toflux(self.magvals_r).T
        self.flux_i = toflux(self.magvals_i).T
        self.flux_z = toflux(self.magvals_z).T

    def get_icl(self):

        self.load_data()
        self.prepare_icl()
        self.colorize_icl()


def int_schechter(l1, l2):
    """
    # int_a^b x x^{alpha} e^{-x} dx ~ (incomplete_Gamma(alpha-1,a)-incomplete_Gamma(alpha-1,b))
    """
    alpha = -1
    return scipy.special.gammainc(alpha + 2, l2) - scipy.special.gammainc(alpha + 2, l1)


class ProcMag(object):
    def __init__(self, galpath, mstarpath):
        self.galpath = galpath
        self.tab = self._load_galtable(self.galpath)
        self.mstar_z, self.mstar_m = self._load_mstar(mstarpath)

    def _load_galtable(self, fname):
        """
        # z, Age, Age_sed, r_des, i_des, r_chft, i_chft
        # for a passively evolving stellar population,
        # BC03 with solar metallicity, no dust, exponentially declining SFH with tau=0.1, Age =10 Gyr at z=0
        """
        _tab = np.genfromtxt(fname)
        columns = ["z", "Age", "Age_sed", "rDES", "iDES", "rCFHT", "iCFHT"]
        tab = pd.DataFrame(data=_tab, columns=columns)
        return tab

    def _load_mstar(self, fname):
        mstar = fio.read(fname)
        mstar_z = [m[0] for m in mstar]
        mstar_m = [m[1] for m in mstar]
        return mstar_z, mstar_m

    def get_mag(self, z, band):
        """
        # magnitudes of a red galaxy template as a function of redshift
        # allows conversion of colors between DES r,i and CFHT r,i and also luminosity evolution of an ageing stellar population
        """
        return np.interp(z, self.tab["z"], self.tab[band])

    def get_mstar(self, z, band):
        """
        # m* magnitude in any of these four filters
        """
        return np.interp(z, self.mstar_z, self.mstar_m) + self.get_mag(z, band) - self.get_mag(z, "iDES")


class ICLProf(object):
    _z_pivot = 0.25
    _m_pivot = 3e14 # M200m, h_72^-1
    _npatch = 40 # number of JK patches
    _nbins = 266
    _pixel_scale = 0.263
    _xarr = np.arange(10,1000,10)
    def __init__(self, jk_profile_root, procmag, H0=72, Om0=0.27):
        self.jk_profile_root = jk_profile_root
        self.procmag = procmag
        self.H0 = H0
        self.Om0 = Om0
        self.cosmo = cosmology.FlatLambdaCDM(H0=self.H0, Om0=self.Om0)
        self._isfitted = False


    def _read_jk_profile(self, patch, suffix=".txt"):
        fname = self.jk_profile_root + str(patch) + suffix
        prof = np.genfromtxt(fname, delimiter=",")
        xarr = prof[:, 0]
        yarr = prof[:, 3] / (self._pixel_scale**2.)
        return xarr, yarr

    def read_icl_raw_prof(self):
        xarr = np.zeros(self._nbins)
        yarr = np.zeros(self._nbins)
        for i in np.arange(self._npatch):
            _xarr, _yarr = self._read_jk_profile(i)
            xarr += _xarr
            yarr += _yarr

        xres = xarr / self._npatch
        yres = yarr / self._npatch

        return xres, yres

    def smooth_vectors(self, xarr, yarr, smooth_min, smooth):

        smoothmatrix = np.zeros((len(xarr), len(yarr)))
        for i in np.arange(smooth):
            smoothmatrix[-(i + 1), - (i + 1)] = 1.
        for i in np.arange(smooth_min + smooth):
            smoothmatrix[i, i] = 1.

        for i in np.arange(smooth_min + smooth, len(xarr) - smooth):
            smoothmatrix[(i - smooth):(i + 1), i] = xarr[(i - smooth):(i + 1 + smooth)] / np.sum(xarr[(i - smooth):(i + 1 + smooth)])

        xarr = np.dot(xarr, smoothmatrix)
        yarr = np.dot(yarr, smoothmatrix)
        return xarr, yarr

    def icl_raw(self, mass=3e14, z=0.25, band="iDES"):
        """
        # get ICL flux (units in counts per arcsec^2 at ZP 30) for clusters at given mass, redshift, and in given filter

        # now assume the same configuration gets put at a different redshift
        # fix the physical stellar surface density

        # (1) re-scale for angular diameter: D_A^{-2}

        # (2) re-scale for (filter,z)-(iCFHT,0.25) color, which contains luminosity distance
        # if color is large, then object is fainter in filter than in rDES, then flux is smaller in filter than in rDES

        # (3) re-scale for mass
        # simple assumption ~valid from Yuanyuan's paper: it's all the same if you look at r in r500 units
        # also blatantly ignoring concentration
        """

        xarr, yarr = self.read_icl_raw_prof()

        cosmo_factor = (self.cosmo.angular_diameter_distance(z).value / self.cosmo.angular_diameter_distance(self._z_pivot).value)**2

        color_factor = self.procmag.get_mag(z, band) - self.procmag.get_mag(self._z_pivot, "iCFHT")
        color_factor = 10**(-0.4 * color_factor)

        mass_factor = (mass / self._m_pivot)**(1. / 3.)

        xarr *= mass_factor
        yarr *= cosmo_factor * color_factor

        return xarr, yarr
    #
    # def icl_smooth(self, mass=3e14, z=0.25, band="iDES", smooth_min=20, smooth=2):
    #     xarr, yarr = self.icl_raw(mass=mass, z=z, band=band)
    #     xarr, yarr = self.smooth_vectors(xarr, yarr, smooth_min, smooth)
    #     return xarr, yarr

    def icl_xgrid(self, xarr, mass=3e14, z=0.25, band="iDES"):
        xin, yin = self.icl_raw(mass, z, band)
        return np.interp(xarr, xin, yin)

    # def icl_xgrid_smooth(self, xarr, mass=3e14, z=0.25, band="iDES", smooth_min=27, smooth=10):
    #     xin, yin = self.icl_smooth(mass, z, band, smooth_min, smooth)
    #     return np.interp(xarr, xin, yin)

    def icl_pixes(self, pixarr, mass=3e14, z=0.25, band="iDES"):
        """X and Y in pix and flux / pix"""
        angarr = pixarr * self._pixel_scale

        # print(self.cosmo.kpc_proper_per_arcmin(0.3).to(u.kpc / u.arcsec))
        arcesc_to_kpc = (self.cosmo.kpc_proper_per_arcmin(0.3).to(u.kpc / u.arcsec)).value
        kpcarr = angarr * arcesc_to_kpc
        # print(kpcarr)
        yarr = self.icl_xgrid(kpcarr, mass, z, band)
        yarr *= self._pixel_scale ** 2.
        return yarr

    def fit_model(self, pixarr, mass=3e14, z=0.25, band="iDES", xpivot=1000):
        self.xpivot = xpivot
        self.mass = mass
        self.z = z
        self.band = band

        yarr = self.icl_pixes(pixarr, mass, z, band)
        self.params = fit_icl_model(pixarr, yarr, xpivot=xpivot)
        self._isfitted = True

    def model(self, xarray, xmin=5):
        if not self._isfitted:
            raise KeyError("Model is not fitted")
        amp = self.params[0]
        alpha = self.params[1]
        beta = self.params[2]

        inds = xarray > xmin
        vals = model_func(xarray[inds], amp=amp, alpha=alpha, beta=beta, xpivot=self.xpivot)
        result = np.zeros(len(xarray))
        result[inds] = vals
        result[np.invert(inds)] = vals.max()

        return result


def fit_icl_model(xarr, yarr, x0=(1, -2.3, -0.9), xpivot=100):
    params = optimize.leastsq(residual, x0, args=(xarr, yarr, xpivot))[0]
    return params

def model_func(xarr, amp=1.16, alpha=-2.27, beta=-1.3, xpivot=100):
    yarr = float(amp) * ((xarr / xpivot)**float(alpha) + (xarr / xpivot)**float(beta))
    return yarr

def residual(x, xarr, yarr, xpivot=100):
    amp = x[0]
    alpha = x[1]
    beta = x[2]

    vals = model_func(xarr, amp, alpha, beta, xpivot)
    res = (yarr - vals) / np.sqrt(np.abs(yarr))
    return res

















