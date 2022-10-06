
import numpy as np
import subprocess
import fitsio as fio
import galsim

from . import shear
from . import render



class Frame(object):
    def __init__(self, catalog, name="canvas", canvas_size=5000, band="i", noise_std=8.36, config_se="config.sex", pixel_scale=0.264):
        self.catalog = catalog
        self.canvas_size = canvas_size
        self.band = band
        self.noise_std = noise_std
        self.name = name
        self.config_se = config_se
        self.pixel_scale = pixel_scale


    def render(self, nprocess=100):
        self.df = render.DrawField(self.canvas_size, self.catalog, band=self.band)
        self.df.prepare()
        self.df.make_wcs()
        self.df.make_infodicts()
        self.df.multi_render(nprocess)
        self.df.collate_stamps()
        self.canvas = self.df.canvas

        self.noise = np.random.normal(scale=self.noise_std, size=(self.canvas_size, self.canvas_size))
        self.canvas += self.noise

        self.write()

    def write(self):
        self.file_name = self.name + ".fits"
        self.canvas.write(self.file_name, clobber=True)
        #fio.write(self.file_name, self.canvas.array, clobber=True)
        self.file_name_psf = self.name + "_epsf.fits"
        #fio.write(self.file_name_psf, self.df.image_epsf.array, clobber=True)
        self.df.image_epsf.write(self.file_name_psf, clobber=True)

    def extract(self):

        self.file_name = self.name + ".fits"
        self.catalog_name = self.name + "_cat.fits"
        self.seg_name = self.name + "_seg.fits"
        #self.weight_name = self.name + "_wmap.fits"

        cmd = "sex " + self.name + ".fits -c " + self.config_se + " -CATALOG_NAME " + self.catalog_name + " -CHECKIMAGE_NAME " + self.seg_name
        print(cmd)
        pp = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = pp.communicate()

        self.scat = fio.read(self.catalog_name)
        self.seg = fio.read(self.seg_name)

    def ksb(self):
        self.file_name = self.name + ".fits"
        self.catalog_name = self.name + "_cat.fits"
        self.seg_name = self.name + "_seg.fits"

        self.scat = fio.read(self.catalog_name)
        self.seg = fio.read(self.seg_name)
        self.file_name_psf = self.name + "_epsf.fits"
        self.epsf = fio.read(self.file_name_psf)

        canvas = fio.read(self.file_name)
        self.canvas = galsim.ImageF(canvas)

        self.ids = self.scat['NUMBER']
        self.cens = np.vstack((self.scat['X_IMAGE'], self.scat['Y_IMAGE'])).T
        self.sizes = np.around(self.scat["FWHM_IMAGE"] * 2)
        self.sizes = np.max((np.zeros(len(self.sizes)) + 8, self.sizes), axis=0)
        self.sc = shear.Shear(self.canvas, self.epsf, self.seg, self.pixel_scale)
        self.sc.extract_stamps(self.cens, imasks=self.ids, sizes=self.sizes)
        self.sc.estimate_shear(sky_var=self.noise_std**2)