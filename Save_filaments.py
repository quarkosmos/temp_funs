"""
Save_Filaments.py is to save each filament in a fits file to use it as an input of FilFinder

The whole FUNS/SerB mapping region will be used without cropping around each filament

The velocity range will be [3, 13] km/s.
The Cubes of filaments will be generated.

"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
import os
from astropy.wcs import WCS
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits

import funstools
from funstools import Cube2map
from run_funstools_Serp import get_filename
from run_funstools_Serp import get_maxmin_maps

from astropy import constants as co
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.visualization import simple_norm

from spectral_cube import SpectralCube
import pyspeckit

import pickle
from astropy.modeling.models import Gaussian1D
from astropy.table import Table, Column
from astropy.io.ascii import read, write
import copy



class SaveFil(Cube2map):

    def __init__(self, pickled_file=None, ff_ver=None, ffpath=None, fntype='finfn'):

        if pickled_file is None:
            pickledir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/', 'Filfinder', 'cube_pickle')
            pickled_file = pickledir+'/SerB_C18O_v06_v0313_snr2.0.pickle'

        if os.path.exists(pickled_file):
            with open(file=pickled_file, mode='rb') as fr:
                self._cube = pickle.load(fr)

        if ff_ver is not None:
            self._ff_ver = ff_ver
        else:
            self._ff_ver = '01252022'

        if ffpath is not None:
            self._ffpath = ffpath
        else:
            self._ffpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/Filfinder/'+self._ff_ver+'/SeqFf_rev/'

        if fntype is not None:
            self._fntype = fntype
        else:
            self._fntype = 'finfn'

        fftable_file = self._ffpath+'SeqFF_C18O_rev5.dat'
        self._ft = read(fftable_file)




    def fil_com_v(self, tp, vc, std):
        gcm = Gaussian1D(tp, vc, std)
        z = gcm(self._cube.x)
        return z


    def reform_data(self, fnum=None):
        ccube = copy.copy(self._cube)
        nz = np.shape(ccube.data)[0]
        not_a_member = np.zeros(nz)
        comp = self._ft[self._ft[self._fntype]==fnum]
        for rr in range(ccube.nr):
            for dd in range(ccube.nd):
                ct = comp[(comp['rp'] == rr) & (comp['dp'] == dd)]
                if len(ct) == 0:
                    cy = not_a_member
                elif len(ct) == 1:
                    cy = self.fil_com_v(ct['tp'][0], ct['vp'][0], ct['sd'][0])
                ccube.data[:, dd, rr] = cy
        return ccube.data





    def write2fits_cube(self, fnum=None):
        if fnum is not None:
            newdata = self.reform_data(fnum)
        else:
            import warnings
            warnings.warn("fnum is not given")

        ccube = copy.copy(self._cube)
        h = ccube.header
        h['FILNUM'] = fnum
        h['COMMENT'] = 'Saved by SaveFil with SeqFF Filaments'
        wcubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/Filfinder/'+self._ff_ver+'/SaveFil_fits/3D/'
        wcubename = 'SerB_C18O_{}_fn_{}.fits'.format(self._ff_ver, fnum)
        fits.writeto(wcubedir+wcubename, data = newdata, header = h, overwrite=True)



    def write2pickle(self, fnum=None):
        if fnum is not None:
            newdata = self.reform_data(fnum)
        else:
            import warnings
            warnings.warn("fnum is not given")

        wcubepickle = '/Users/khkim/Dropbox/Work5/FUNS-Serp/Filfinder/'+self._ff_ver+'/SaveFil_pickle/'
        wcubename = 'SerB_C18O_{}_fn_{}.pickle'.format(self._ff_ver, fnum)
        with open(file=wcubepickle+wcubename, mode='wb') as fa:
            pickle.dump(newdata, fa)



    def write2fits_m0(self, fnum=None):
        data2d = np.zeros((1, self._cube.nd, self._cube.nr))
        comp = self._ft[self._ft[self._fntype]==fnum]
        for rr in range(self._cube.nr):
            for dd in range(self._cube.nd):
                ct = comp[(comp['rp'] == rr) & (comp['dp'] == dd)]
                if len(ct) == 0:
                    area = float('nan')#0
                elif len(ct) == 1:
                    area = ct['area'][0]
                data2d[0, dd, rr] = area

        h = self._cube.header2d
        h['FILNUM'] = fnum
        h['COMMENT'] = 'Saved by SaveFil with SeqFF Filaments'
        h['COMMENT'] = 'm0 (area of a gaussian component)'
        wcubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/Filfinder/'+self._ff_ver+'/SaveFil_fits/2D/'
        wcubename = 'SerB_C18O_{}_fn_{}_m0.fits'.format(self._ff_ver, fnum)
        fits.writeto(wcubedir+wcubename, data = data2d, header = h, overwrite=True)


    #
    # def plot2D():
    #
    # def plot3D():
