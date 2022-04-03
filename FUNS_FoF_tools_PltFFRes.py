
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io.ascii import read, write
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.wcs import WCS
from astropy.table import Table
from mpl_toolkits.mplot3d import axes3d

from mpl_toolkits.mplot3d.axes3d import Axes3D

from astropy.modeling.models import Gaussian1D
from spectral_cube import SpectralCube
from astropy import units as u
import os


import funstools
from funstools import get_rms
#from funstools import Cube2map
from funstools import save_fits

from SerB_Cube import SerBCube2map as Cube2map
from run_funstools_Serp import get_crop_filename
from matplotlib.backends.backend_pdf import PdfPages

class FUNS_FoF:

    """
    FUNS_FoF is the collection of functions which are useful for running Friends-of-Friends and examining the resulted filaments structures.
    The functions were originally written by Hyunju Yu.

    Kyoung Hee is trying to collect the fuctions in a class, so it can be used in some convenent way...

    run_FoF:
    FoF_line_profiles
    FoF_3D_plot
    Fil_in_2D

    """

    def __init__(self, decompose_key='', decompose_path='', decompose_result='', region='', line='C18O', vtype='v06', rmssize=None, snr=None, ssm_find=None, vsm_find=None, ssm_fit=None, vsm_fit=None, sigma=None, v0=None, rd0=None, fof_wfile=''):
        #region, line, vtype, snr, ssm_find, vsm_find, ssm_fit, vsm_fit, mww, sigma, v0, rd0 = input("Enter parameters: region line vtype snr ssm_find vsm_find ssm_fit vsm_fit fof_sigma_cut dist_v0 dist_rd0:").split()
        #EX# SerB_A_ C18O v06 3 3 3 1 1 mww062 4.0 0.06 2.0
        self._region = region
        self._snr = float(snr)
        self._ssm_find = int(ssm_find)
        self._vsm_find = int(vsm_find)
        self._ssm_fit = int(ssm_fit)
        self._vsm_fit = int(vsm_fit)
        self._sigma = float(sigma)
        self._dist_v0 = float(v0)
        self._dist_rd0 = float(rd0)
        self._line = line
        self._vtype = vtype
        self._rmssize = int(rmssize)


        if self._region=='SerB_':
            self._region_path = 'Main'
        elif self._region == 'SerB_A_':
            self._region_path = 'A'
        elif self._region == 'SerB_B_':
            self._region_path = 'B'



        # self._decompose_key = str(ssm_find)+str(vsm_find)+str(ssm_fit)+str(vsm_fit)+vtype+'snr'+str(self._snr)
        # self._decompose_path = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose2/'+self._decompose_key+'/' #changed the path decompose2
        # self._decompose_result = self._decompose_path+'SerB_'+self._line+'_decompose_result_'+self._decompose_key+'.dat'
        # self._FoF_save='/Users/khkim/Dropbox/Work5/FUNS-Serp/FOF/FoF_dc2/'
        # self._fof_wfile = self._FoF_save+'FoF_result_'+self._line+'_'+self._decompose_key+'_'+str(self._sigma)+'sig_dv0_'+str(self._dist_v0)+'_rd0_'+str(self._dist_rd0)+'.dat'
        # self._plot_note = self._line+'/ '+self._decompose_key+'/ '+str(self._sigma)+'sig_dv0_'+str(self._dist_v0)+'_rd0_'+str(self._dist_rd0)

        #self._decompose_key = str(ssm_find)+str(vsm_find)+str(ssm_fit)+str(vsm_fit)+vtype+'snr'+str(self._snr)+mww
        #self._decompose_path = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose5_AB/'+self._decompose_key+'/' #changed the path decompose2
        #self._decompose_result = self._decompose_path+self._region+self._line+'_decompose_result_'+self._decompose_key+'.dat'
        #self._FoF_save='/Users/khkim/Dropbox/Work5/FUNS-Serp/FOF/FoF_dc5_AB/'
        #self._fof_wfile = self._FoF_save+'FoF_result_'+self._region+self._line+'_'+self._decompose_key+'_'+str(self._sigma)+'sig_dv0_'+str(self._dist_v0)+'_rd0_'+str(self._dist_rd0)+'.dat'


        if decompose_key is not None:
            self._decompose_key = decompose_key
        else:
            self._decompose_key = 'sn2.5d3s1'

        if decompose_path is not None:
            self._decompose_path = decompose_path
        else:
            self._decompose_path = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+self._decompose_key+'/output/gpy_decomposed/'

        if decompose_result is not None:
            self._decompose_result = decompose_result
        else:
            self._decompose_result = self._decompose_path+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(self._decompose_key)


        self._FoF_save=os.path.join(self._decompose_path, 'fof')
        if fof_wfile is not None:
            self._fof_wfile = fof_wfile
        else:
            self._fof_wfile = self._FoF_save+'/FoF_result_'+self._region+self._line+'_'+self._decompose_key+'_'+str(self._sigma)+'sig_dv0_'+str(self._dist_v0)+'_rd0_'+str(self._dist_rd0)+'.dat'

        self._plot_note = self._region+self._line+'/ '+self._decompose_key+'/ '+str(self._sigma)+'sig_dv0_'+str(self._dist_v0)+'_rd0_'+str(self._dist_rd0)


        if not os.path.exists(self._decompose_path):
            raise TypeError('No decomposed result table exist: Please run decompose first.')
        else:
            pass

    _tdet = None


############################
# run FoF by Dr. Yu's code
############################

    def run_fof(self):

        result = read(self._decompose_result)

        sigma = self._sigma
        dist_v0 = self._dist_v0
        dist_rd0 = self._dist_rd0


        self._tdet = result[result['snr'] > sigma]
        tdet = self._tdet
        tdet['fn'] = 0
        tdet['in'] = 0
        #tdet['color'] = 0


        ifil = 1
        while np.sum(tdet['fn'] == 0) > 0:
            remain = tdet[tdet['fn'] == 0]
            remain_max = np.max(remain['tp'])
            print('n_remain = {}'.format(len(remain)))
            kseed = np.argwhere((tdet['fn'] == 0) & (tdet['tp'] == remain_max))[0]
            ir, id, iv = tdet['rp'][kseed], tdet['dp'][kseed], tdet['vp'][kseed]
            print('kseed = {}'.format(kseed))
            print('kseed info = {}, {}, {}'.format(ir, id, iv))
            k = 1
            tdet['fn'][kseed] = ifil
            tdet['in'][kseed] = k
            nfc = 1
            while nfc > 0:
                nfc = 0
                k += 1
                if k == 2 :
                    fcond = ((tdet['rp'] - ir)**2 + (tdet['dp'] - id)**2 <= dist_rd0**2) \
                            & ~((tdet['rp'] == ir) & (tdet['dp'] == id)) & (np.abs(tdet['vp'] - iv) <= dist_v0)
                    nfc += np.sum(fcond)
                    tdet['fn'][fcond] = ifil
                    tdet['in'][fcond] = k
                else :
                    seeds = tdet[(tdet['fn'] == ifil) & (tdet['in'] == k-1)]
                    for j in range(len(seeds)):
                        ir, id, iv = seeds['rp'][j], seeds['dp'][j], seeds['vp'][j]
                        fcond = (tdet['fn'] == 0) & ((tdet['rp'] - ir) ** 2 + (tdet['dp'] - id) ** 2 <= dist_rd0 ** 2) \
                                & ~((tdet['rp'] == ir) & (tdet['dp'] == id)) & (np.abs(tdet['vp'] - iv) <= dist_v0)
                        nfc += np.sum(fcond)
                        tdet['fn'][fcond] = ifil
                        tdet['in'][fcond] = k
                print('nfc at {} fil {} loop = {}'.format(ifil, k, nfc))
                # print kseed, ifil, k, nfc
            ifil += 1

#        self._fof_wfile = FoF_save+'FoF_result_'+self._line+'_'+self._decompose_key+'_'+str(float(sigma))+'sig_dv0_'+str(dist_v0)+'_rd0_'+str(dist_rd0)+'.dat'
        write(tdet, self._fof_wfile)#, overwrite=True)
        self._tdet = tdet
        return self._tdet

############################
# search FoF results exist, if not run FoF for the requested specific conditions
############################
    def get_tdet(self):
         if self._tdet is None:
             if os.path.exists(self._fof_wfile):
                 self._tdet = read(self._fof_wfile)
             else:
                 self._tdet = self.run_FoF()
         return self._tdet

############################
# 3D plot for FoF components by Dr. Yu's code
############################
## Dr. Yu’s 3D plots code
    def fof_3D_plot(self, basic_m0=None, tdet=None, nbeam=None):
        tdet = self.get_tdet()

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #basic_m0 = fits.open('/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v-520/sn3vs2ss2/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits')
        #hdu = basic_m0[0]
        #w_line = WCS(hdu.header)

        #tdet = read('/Users/Hyunju/Dropbox/FUNS_OrionB/Individual_Filaments/2step_modify/Merged_FoF_result_2nd_mod.dat')

        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        #line, vtype, snrchoice, ssmo, vsmo, v1, v2 = input("Enter parameters: 'line' vtype snr ssmo vsmo v1 v2:").split()
        # ex) Enter parameters: 'line' snr vtype ssmo vsmo v1 v2:N2HP v10 3 2 2 3 13
        #v1 = float(v1)
        #v2 = float(v2)
        #snrchoice = float(snrchoice)
        v1, v2 = input("Enter velocity range: v1 v2:").split()

        #ssmo = int(ssmo)
        #vsmo = int(vsmo)

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize
        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)

        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        # if not os.file.exists(m0filename):
        #     save_fits(m0filename, cube.m, cube.header2d, overwrite=True)
        # else:
        #     pass

        #m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
        #w_line = basic_m0.wcs2d




        #onebeam = 5.171479 # (47 × 47 × π) ÷ (4 ln(2) × 22 × 22) for OriB
        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        nfil = np.max(tdet['fn'])

        f_count = {}
        for lst in tdet['fn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1

        print(f_count)


        fig = plt.figure(figsize=(5/xyr, 5/xyr))
        ax = fig.add_subplot(111, projection='3d')  # Axe3D object
        ax.set_xlim3d(np.min(tdet['rp'])-20, np.max(tdet['rp'])+20)
        ax.set_ylim3d(np.min(tdet['dp'])-20, np.max(tdet['dp'])+20)
        ax.set_zlim3d(v1, v2)

        # line_image = basic_m0[0]
        # iy, ix = np.mgrid[0:line_image.shape[1], 0:line_image.shape[2]]

        cindex = 0
        a = 0
        for fil, cnt in f_count.items():

            if cnt > npix0:

                cindex += 1

                print(fil, cnt)
                fof_result = tdet[tdet['fn'] == fil]

                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.
                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # central scatter
                ax.scatter(fof_result['rp'], fof_result['dp'], fof_result['vp'], marker='.', s=10, color = clr_list[cindex], alpha=0.3)


                xflat = np.full_like(fof_result['rp'], min(ax.get_xlim()))
                yflat = np.full_like(fof_result['dp'], max(ax.get_ylim()))
                zflat = np.full_like(fof_result['vp'], min(ax.get_zlim()))

                ax.set_xlabel('RA')
                ax.set_ylabel('Dec.')
                ax.set_zlabel('Velocity')

                a += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'num. of filaments = '+str(cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)
        plt.show()
        # fig.savefig('/Users/Hyunju/Dropbox/FUNS_OrionB/Plot/HFoF_3d_view.pdf', dpi=300, max_dpi=500)


#######################################################################
## Dr. Yu’s code for mapping the filaments in 2D (RA & Dec) plane (she used > 5 beams )
 #######################################################################

    def fof_2D_plot(self, basic_m0=None, tdet=None, nbeam=None, note='', overplot_m0=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
    'olive', 'peru', 'gold', 'coral', 'navy',\
    'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
    'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
   'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
   'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
   'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
   'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
   'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
   'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        #v1, v2 = input("Enter velocity range: v1 v2:").split()
        v1='2'
        v2='14'

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'

        rs = self._rmssize
        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        fig1 = plt.figure(figsize=(8, 8/xyr))
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['fn']!=0]

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        coloridx = np.zeros(len(tdet['fn']))
        nfil = np.max(tdet['fn'])
        print('nfil =', nfil)

        f_count = {}
        for lst in tdet['fn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1



        f_count_sort = sorted(f_count.items())

        #print(f_count)

        ax1.set_xlabel('R.A.', fontsize=15)
        ax1.set_ylabel('Dec.', fontsize=15)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=15)
        ax1.set_xlim(-10, cube.nr+10)
        ax1.set_ylim(-10, cube.nd+10)
        ax1.set_xticks(np.arange(-10, cube.nr+10, 10))
        ax1.set_yticks(np.arange(-10, cube.nd+10, 10))

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
            cindex = ifil+1
            if cnt > npix0 :
                print(fil, cnt)
                fof_result = tdet[tdet['fn'] == fil]
                #value_index = tdet['fn'].index(fil)
                ctemp_index = np.where(tdet['fn']==fil)
                coloridx[ctemp_index] = cindex
                cindex += 1

                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # for c in range(0, len(tdet)):
                #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])



                a += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FoF: num. of filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)

        plt.show()
        dirname = str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'mww'
        plt.savefig(self._fof_wfile+'_nb{}.pdf'.format(nbeam))
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)


############################
# statistics of decompose results and FoF results
############################
    def dcstat(self):
        dr_data = read(self._decompose_result)
        dr_rp = dr_data[0][:]
        dr_dp = dr_data[1][:]
        #drst =  Table(names=('rp', 'dp', 'tn'), dtype=('i4', 'i4', 'i4'))
        statfile = open(self._decompose_path+'decompose_stat.txt', 'a')
        arp = np.array([], dtype='i4')
        adp = np.array([], dtype='i4')
        atn = np.array([], dtype='i4')
        for ira in range(118): #ra direction
            for idec in range(235): #dec direction
                w = np.where((dr_rp == ira) & (dr_dp == idec))
                if np.size(w)!=0:
                    arp = np.append(arp, ira)
                    adp = np.append(adp, idec)
                    atn = np.append(atn, np.size(w))
                else:
                    pass
        max_num_dc = np.max(atn)
        npixel = np.size(arp)
        #print(max_num_dc)
        statfile.write(self._decompose_result)
        statfile.write("\n")
        statfile.write("the total number of pixels with decomposed components: "+str(npixel))
        statfile.write("\n")
        statfile.write("the maximum number of components: "+str(max_num_dc))
        statfile.write("\n")
        statfile.write("N(compentns),    N(pixels with the components),    rate of pixels with the num. of components")
        nsum = 0
        for i in range(max_num_dc):
            c_i = i+1
            nc_i = np.where(atn == c_i)
            nc = np.size(nc_i)
            rcrate = nc / npixel
            print(c_i, nc)
            statfile.write("\n")
            statfile.write(str(c_i)+'    '+str(nc)+'    '+format(rcrate, '.3f'))
            nsum = nsum + nc
            #statfile.write("the num. of pixels with number of components of "+str(c_i)+": "+str(nc)+" ("+str(ncrate)+")")
        statfile.write("\n")
        statfile.write('total'+'    '+str(nsum)+'    '+format(nsum/npixel, '.3f'))
        statfile.write("\n")
        statfile.close()



############################
# Save each filament in a
############################
    def save_fil_by_fil(self, tdet=None, nbeam=None):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        tdet = self.get_tdet()

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        f_count = {}
        for lst in tdet['fn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1

        ff_save_dic = {}
        large_fil_list = [k for k, v in f_count.items() if v > npix0]
        for i in range(np.size(large_fil_list)):
            ifnum = large_fil_list[i]
            fil_ord_num = i+1
            ff_save_dic[fil_ord_num] = tdet[tdet['fn']==ifnum]
            temp_fil_info = ff_save_dic[fil_ord_num]
            temp_fil_info['color'] = clr_list[fil_ord_num]
            ff_save_dic[fil_ord_num] = temp_fil_info

        import pickle
        with open(file=self._fof_wfile+'_nb{}.pickle'.format(nbeam), mode='wb') as f:
            pickle.dump(ff_save_dic, f)


############################
# Save each filament in a
############################
    def SeqFF_save_fil_by_fil(self, tdet=None, nbeam=None):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        tdet = self.get_tdet()

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        f_count = {}
        for lst in tdet['finfn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1

        ff_save_dic = {}
        large_fil_list = [k for k, v in f_count.items() if v > npix0]

        for i in range(np.size(large_fil_list)):
            ifnum = large_fil_list[i]
            ff_save_dic[ifnum] = tdet[tdet['finfn']==ifnum]
            temp_fil_info = ff_save_dic[ifnum]
            temp_fil_info['color'] = clr_list[ifnum]
            ff_save_dic[ifnum] = temp_fil_info

        import pickle
        with open(file=self._fof_wfile+'_nb{}.pickle'.format(nbeam), mode='wb') as f:
            pickle.dump(ff_save_dic, f)



#########################################################
# load pickled filament data which is larger than npix0
#########################################################
    def load_filaments_info(self, nbeam=None):
        import pickle
        pickled_file = self._fof_wfile+'_nb{}.pickle'.format(nbeam)
        if os.path.exists(pickled_file):
            with open(file=pickled_file, mode='rb') as f:
                fil_dic = pickle.load(f)
        else:
            self.SeqFF_save_fil_by_fil(nbeam=None)
            with open(file=pickled_file, mode='rb') as f:
                fil_dic = pickle.load(f)

        return fil_dic






#################################
# Plot Filaments on channel maps
#################################
    def plot_all_fils_on_chmap(self, basic_m0=None, tdet=None, nbeam=None, note=''):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurplue', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateble', 'mediumvioletred', 'chartreuse', 'burlywood']


        t_fil = self.load_filaments_info()
        n_fil = len(t_fil)


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        from multiaxes import Multiaxes


        #v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # set map value_range, colormap
        int_vr = [0., 0.13]
        #cmap = plt.get_cmap('nipy_spectral')
        #cmap.set_bad('grey')
        cmap = plt.get_cmap('Greys')
        cmap.set_bad('white')

        maps, labels = cube.chmap(50, (5.5, 10.5))

        pdfname = self._fof_wfile+'_chmap.pdf'
        with PdfPages(pdfname) as pdf:
            for i in range(5):
                # draw figure using Multiaxes
                mx = Multiaxes(2, 5, 2, xyr, 0.36, 0.56, margin=(0.03, 0.03, 0.75, 0.05), scale=0.7, proj=cube.wcs2d)
                mx.shareaxes((True, True), 0.1)
                fig, ax, _ = mx.drawfig()
                cax = mx.sharecolorbar('right', 0.15, 0.13)
                for j in range(10):
                    iax = ax[j//5, j%5]
                    k = i*10+j
                    cs = iax.imshow(maps[k], cmap=cmap, vmin=int_vr[0], vmax=int_vr[1])
                    iax.annotate(labels[k]+' km/s', xy=(0.55, 0.05), xycoords='axes fraction', \
                                 ha='center', va='baseline', color='black', fontsize='x-small')
                    iax.tick_params(axis='x', labelsize=5)
                    iax.tick_params(axis='y', labelsize=5)

                    if j == 5:
                        iax.coords[0].set_axislabel('R.A. (J2000)', fontsize=7)
                        iax.coords[1].set_axislabel('Dec. (J2000)', fontsize=7)
                    else:
                        iax.coords[0].set_axislabel(' ')
                        iax.coords[1].set_axislabel(' ')

                    for kk in range(n_fil):
                        fof_result = t_fil[kk+1]
                        fil_array = np.zeros((cube.nd, cube.nr))
                        for ii in range(len(fof_result)):
                            fil_array[fof_result['dp'][ii], fof_result['rp'][ii]] = 1.

                        kernel = Gaussian2DKernel(0.4)
                        fil_array = convolve(fil_array, kernel)
                        iax.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=0.5, alpha=0.1, color = fof_result['color'])
                        iax.contour(fil_array, levels=[0.5], alpha=0.5, colors=fof_result['color'])

                plt.colorbar(cs, cax)
                cax.tick_params(labelsize=5)
                cax.set_ylabel(r'Integrated intensity (K km s$^{-1}$)', fontsize=7)

                pdf.savefig()
                plt.close()




##################################################################################
# Plot Filaments which are in each correponding velocity range on channel maps
#################################################################################
    def plot_vc_fils_on_chmap(self, basic_m0=None, tdet=None, nbeam=None, note=''):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurplue', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateble', 'mediumvioletred', 'chartreuse', 'burlywood']


        t_fil = self.load_filaments_info()
        n_fil = len(t_fil)


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        from multiaxes import Multiaxes


        #v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # set map value_range, colormap
        int_vr = [0., 0.13]
        #cmap = plt.get_cmap('nipy_spectral')
        #cmap.set_bad('grey')
        cmap = plt.get_cmap('Greys')
        cmap.set_bad('white')

        maps, labels = cube.chmap(50, (5.5, 10.5))

        pdfname = self._fof_wfile+'_chmap_vp.pdf'
        with PdfPages(pdfname) as pdf:
            for i in range(5):
                # draw figure using Multiaxes
                mx = Multiaxes(2, 5, 2, xyr, 0.36, 0.56, margin=(0.03, 0.03, 0.75, 0.03), scale=0.7, proj=cube.wcs2d)
                mx.shareaxes((True, True), 0.1)
                fig, ax, _ = mx.drawfig()
                cax = mx.sharecolorbar('right', 0.15, 0.1)
                for j in range(10):
                    iax = ax[j//5, j%5]
                    k = i*10+j
                    cs = iax.imshow(maps[k], cmap=cmap, vmin=int_vr[0], vmax=int_vr[1])
                    iax.annotate(labels[k]+' km/s', xy=(0.55, 0.05), xycoords='axes fraction', \
                                 ha='center', va='baseline', color='black', fontsize='x-small')
                    iax.tick_params(axis='x', labelsize=5)
                    iax.tick_params(axis='y', labelsize=5)

                    if j == 5:
                        iax.coords[0].set_axislabel('R.A. (J2000)', fontsize=7)
                        iax.coords[1].set_axislabel('Dec. (J2000)', fontsize=7)
                    else:
                        iax.coords[0].set_axislabel(' ')
                        iax.coords[1].set_axislabel(' ')

                    for kk in range(n_fil):
                        fof_result = t_fil[kk+1]
                        fvp = fof_result['vp']
                        vp_min = np.quantile(fvp, 0., axis=0)
                        vp_max = np.quantile(fvp, 1., axis=0)
                        ch_vr = [5.5+(k*0.1), 5.5+(k+1)*0.1]
                        if (vp_min <= ch_vr[0]) and (vp_max >= ch_vr[1]):
                            fil_array = np.zeros((cube.nd, cube.nr))
                            for ii in range(len(fof_result)):
                                fil_array[fof_result['dp'][ii], fof_result['rp'][ii]] = 1.

                            kernel = Gaussian2DKernel(0.4)
                            fil_array = convolve(fil_array, kernel)
                            iax.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=0.5, alpha=0.1, color = fof_result['color'])
                            iax.contour(fil_array, levels=[0.5], alpha=0.5, colors=fof_result['color'])

                plt.colorbar(cs, cax)
                cax.tick_params(labelsize=5)
                cax.set_ylabel(r'Integrated intensity (K km s$^{-1}$)', fontsize=7)

                pdf.savefig()
                plt.close()



##################################################################################
# Plot Filaments which are in each correponding velocity range on channel maps
#################################################################################
    def plot_chmap(self, basic_m0=None, tdet=None, nbeam=None, note=''):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurplue', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateble', 'mediumvioletred', 'chartreuse', 'burlywood']


        t_fil = self.load_filaments_info()
        n_fil = len(t_fil)


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        from multiaxes import Multiaxes


        #v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # set map value_range, colormap
        int_vr = [0., 0.3]
        cmap = plt.get_cmap('nipy_spectral')
        cmap.set_bad('grey')
        #cmap = plt.get_cmap('Greys')
        #cmap.set_bad('white')

        #maps, labels = cube.chmap(50, (5.5, 10.5))
        maps, labels = cube.chmap(25, (5.5, 10.5))

        pdfname = self._fof_wfile+'_chmap_only_clr.pdf'
        with PdfPages(pdfname) as pdf:
            for i in range(5):
                # draw figure using Multiaxes
                mx = Multiaxes(2, 5, 2, xyr, 0.36, 0.56, margin=(0.03, 0.03, 0.75, 0.03), scale=0.7, proj=cube.wcs2d)
                mx.shareaxes((True, True), 0.1)
                fig, ax, _ = mx.drawfig()
                cax = mx.sharecolorbar('right', 0.15, 0.1)
                for j in range(10):
                    iax = ax[j//5, j%5]
                    k = i*10+j
                    cs = iax.imshow(maps[k], cmap=cmap, vmin=int_vr[0], vmax=int_vr[1])
                    iax.annotate(labels[k]+' km/s', xy=(0.55, 0.05), xycoords='axes fraction', \
                                 ha='center', va='baseline', color='black', fontsize='x-small')
                    iax.tick_params(axis='x', labelsize=5)
                    iax.tick_params(axis='y', labelsize=5)

                    if j == 5:
                        iax.coords[0].set_axislabel('R.A. (J2000)', fontsize=7)
                        iax.coords[1].set_axislabel('Dec. (J2000)', fontsize=7)
                    else:
                        iax.coords[0].set_axislabel(' ')
                        iax.coords[1].set_axislabel(' ')


                plt.colorbar(cs, cax)
                cax.tick_params(labelsize=5)
                cax.set_ylabel(r'Integrated intensity (K km s$^{-1}$)', fontsize=7)

                pdf.savefig()
                plt.close()


##################################################################################
# Plot all decomposed components in a PPV space
#################################################################################

    def dc_3D_plot(self, basic_m0=None, tdet=None, nbeam=None):
        tdet = self.get_tdet()

        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize
        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)

        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)

        fig = plt.figure(figsize=(5/xyr, 5/xyr))
        ax = fig.add_subplot(111, projection='3d')  # Axe3D object
        ax.set_xlim3d(np.min(tdet['rp'])-20, np.max(tdet['rp'])+20)
        ax.set_ylim3d(np.min(tdet['dp'])-20, np.max(tdet['dp'])+20)
        ax.set_zlim3d(v1, v2)

        ax.scatter(tdet['rp'], tdet['dp'], tdet['vp'], marker='.', s=10, color = 'gray', alpha=0.3)
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec.')
        ax.set_zlabel('Velocity')


        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'all decomposed components', fontsize=11)
        plt.show()
        # fig.savefig('/Users/Hyunju/Dropbox/FUNS_OrionB/Plot/HFoF_3d_view.pdf', dpi=300, max_dpi=500)




#######################################################################
## To plot seeds from FIVE step 1
 #######################################################################

    def five_step1_seed_plot(self, basic_m0=None, tdet=None, note='', nbeam=None, overplot_m0=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['black', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        fig1 = plt.figure(figsize=(8, 8/xyr))
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['seed']!=0]

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam


        coloridx = np.zeros(len(tdet['seed']))
        nseed = np.max(tdet['seed'])
        print('seed =', nseed)

        f_count = {}
        for lst in tdet['seed']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1

        print(f_count)

        a = 0
        ax1.set_xlabel('R.A.', fontsize=15)
        ax1.set_ylabel('Dec.', fontsize=15)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=15)
        ax1.set_xlim(-10, cube.nr+10)
        ax1.set_ylim(-10, cube.nd+10)
        ax1.set_xticks(np.arange(-10, cube.nr+10, 10))
        ax1.set_yticks(np.arange(-10, cube.nd+10, 10))
#
        cindex = 0
        for seed, cnt in f_count.items():
            if cnt > npix0 :
                print(seed, cnt)
                fof_result = tdet[tdet['seed'] == seed]
                #value_index = tdet['fn'].index(fil)
                ctemp_index = np.where(tdet['seed']==seed)
                coloridx[ctemp_index] = cindex


                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # for c in range(0, len(tdet)):
                #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                print('color index = {}, {}'.format(cindex, clr_list[cindex]))
                ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
            cindex += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FIVE step 1: num. of seeds = '+str(cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)


        plt.show()
        plt.savefig(self._fof_wfile+'_five_step1_seeds_nb{}.pdf'.format(nbeam))
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)

##########################################################################
###### To plot 2D plot for filaments from Step 1 in FIVE algorightm ######
##########################################################################


    def five_step1_fils_plot(self, basic_m0=None, tdet=None, note='', nbeam=None, overplot_m0=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['black', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        fig1 = plt.figure(figsize=(8, 8/xyr))
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['fn']!=0]

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam


        coloridx = np.zeros(len(tdet['fn']))
        nseed = np.max(tdet['fn'])
        print('fils =', nseed)

        f_count = {}
        for lst in tdet['fn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1


        f_count_sort = sorted(f_count.items())

        #print(f_count)

        ax1.set_xlabel('R.A.', fontsize=15)
        ax1.set_ylabel('Dec.', fontsize=15)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=15)
        ax1.set_xlim(-10, cube.nr+10)
        ax1.set_ylim(-10, cube.nd+10)
        ax1.set_xticks(np.arange(-10, cube.nr+10, 10))
        ax1.set_yticks(np.arange(-10, cube.nd+10, 10))

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
# #
            cindex = ifil+1
#         for fil, cnt in f_count.items():

            if cnt > npix0 :
                print(fil, cnt)
                fof_result = tdet[tdet['fn'] == fil]
                #value_index = tdet['fn'].index(fil)
                ctemp_index = np.where(tdet['fn']==fil)
                coloridx[ctemp_index] = cindex


                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # for c in range(0, len(tdet)):
                #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                print('color index = {}, {}'.format(cindex, clr_list[cindex]))
                ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
            cindex += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FIVE: num. of step 1 filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)

        plt.show()
        plt.savefig(self._fof_wfile+'_five_step1_fils_nb{}.pdf'.format(nbeam))
        #plt.savefig(self._fof_wfile+'_five_step1_fils.pdf')
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)


##########################################################################
###### To plot 2D plot for filaments from Step 2 in FIVE algorightm ######
##########################################################################


    def five_step2_fils_plot(self, basic_m0=None, tdet=None, note='', nbeam=None, overplot_m0=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['black', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        fig1 = plt.figure(figsize=(8, 8/xyr))
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['fn2']!=0]

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam


        coloridx = np.zeros(len(tdet['fn2']))
        nseed = np.max(tdet['fn2'])
        print('fils =', nseed)

        f_count = {}
        for lst in tdet['fn2']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1


        f_count_sort = sorted(f_count.items())

        #print(f_count)

        a = 0
        ax1.set_xlabel('R.A.', fontsize=15)
        ax1.set_ylabel('Dec.', fontsize=15)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=15)
        ax1.set_xlim(-10, cube.nr+10)
        ax1.set_ylim(-10, cube.nd+10)
        ax1.set_xticks(np.arange(-10, cube.nr+10, 10))
        ax1.set_yticks(np.arange(-10, cube.nd+10, 10))

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            #if f_count_sort[fil][0] == selfil[isf]:
            #    cnt = f_count_sort[fil][1]
            #    cindex = fil+1
        #for fil, cnt in f_count.items():
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
#
#             if cnt > npix0:
#
#
#
#
# #
            cindex = ifil+1
#         for fil, cnt in f_count.items():

            if cnt > npix0 :
                print(fil, cnt)
                fof_result = tdet[tdet['fn2'] == fil]
                #value_index = tdet['fn'].index(fil)
                ctemp_index = np.where(tdet['fn2']==fil)
                coloridx[ctemp_index] = cindex


                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # for c in range(0, len(tdet)):
                #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                print('color index = {}, {}'.format(cindex, clr_list[cindex]))
                ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                a+=1


        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FIVE: num. of step 2 filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)

        plt.show()
        plt.savefig(self._fof_wfile+'_five_step2_fils_nb{}.pdf'.format(nbeam))
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)


##########################################################################
###### To plot 2D plot for filaments from Step 3 in FIVE algorightm ######
##########################################################################


    def five_step3_fils_plot(self, basic_m0=None, tdet=None, note='', nbeam=None, overplot_m0=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['black', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        fig1 = plt.figure(figsize=(8, 8/xyr))
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['fn3']!=0]

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam


        coloridx = np.zeros(len(tdet['fn3']))
        nseed = np.max(tdet['fn3'])
        print('fils =', nseed)

        f_count = {}
        for lst in tdet['fn3']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1


        #print(f_count)
        f_count_sort = sorted(f_count.items())

        #print(f_count)

        a = 0
        ax1.set_xlabel('R.A.', fontsize=15)
        ax1.set_ylabel('Dec.', fontsize=15)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=15)
        ax1.set_xlim(-10, cube.nr+10)
        ax1.set_ylim(-10, cube.nd+10)
        ax1.set_xticks(np.arange(-10, cube.nr+10, 10))
        ax1.set_yticks(np.arange(-10, cube.nd+10, 10))

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            #if f_count_sort[fil][0] == selfil[isf]:
            #    cnt = f_count_sort[fil][1]
            #    cindex = fil+1
        #for fil, cnt in f_count.items():
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
#
#             if cnt > npix0:
#
#
#
#
# #
            cindex = ifil+1
#         for fil, cnt in f_count.items():

            if cnt > npix0 :
                print(fil, cnt)
                fof_result = tdet[tdet['fn3'] == fil]
                #value_index = tdet['fn'].index(fil)
                ctemp_index = np.where(tdet['fn3']==fil)
                coloridx[ctemp_index] = cindex


                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # for c in range(0, len(tdet)):
                #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                a += 1
                ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
            print('color index = {}, {}'.format(cindex, clr_list[cindex]))

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FIVE: num. of step 3 filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)

        plt.show()
        plt.savefig(self._fof_wfile+'_five_step3_fils_nb{}.pdf'.format(nbeam))
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)




############################
# 3D plot for FIVE components
############################
## Dr. Yu’s 3D plots code
    def five_3D_plot(self, basic_m0=None, tdet=None, nbeam=None, five_fn=None):
        tdet = self.get_tdet()

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #basic_m0 = fits.open('/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v-520/sn3vs2ss2/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits')
        #hdu = basic_m0[0]
        #w_line = WCS(hdu.header)

        #tdet = read('/Users/Hyunju/Dropbox/FUNS_OrionB/Individual_Filaments/2step_modify/Merged_FoF_result_2nd_mod.dat')

        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        #line, vtype, snrchoice, ssmo, vsmo, v1, v2 = input("Enter parameters: 'line' vtype snr ssmo vsmo v1 v2:").split()
        # ex) Enter parameters: 'line' snr vtype ssmo vsmo v1 v2:N2HP v10 3 2 2 3 13
        #v1 = float(v1)
        #v2 = float(v2)
        #snrchoice = float(snrchoice)
        v1, v2 = input("Enter velocity range: v1 v2:").split()

        #ssmo = int(ssmo)
        #vsmo = int(vsmo)

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)

        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        # if not os.file.exists(m0filename):
        #     save_fits(m0filename, cube.m, cube.header2d, overwrite=True)
        # else:
        #     pass

        #m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
        #w_line = basic_m0.wcs2d




        #onebeam = 5.171479 # (47 × 47 × π) ÷ (4 ln(2) × 22 × 22) for OriB
        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        if five_fn is not None:
            nfil = np.max(tdet[five_fn])
        else:
            five_fn = 'fn'
            nfil = np.max(tdet[five_fn])

        f_count = {}
        for lst in tdet[five_fn]:
            try: f_count[lst] += 1
            except: f_count[lst] = 1

        #print(f_count)
        f_count_sort = sorted(f_count.items())


        fig = plt.figure(figsize=(5/xyr, 5/xyr))
        ax = fig.add_subplot(111, projection='3d')  # Axe3D object
        ax.set_xlim3d(np.min(tdet['rp'])-20, np.max(tdet['rp'])+20)
        ax.set_ylim3d(np.min(tdet['dp'])-20, np.max(tdet['dp'])+20)
        ax.set_zlim3d(v1, v2)

        # line_image = basic_m0[0]
        # iy, ix = np.mgrid[0:line_image.shape[1], 0:line_image.shape[2]]

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            #if f_count_sort[fil][0] == selfil[isf]:
            #    cnt = f_count_sort[fil][1]
            #    cindex = fil+1
        #for fil, cnt in f_count.items():
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
            cindex = ifil+1

            if cnt > npix0:
                print(fil, cnt)
                fof_result = tdet[tdet[five_fn] == fil]

                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.
                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # central scatter
                ax.scatter(fof_result['rp'], fof_result['dp'], fof_result['vp'], marker='.', s=10, color = clr_list[cindex], alpha=0.3)


                xflat = np.full_like(fof_result['rp'], min(ax.get_xlim()))
                yflat = np.full_like(fof_result['dp'], max(ax.get_ylim()))
                zflat = np.full_like(fof_result['vp'], min(ax.get_zlim()))

                ax.set_xlabel('RA')
                ax.set_ylabel('Dec.')
                ax.set_zlabel('Velocity')

                a += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FIVE: num. of step 3 filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)
        plt.show()
        # fig.savefig('/Users/Hyunju/Dropbox/FUNS_OrionB/Plot/HFoF_3d_view.pdf', dpi=300, max_dpi=500)



############################
# 3D plot for FIVE components
############################
## Dr. Yu’s 3D plots code
    def five_3D_plot_selfil(self, basic_m0=None, tdet=None, nbeam=None, five_fn=None, selfil=None):
        tdet = self.get_tdet()

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #basic_m0 = fits.open('/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v-520/sn3vs2ss2/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits')
        #hdu = basic_m0[0]
        #w_line = WCS(hdu.header)

        #tdet = read('/Users/Hyunju/Dropbox/FUNS_OrionB/Individual_Filaments/2step_modify/Merged_FoF_result_2nd_mod.dat')

        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        #line, vtype, snrchoice, ssmo, vsmo, v1, v2 = input("Enter parameters: 'line' vtype snr ssmo vsmo v1 v2:").split()
        # ex) Enter parameters: 'line' snr vtype ssmo vsmo v1 v2:N2HP v10 3 2 2 3 13
        #v1 = float(v1)
        #v2 = float(v2)
        #snrchoice = float(snrchoice)
        v1, v2 = input("Enter velocity range: v1 v2:").split()

        #ssmo = int(ssmo)
        #vsmo = int(vsmo)

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize
        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)

        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        # if not os.file.exists(m0filename):
        #     save_fits(m0filename, cube.m, cube.header2d, overwrite=True)
        # else:
        #     pass

        #m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
        #w_line = basic_m0.wcs2d

        sf_ra = []
        sf_dec = []
        for isf in selfil:
            sfmask = tdet['fn3'] == isf
            sf_ra.append(np.min(np.unique(tdet['rp'][sfmask])))
            sf_ra.append(np.max(np.unique(tdet['rp'][sfmask])))
            sf_dec.append(np.min(np.unique(tdet['dp'][sfmask])))
            sf_dec.append(np.max(np.unique(tdet['dp'][sfmask])))


        #onebeam = 5.171479 # (47 × 47 × π) ÷ (4 ln(2) × 22 × 22) for OriB
        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        if five_fn is not None:
            nfil = np.max(tdet[five_fn])
        else:
            five_fn = 'fn'
            nfil = np.max(tdet[five_fn])

        f_count = {}
        for lst in tdet[five_fn]:
            try: f_count[lst] += 1
            except: f_count[lst] = 1

        #print(f_count_temp)
        f_count_sort = sorted(f_count.items())





        fig = plt.figure(figsize=(5/xyr, 5/xyr))
        ax = fig.add_subplot(111, projection='3d')  # Axe3D object
        ax.set_xlim3d(np.min(sf_ra)-20, np.max(sf_ra)+20)
        ax.set_ylim3d(np.min(sf_dec)-20, np.max(sf_dec)+20)
        ax.set_zlim3d(v1, v2)

        # line_image = basic_m0[0]
        # iy, ix = np.mgrid[0:line_image.shape[1], 0:line_image.shape[2]]

        cindex = 0

        for isf in range(np.size(selfil)):
            for fil in range(np.shape(f_count_sort)[0]):
                if f_count_sort[fil][0] == selfil[isf]:
                    cnt = f_count_sort[fil][1]
                    cindex = fil+1
            #
            #
            # # cnt = f_count.get(selfil)
            # if cnt > npix0:
            #
            #     cindex += 1


                fof_result = tdet[tdet[five_fn] == selfil[isf]]

                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.
                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # central scatter
                ax.scatter(fof_result['rp'], fof_result['dp'], fof_result['vp'], marker='.', s=10, color = clr_list[cindex], alpha=0.3)


                xflat = np.full_like(fof_result['rp'], min(ax.get_xlim()))
                yflat = np.full_like(fof_result['dp'], max(ax.get_ylim()))
                zflat = np.full_like(fof_result['vp'], min(ax.get_zlim()))

                ax.set_xlabel('RA')
                ax.set_ylabel('Dec.')
                ax.set_zlabel('Velocity')

            print(selfil[isf], cnt)

        #plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'filament number: {}'.format(selfil), fontsize=11)
        #plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)
        plt.show()
        # fig.savefig('/Users/Hyunju/Dropbox/FUNS_OrionB/Plot/HFoF_3d_view.pdf', dpi=300, max_dpi=500)




#######################################################################
## Dr. Yu’s code for mapping the filaments in 2D (RA & Dec) plane (she used > 5 beams )
 #######################################################################

    def fof_2D_plot_fnum(self, basic_m0=None, tdet=None, nbeam=None, note='', overplot_m0=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
    'olive', 'peru', 'gold', 'coral', 'navy',\
    'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
    'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
   'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
   'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
   'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
   'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
   'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
   'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        #v1, v2 = input("Enter velocity range: v1 v2:").split()
        v1='2'
        v2='14'

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        fig1 = plt.figure(figsize=(8, 8/xyr))
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['fn']!=0]

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        coloridx = np.zeros(len(tdet['fn']))
        nfil = np.max(tdet['fn'])
        print('nfil =', nfil)

        f_count = {}
        for lst in tdet['fn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1



        f_count_sort = sorted(f_count.items())

        #print(f_count)

        ax1.set_xlabel('R.A.', fontsize=15)
        ax1.set_ylabel('Dec.', fontsize=15)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=15)
        ax1.set_xlim(-10, cube.nr+10)
        ax1.set_ylim(-10, cube.nd+10)
        ax1.set_xticks(np.arange(-10, cube.nr+10, 10))
        ax1.set_yticks(np.arange(-10, cube.nd+10, 10))

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
            cindex = ifil+1
            if cnt > npix0 :
                print(fil, cnt)
                fof_result = tdet[tdet['fn'] == fil]
                #value_index = tdet['fn'].index(fil)
                ctemp_index = np.where(tdet['fn']==fil)
                coloridx[ctemp_index] = cindex
                cindex += 1

                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # for c in range(0, len(tdet)):
                #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                ax1.annotate('c[{}]f[{}]'.format(cindex, fil),(fof_result['rp'][0], fof_result['dp'][0]), fontsize=9, color=clr_list[cindex], weight='bold')
                # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])

                a += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FoF: num. of filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)

        plt.show()
        dirname = str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'mww'
        plt.savefig(self._fof_wfile+'_fnum_nb{}.pdf'.format(nbeam))
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)




#######################################################################
## Dr. Yu’s code for mapping the filaments in 2D (RA & Dec) plane (she used > 5 beams )
 #######################################################################

    def SeqFF_2D_plot_fnum(self, basic_m0=None, tdet=None, nbeam=None, note='', overplot_m0=True, fnum_note=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
    'olive', 'peru', 'gold', 'coral', 'navy',\
    'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
    'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
   'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
   'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
   'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
   'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
   'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
   'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        #v1, v2 = input("Enter velocity range: v1 v2:").split()
        v1='2'
        v2='14'
        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        #filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        import pickle
        pickled_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/4_121921/cube_pickle/SerB_C18O_v06_snr2.5.pickle'
        if os.path.exists(pickled_file):
            with open(file=pickled_file, mode='rb') as fr:
                cube = pickle.load(fr)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        #cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        fig1 = plt.figure(figsize=(8, 8/xyr))
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['finfn']!=0]

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        coloridx = np.zeros(len(tdet['finfn']))
        nfil = np.max(tdet['finfn'])
        print('nfil =', nfil)

        f_count = {}
        for lst in tdet['finfn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1



        f_count_sort = sorted(f_count.items())

        #print(f_count)

        ax1.set_xlabel('R.A.', fontsize=15)
        ax1.set_ylabel('Dec.', fontsize=15)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=15)
        ax1.set_xlim(-10, cube.nr+10)
        ax1.set_ylim(-10, cube.nd+10)
        ax1.set_xticks(np.arange(-10, cube.nr+10, 10))
        ax1.set_yticks(np.arange(-10, cube.nd+10, 10))

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
            #cindex = ifil+1
            cindex = fil
            if cnt > npix0 :
                print(fil, cnt)
                fof_result = tdet[tdet['finfn'] == fil]
                #value_index = tdet['fn'].index(fil)
                ctemp_index = np.where(tdet['finfn']==fil)
                coloridx[ctemp_index] = cindex
                #cindex += 1

                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # for c in range(0, len(tdet)):
                #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                if fnum_note == True:
                    ax1.annotate('c[{}]f[{}]'.format(cindex, fil),(fof_result['rp'][0], fof_result['dp'][0]), fontsize=9, color=clr_list[cindex], weight='bold')
                # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])

                a += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FoF: num. of filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)

        plt.show()
        dirname = str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'mww'
        if fnum_note == True:
            plt.savefig(self._fof_wfile+'_fnum_nb{}.pdf'.format(nbeam))
        else:
            plt.savefig(self._fof_wfile+'_nb{}.pdf'.format(nbeam))
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)




#######################################################################
## Dr. Yu’s code for mapping the filaments in 2D (RA & Dec) plane (she used > 5 beams )
 #######################################################################

    def SeqFF_2D_plot_fnum_selfil(self, basic_m0=None, tdet=None, nbeam=None, note='', overplot_m0=True, fnum_note=True, selfil=None):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
    'olive', 'peru', 'gold', 'coral', 'navy',\
    'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
    'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
   'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
   'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
   'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
   'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
   'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
   'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        #v1, v2 = input("Enter velocity range: v1 v2:").split()
        v1='2'
        v2='14'
        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        #filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        import pickle
        pickled_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/4_121921/cube_pickle/SerB_C18O_v06_snr2.5.pickle'
        if os.path.exists(pickled_file):
            with open(file=pickled_file, mode='rb') as fr:
                cube = pickle.load(fr)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        #cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        #print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)

        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            #ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            #ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            #ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            #ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['finfn']!=0]

        sf_ra = []
        sf_dec = []
        for isf in selfil:
            sfmask = tdet['finfn'] == isf
            sf_ra.append(np.min(np.unique(tdet['rp'][sfmask])))
            sf_ra.append(np.max(np.unique(tdet['rp'][sfmask])))
            sf_dec.append(np.min(np.unique(tdet['dp'][sfmask])))
            sf_dec.append(np.max(np.unique(tdet['dp'][sfmask])))


        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam


        #nfil = np.max(tdet['finfn'])

        # coloridx = np.zeros(len(tdet['finfn']))
        # nfil = np.max(tdet['finfn'])
        # print('nfil =', nfil)

        f_count = {}
        for lst in tdet['finfn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1



        f_count_sort = sorted(f_count.items())

        #print(f_count)
        xyr = np.size(sf_dec)/np.size(sf_ra)
        fig1 = plt.figure(figsize=(5, 5*xyr))
        ax1 = fig1.add_subplot(111)
        ax1.set_xlabel('R.A.', fontsize=12)
        ax1.set_ylabel('Dec.', fontsize=12)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=11)
        ax1.set_xlim(np.min(sf_ra)-8, np.max(sf_ra)+8)
        ax1.set_ylim(np.min(sf_dec)-8, np.max(sf_dec)+8)
        ax1.set_xticks(np.arange(np.min(sf_ra)-8, np.max(sf_ra)+8, 5))
        ax1.set_yticks(np.arange(np.min(sf_dec)-8, np.max(sf_dec)+8, 5))

        cindex = 0
        #a = 0
        for isf in selfil:
            for ifil in range(np.shape(f_count_sort)[0]):
                if f_count_sort[ifil][0] == isf:
                    fil = f_count_sort[ifil][0]
                    cnt = f_count_sort[ifil][1]
                    #cindex = ifil+1
                    cindex = fil
                    if cnt > npix0 :
                        print(fil, cnt)
                        fof_result = tdet[tdet['finfn'] == fil]
                        #value_index = tdet['fn'].index(fil)
                        #ctemp_index = np.where(tdet['finfn']==fil)
                        #coloridx[ctemp_index] = cindex
                        #cindex += 1

                        fil_array = np.zeros((cube.nd, cube.nr))
                        for i in range(len(fof_result)):
                            fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                        kernel = Gaussian2DKernel(0.4)
                        fil_array = convolve(fil_array, kernel)

                        # for c in range(0, len(tdet)):
                        #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                        #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                        # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                        # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                        ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.35, color = clr_list[cindex])
                        ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                        if fnum_note == True:
                            ax1.annotate('f[{}]'.format(fil),(fof_result['rp'][0], fof_result['dp'][0]-2), fontsize=9, color=clr_list[cindex], weight='bold')
                        # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])

                #a += 1

        #plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        #plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        #plt.gcf().text(0.25, 0.94, r'FoF: num. of filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        #plt.gcf().text(0.1, 0.96, r'nbeam:'+str(nbeam), fontsize=11)
        plt.gcf().text(0.1, 0.96, r'filament number: {}'.format(selfil), fontsize=11)

        plt.show()
        dirname = str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'mww'
        if fnum_note == True:
            plt.savefig(self._fof_wfile+'_fnum_nb{}_selfil_{}.pdf'.format(nbeam, selfil))
        else:
            plt.savefig(self._fof_wfile+'_nb{}_selfil.pdf'.format(nbeam))
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)




############################
# 3D plot for SeqFF components
############################
## Dr. Yu’s 3D plots code
    def SeqFF_3D_plot(self, basic_m0=None, tdet=None, nbeam=None, five_fn=None):
        tdet = self.get_tdet()

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #basic_m0 = fits.open('/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v-520/sn3vs2ss2/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits')
        #hdu = basic_m0[0]
        #w_line = WCS(hdu.header)

        #tdet = read('/Users/Hyunju/Dropbox/FUNS_OrionB/Individual_Filaments/2step_modify/Merged_FoF_result_2nd_mod.dat')

        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        #line, vtype, snrchoice, ssmo, vsmo, v1, v2 = input("Enter parameters: 'line' vtype snr ssmo vsmo v1 v2:").split()
        # ex) Enter parameters: 'line' snr vtype ssmo vsmo v1 v2:N2HP v10 3 2 2 3 13
        #v1 = float(v1)
        #v2 = float(v2)
        #snrchoice = float(snrchoice)
        v1, v2 = input("Enter velocity range: v1 v2:").split()

        #ssmo = int(ssmo)
        #vsmo = int(vsmo)
        import pickle
        pickled_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/4_121921/cube_pickle/SerB_C18O_v06_snr2.5.pickle'
        if os.path.exists(pickled_file):
            with open(file=pickled_file, mode='rb') as fr:
                cube = pickle.load(fr)



        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize
        #filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        #cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)

        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        # if not os.file.exists(m0filename):
        #     save_fits(m0filename, cube.m, cube.header2d, overwrite=True)
        # else:
        #     pass

        #m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
        #w_line = basic_m0.wcs2d




        #onebeam = 5.171479 # (47 × 47 × π) ÷ (4 ln(2) × 22 × 22) for OriB
        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        if five_fn is not None:
            nfil = np.max(tdet[five_fn])
        else:
            five_fn = 'fn'
            nfil = np.max(tdet[five_fn])

        f_count = {}
        for lst in tdet[five_fn]:
            try: f_count[lst] += 1
            except: f_count[lst] = 1

        #print(f_count)
        f_count_sort = sorted(f_count.items())


        fig = plt.figure(figsize=(5/xyr, 5/xyr))
        ax = fig.add_subplot(111, projection='3d')  # Axe3D object
        ax.set_xlim3d(np.min(tdet['rp'])-20, np.max(tdet['rp'])+20)
        ax.set_ylim3d(np.min(tdet['dp'])-20, np.max(tdet['dp'])+20)
        ax.set_zlim3d(v1, v2)

        # line_image = basic_m0[0]
        # iy, ix = np.mgrid[0:line_image.shape[1], 0:line_image.shape[2]]

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            #if f_count_sort[fil][0] == selfil[isf]:
            #    cnt = f_count_sort[fil][1]
            #    cindex = fil+1
        #for fil, cnt in f_count.items():
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
            cindex = fil
            if fil != 0:
                if cnt > npix0:
                    print(fil, cnt)
                    fof_result = tdet[tdet[five_fn] == fil]

                    fil_array = np.zeros((cube.nd, cube.nr))
                    for i in range(len(fof_result)):
                        fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.
                    kernel = Gaussian2DKernel(0.4)
                    fil_array = convolve(fil_array, kernel)

                    # central scatter
                    ax.scatter(fof_result['rp'], fof_result['dp'], fof_result['vp'], marker='.', s=10, color = clr_list[cindex], alpha=0.3)


                    xflat = np.full_like(fof_result['rp'], min(ax.get_xlim()))
                    yflat = np.full_like(fof_result['dp'], max(ax.get_ylim()))
                    zflat = np.full_like(fof_result['vp'], min(ax.get_zlim()))

                    ax.set_xlabel('RA')
                    ax.set_ylabel('Dec.')
                    ax.set_zlabel('Velocity')

                    a += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FIVE: num. of step 3 filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)
        plt.show()
        # fig.savefig('/Users/Hyunju/Dropbox/FUNS_OrionB/Plot/HFoF_3d_view.pdf', dpi=300, max_dpi=500)



############################
# 3D plot for SeqFF components
############################
## Dr. Yu’s 3D plots code
    def SeqFF_3D_plot_selfil(self, basic_m0=None, tdet=None, nbeam=None, five_fn=None, selfil=None, fnum_note=False):
        tdet = self.get_tdet()

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #basic_m0 = fits.open('/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v-520/sn3vs2ss2/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits')
        #hdu = basic_m0[0]
        #w_line = WCS(hdu.header)

        #tdet = read('/Users/Hyunju/Dropbox/FUNS_OrionB/Individual_Filaments/2step_modify/Merged_FoF_result_2nd_mod.dat')

        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        #line, vtype, snrchoice, ssmo, vsmo, v1, v2 = input("Enter parameters: 'line' vtype snr ssmo vsmo v1 v2:").split()
        # ex) Enter parameters: 'line' snr vtype ssmo vsmo v1 v2:N2HP v10 3 2 2 3 13
        #v1 = float(v1)
        #v2 = float(v2)
        #snrchoice = float(snrchoice)
        v1, v2 = input("Enter velocity range: v1 v2:").split()

        #ssmo = int(ssmo)
        #vsmo = int(vsmo)

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize
        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'

        import pickle
        pickled_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/4_121921/cube_pickle/SerB_C18O_v06_snr2.5.pickle'
        if os.path.exists(pickled_file):
            with open(file=pickled_file, mode='rb') as fr:
                cube = pickle.load(fr)

        #cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)

        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        # if not os.file.exists(m0filename):
        #     save_fits(m0filename, cube.m, cube.header2d, overwrite=True)
        # else:
        #     pass

        #m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
        #w_line = basic_m0.wcs2d

        sf_ra = []
        sf_dec = []
        for isf in selfil:
            sfmask = tdet['finfn'] == isf
            sf_ra.append(np.min(np.unique(tdet['rp'][sfmask])))
            sf_ra.append(np.max(np.unique(tdet['rp'][sfmask])))
            sf_dec.append(np.min(np.unique(tdet['dp'][sfmask])))
            sf_dec.append(np.max(np.unique(tdet['dp'][sfmask])))


        #onebeam = 5.171479 # (47 × 47 × π) ÷ (4 ln(2) × 22 × 22) for OriB
        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        if five_fn is not None:
            nfil = np.max(tdet[five_fn])
        else:
            five_fn = 'fn'
            nfil = np.max(tdet[five_fn])

        f_count = {}
        for lst in tdet[five_fn]:
            try: f_count[lst] += 1
            except: f_count[lst] = 1

        #print(f_count_temp)
        f_count_sort = sorted(f_count.items())





        fig = plt.figure(figsize=(5/xyr, 5/xyr))
        ax = fig.add_subplot(111, projection='3d')  # Axe3D object
        ax.set_xlim3d(np.min(sf_ra)-20, np.max(sf_ra)+20)
        ax.set_ylim3d(np.min(sf_dec)-20, np.max(sf_dec)+20)
        ax.set_zlim3d(v1, v2)

        # line_image = basic_m0[0]
        # iy, ix = np.mgrid[0:line_image.shape[1], 0:line_image.shape[2]]

        cindex = 0

        for isf in selfil:
            for fil in range(np.shape(f_count_sort)[0]):
                if f_count_sort[fil][0] == isf:
                    cnt = f_count_sort[fil][1]
                    #cindex = fil+1
                    cindex = f_count_sort[fil][0]
            #
            #
            # # cnt = f_count.get(selfil)
            # if cnt > npix0:
            #
            #     cindex += 1


                    fof_result = tdet[tdet[five_fn] == isf]

                    fil_array = np.zeros((cube.nd, cube.nr))
                    for i in range(len(fof_result)):
                        fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.
                    kernel = Gaussian2DKernel(0.4)
                    fil_array = convolve(fil_array, kernel)

                    # central scatter
                    ax.scatter(fof_result['rp'], fof_result['dp'], fof_result['vp'], marker='.', s=10, color = clr_list[cindex], alpha=0.45)
                    if fnum_note == True:
                        ax.text(fof_result['rp'][0]-5, fof_result['dp'][0]-5, fof_result['vp'][0]-0.5, "f[{}]".format(str(cindex)), fontsize=11, color=clr_list[cindex], weight='bold')



                    xflat = np.full_like(fof_result['rp'], min(ax.get_xlim()))
                    yflat = np.full_like(fof_result['dp'], max(ax.get_ylim()))
                    zflat = np.full_like(fof_result['vp'], min(ax.get_zlim()))

                    ax.set_xlabel('RA')
                    ax.set_ylabel('Dec.')
                    ax.set_zlabel('Velocity')

            print(isf, cnt)


        #plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'filament number: {}'.format(selfil), fontsize=11)
        #plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)
        plt.show()
        # fig.savefig('/Users/Hyunju/Dropbox/FUNS_OrionB/Plot/HFoF_3d_view.pdf', dpi=300, max_dpi=500)





##################################################################################
# Plot Filaments which are in each correponding velocity range on channel maps
#################################################################################
    def SeqFF_plot_vc_fils_on_chmap(self, basic_m0=None, tdet=None, nbeam=None, note=''):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        t_fil = self.load_filaments_info()
        n_fil = len(t_fil)


        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename
        from multiaxes import Multiaxes


        #v1, v2 = input("Enter velocity range: v1 v2:").split()

        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)


        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))


        # set map value_range, colormap
        int_vr = [0., 0.13]
        #cmap = plt.get_cmap('nipy_spectral')
        #cmap.set_bad('grey')
        cmap = plt.get_cmap('Greys')
        cmap.set_bad('white')

        maps, labels = cube.chmap(50, (5.5, 10.5))

        pdfname = self._fof_wfile+'_nb{}_chmap_vp.pdf'.format(nbeam)
        with PdfPages(pdfname) as pdf:
            for i in range(5):
                # draw figure using Multiaxes
                mx = Multiaxes(2, 5, 2, xyr, 0.36, 0.56, margin=(0.03, 0.03, 0.75, 0.03), scale=0.7, proj=cube.wcs2d)
                mx.shareaxes((True, True), 0.1)
                fig, ax, _ = mx.drawfig()
                cax = mx.sharecolorbar('right', 0.15, 0.1)
                for j in range(10):
                    iax = ax[j//5, j%5]
                    k = i*10+j
                    cs = iax.imshow(maps[k], cmap=cmap, vmin=int_vr[0], vmax=int_vr[1])
                    iax.annotate(labels[k]+' km/s', xy=(0.55, 0.05), xycoords='axes fraction', \
                                 ha='center', va='baseline', color='black', fontsize='x-small')
                    iax.tick_params(axis='x', labelsize=5)
                    iax.tick_params(axis='y', labelsize=5)

                    if j == 5:
                        iax.coords[0].set_axislabel('R.A. (J2000)', fontsize=7)
                        iax.coords[1].set_axislabel('Dec. (J2000)', fontsize=7)
                    else:
                        iax.coords[0].set_axislabel(' ')
                        iax.coords[1].set_axislabel(' ')

                    for kk in range(n_fil):
                        fof_result = t_fil[kk+1]
                        if len(fof_result) > npix0:
                            fvp = fof_result['vp']
                            vp_min = np.quantile(fvp, 0., axis=0)
                            vp_max = np.quantile(fvp, 1., axis=0)
                            ch_vr = [5.5+(k*0.1), 5.5+(k+1)*0.1]
                            if (vp_min <= ch_vr[0]) and (vp_max >= ch_vr[1]):
                                fil_array = np.zeros((cube.nd, cube.nr))
                                for ii in range(len(fof_result)):
                                    fil_array[fof_result['dp'][ii], fof_result['rp'][ii]] = 1.

                                kernel = Gaussian2DKernel(0.4)
                                fil_array = convolve(fil_array, kernel)
                                iax.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=0.5, alpha=0.1, color = clr_list[fof_result['finfn']])
                                iax.contour(fil_array, levels=[0.5], alpha=0.5, colors=clr_list[fof_result['finfn']])

                plt.colorbar(cs, cax)
                cax.tick_params(labelsize=5)
                cax.set_ylabel(r'Integrated intensity (K km s$^{-1}$)', fontsize=7)

                pdf.savefig()
                plt.close()




#######################################################################
## Dr. Yu’s code for mapping the filaments in 2D (RA & Dec) plane (she used > 5 beams )
 #######################################################################

    def SeqFF_2D_plot_cores_fnum(self, basic_m0=None, tdet=None, nbeam=None, note='', overplot_m0=True, fnum_note=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
              'olive', 'peru', 'gold', 'coral', 'navy',\
              'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
              'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
             'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
             'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
             'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
             'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
             'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
             'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
           'olive', 'peru', 'gold', 'coral', 'navy',\
           'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
           'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
          'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
          'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
          'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
          'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
          'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
          'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'olive', 'peru', 'gold', 'coral', 'navy',\
                  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                 'olive', 'peru', 'gold', 'coral', 'navy',\
                 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
                'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
                'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
                'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
                'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
                'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
                'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
                'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
               'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
               'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
             'olive', 'peru', 'gold', 'coral', 'navy',\
             'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
             'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
            'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
            'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
            'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
            'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
            'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
            'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
    'olive', 'peru', 'gold', 'coral', 'navy',\
    'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
    'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
   'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
   'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
   'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
   'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
   'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
   'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
       'olive', 'peru', 'gold', 'coral', 'navy',\
       'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
       'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
      'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
      'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
      'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
      'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
      'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
      'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
      'olive', 'peru', 'gold', 'coral', 'navy',\
      'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
      'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
     'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
     'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
     'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
     'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
     'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
     'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
   'olive', 'peru', 'gold', 'coral', 'navy',\
   'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
   'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
  'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
  'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
  'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
  'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
  'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
  'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'olive', 'peru', 'gold', 'coral', 'navy',\
          'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
          'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
         'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
         'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
         'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
         'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
         'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
         'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
         'olive', 'peru', 'gold', 'coral', 'navy',\
         'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
         'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
        'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
        'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
        'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
        'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
        'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
        'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
        'olive', 'peru', 'gold', 'coral', 'navy',\
        'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
        'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
       'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
       'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
       'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
       'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
       'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
       'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
     'olive', 'peru', 'gold', 'coral', 'navy',\
     'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
     'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
    'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
    'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
    'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
    'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
    'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
    'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
  'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
  'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
 'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
 'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
 'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
 'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
 'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
 'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
 'olive', 'peru', 'gold', 'coral', 'navy',\
 'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
 'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood', \
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood',
'olive', 'peru', 'gold', 'coral', 'navy',\
'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
'darkslateblue', 'mediumslateblue', 'mediumvioletred', 'chartreuse', 'burlywood']


        #to use Cube2map
        from SerB_Cube import SerBCube2map as Cube2map
        from run_funstools_Serp import get_crop_filename

        #v1, v2 = input("Enter velocity range: v1 v2:").split()
        v1='2'
        v2='14'
        if self._vtype == 'v10':
            rs = 200
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
        elif self._vtype == 'v06':
            #rs = 300
            rs = 300
            #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        rs = self._rmssize

        #filename = get_crop_filename(line=self._line, vtype=self._vtype, region=self._region)

        import pickle
        pickled_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/4_121921/cube_pickle/SerB_C18O_v06_snr2.5.pickle'
        if os.path.exists(pickled_file):
            with open(file=pickled_file, mode='rb') as fr:
                cube = pickle.load(fr)

        #### by using SerBCube2map.line_scan to easily check line profile in nxn map: n is odd number greater than 3.
        #filename = get_filename(line=self._line, vtype=self._vtype, ctype='matched')
        #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_match_cube.fits'
        #cube = Cube2map(filename, rmssize=rs, snr=self._snr, velocity_smo=self._vsm_find, spatial_smo=self._ssm_find)
        xyr = float(cube.nr/cube.nd)
        ssnr=str(int(self._snr))
        svs=str(int(self._vsm_find))
        ssm=str(int(self._ssm_find))

        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_temp/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
        print(m0filename)
        # m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region+self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/m0_C18O/'+self._region++self._line+'_'+str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'vr'+v1+v2+'_m0.fits'
#        m0filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir_v'+v1+v2+'/sn'+ssnr+'vs'+svs+'ss'+ssm+'/'+self._vtype+'/SerB_'+self._line+'_'+self._vtype+'_match_snr3_vsmo2_ssmo2_m0.fits'
        v1 = float(v1)
        v2 = float(v2)
        fig1 = plt.figure(figsize=(8, 8/xyr))
        if basic_m0 is not None:
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        elif os.path.exists(m0filename):
            basic_m0 = fits.open(m0filename)
            hdu = basic_m0[0]
            w_line = WCS(hdu.header)
            ax1 = fig1.add_subplot(111, projection = w_line)
        else:
            basic_m0 = cube.moment0(vr=[v1,v2])
            save_fits(m0filename, cube.m0, cube.header2d, overwrite=True)
            hdu = cube.header
            w_line = WCS(hdu)
            ax1 = fig1.add_subplot(111, projection = cube.wcs2d)

        if overplot_m0==True:
            imtemp = hdu.data
            imtemp[np.isnan(imtemp)]=0
            #tim = np.transpose(imtemp)
            ax1.imshow(imtemp, origin='lower', vmin=0, cmap=plt.get_cmap('binary'))




        #w_line = basic_m0.wcs2d

        tdet_all = self.get_tdet()
        tdet = tdet_all[tdet_all['finfn']!=0]

        onebeam = 6.324 # (47.26 × 47.26 × π) ÷ (4 ln(2) × 19.999122 × 19.999122) for SerB
        if nbeam is None:
            nbeam = 5
        else:
            nbeam = nbeam
        npix0 = nbeam * onebeam

        coloridx = np.zeros(len(tdet['finfn']))
        nfil = np.max(tdet['finfn'])
        print('nfil =', nfil)

        f_count = {}
        for lst in tdet['finfn']:
            try: f_count[lst] += 1
            except: f_count[lst] = 1



        f_count_sort = sorted(f_count.items())

        #print(f_count)

        ax1.set_xlabel('R.A.', fontsize=15)
        ax1.set_ylabel('Dec.', fontsize=15)

        # ax1.annotate(fil, xy=(fof_result[0][0], fof_result[0][1]), fontsize=10, color=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])

        #ax1.set_title(note+', nfil='+str(nfil), fontsize=15)
        ax1.set_title(note, fontsize=15)
        ax1.set_xlim(-10, cube.nr+10)
        ax1.set_ylim(-10, cube.nd+10)
        ax1.set_xticks(np.arange(-10, cube.nr+10, 10))
        ax1.set_yticks(np.arange(-10, cube.nd+10, 10))

        cindex = 0
        a = 0
        for ifil in range(np.shape(f_count_sort)[0]):
            fil = f_count_sort[ifil][0]
            cnt = f_count_sort[ifil][1]
            #cindex = ifil+1
            cindex = fil
            if cnt > npix0 :
                print(fil, cnt)
                fof_result = tdet[tdet['finfn'] == fil]
                #value_index = tdet['fn'].index(fil)
                ctemp_index = np.where(tdet['finfn']==fil)
                coloridx[ctemp_index] = cindex
                #cindex += 1

                fil_array = np.zeros((cube.nd, cube.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.

                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)

                # for c in range(0, len(tdet)):
                #     ax1.scatter(tdet['col1'][c], tdet['col2'][c], marker='s', s=10, alpha=0.2,
                #                 color=clr_list[np.mod(tdet['col6'][c], len(clr_list))])


                # ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[np.mod(fof_result['fn'][0],len(clr_list))])
                # ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[np.mod(fof_result['fn'][0], len(clr_list))])
                ax1.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                ax1.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])
                if fnum_note == True:
                    ax1.annotate('c[{}]f[{}]'.format(cindex, fil),(fof_result['rp'][0], fof_result['dp'][0]), fontsize=9, color=clr_list[cindex], weight='bold')
                # print(clr_list[np.mod(fof_result['fn'][0],len(clr_list))])

                a += 1

        plt.gcf().text(0.1, 0.96, self._plot_note, fontsize=11)
        plt.gcf().text(0.1, 0.94, r'vr=['+str(v1)+', '+str(v2)+']', fontsize=11)
        plt.gcf().text(0.25, 0.94, r'FoF: num. of filaments > nbeam = {} /{} '.format(a, cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)




        ### overplot Cores
        from astropy import constants as co
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        map_ws = SkyCoord('18:27:40 +0:05:15', unit=(u.hourangle, u.deg), frame='fk5')
        map_en = SkyCoord('18:31:30 +1:34:25', unit=(u.hourangle, u.deg), frame='fk5')

        YSOs_table = '/Users/khkim/Dropbox/DATA/FUNS/Serp/Spitzer_c2d_GBS_YSOs/apjs519129t2_mrt.txt'
        Cores_table = '/Users/khkim/Dropbox/DATA/FUNS/HGBS/HGBS_serpens_derived_core_catalog.txt'
        Coresdata = np.loadtxt(Cores_table, skiprows = 0+1+2+3+4+5+6+7+6, usecols=(2, 3, 17), dtype='str', encoding = "ISO-8859-1")
        Core_pos = list()
        Core_class = list()
        Core_marker_color = list()
        cc1 = 'red'
        cc2 = 'red'#'goldenrod'
        cc3 = 'violet'
        cc4 = 'blue'

        for i in range(len(Coresdata)):
            if (int(Coresdata[i][2]) == 2) or (int(Coresdata[i][2]) == 4): # to select robust prestellar cores and protostellar core
                aaa = Coresdata[i][0]+' '+Coresdata[i][1]
                pos_temp = SkyCoord(aaa, unit=(u.hourangle, u.deg), frame='fk5')
                if (float(pos_temp.ra/u.deg) >= float(map_ws.ra/u.deg)) and \
                    (float(pos_temp.ra/u.deg) <= float(map_en.ra/u.deg)) and \
                    (float(pos_temp.dec/u.deg) >= float(map_ws.dec/u.deg)) and \
                    (float(pos_temp.dec/u.deg) <= float(map_en.dec/u.deg)):

                    Core_pos.append(aaa)
                    cclass = int(Coresdata[i][2])
                    Core_class.append(cclass)
                    if cclass == 1:
                        Core_marker_color.append(cc1) #unbound starless core
                    elif cclass == 2:
                        Core_marker_color.append(cc2) # prestellar core (robust)
                    elif cclass == 3:
                        Core_marker_color.append(cc3) # candidate prestellar (non-robust)
                    elif cclass == 4:
                        Core_marker_color.append(cc4) # protostellar core

        cpos = SkyCoord(Core_pos, unit=(u.hourangle, u.deg), frame='fk5')

        for j in range(len(cpos)):
            ax1.scatter(float(cpos[j].ra/u.deg), float(cpos[j].dec/u.deg), marker='o', edgecolor='k', \
            facecolor=Core_marker_color[j], transform=ax1.get_transform('fk5'), alpha=0.7)



        plt.show()
        dirname = str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'mww'
        if fnum_note == True:
            plt.savefig(self._fof_wfile+'_cores_fnum_nb{}.pdf'.format(nbeam))
        else:
            plt.savefig(self._fof_wfile+'_cores_nb{}.pdf'.format(nbeam))
        #u/Dropbox/FUNS_OrionB/Plot/HFoF_final.pdf', dpi=300, max_dpi=500)
