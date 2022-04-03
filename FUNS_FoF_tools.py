
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

    def __init__(self, decompose_key='', decompose_path='', decompose_result='', region='', line='C18O', vtype='v06', snr=None, ssm_find=None, vsm_find=None, ssm_fit=None, vsm_fit=None, sigma=None, v0=None, rd0=None, fof_wfile=''):
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
               'darkslateblue', 'mediumslateble', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateble', 'mediumvioletred', 'chartreuse', 'burlywood']


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


        fig = plt.figure(figsize=(7/xyr, 7/xyr))
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
               'darkslateblue', 'mediumslateble', 'mediumvioletred', 'chartreuse', 'burlywood', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'cornflowerblue', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lime', 'slateblue', 'olivedrab', 'firebrick', 'teal',  \
              'violet', 'darkviolet', 'indigo', 'rosybrown', 'cadetblue', \
              'darkslateblue', 'mediumslateble', 'mediumvioletred', 'chartreuse', 'burlywood']


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

        tdet = self.get_tdet()

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
        for fil, cnt in f_count.items():

            # if cnt > 15 and cnt < 100 :

            if cnt > npix0 :
            # if 5 < cnt < npix0 :

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
        plt.gcf().text(0.25, 0.94, r'num. of filaments = '+str(cindex), fontsize=11)
        plt.gcf().text(0.1, 0.92, r'nbeam:'+str(nbeam), fontsize=11)

        plt.show()
        dirname = str(self._ssm_find)+str(self._vsm_find)+str(self._ssm_fit)+str(self._vsm_fit)+self._vtype+'snr'+str(self._snr)+'mww'
        plt.savefig(self._fof_wfile+'.pdf')
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
