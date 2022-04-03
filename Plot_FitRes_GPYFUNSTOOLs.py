
"""
To plot spectrum and the decomposed components and/or Filament assingend components \
option: spatial smoothing, velocity smoothing parameter
        decomposed spectral componenets
        if FoF applied, FoF results
"""


"""
read a cube
check smoothing factor. smoothed fits are not exist, smooth and save to a new fits cube. (for later use, i.e. input cube to GPY+)

"""

import funstools
from funstools import get_rms, get_mask, get_det, wcs2d, header2d, save_fits
from funstools import smooth3d, radsmo2d, boxsmo1d, make_velo
from funstools import Cube2map
import numpy as np
from scipy.signal import find_peaks
from astropy.modeling.models import Gaussian1D
from astropy.table import Table, Column
from astropy.io.ascii import read, write
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
#import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, Circle
import matplotlib.colors as mcolors
import os
from astropy.convolution import convolve, Gaussian2DKernel

class PlotFitRes(Cube2map):

    def  __init__(self, cube=None, ext=None, region='SerB_A_', getrms='both', rmssize=300, max_rms=None,
                 snr=None, velocity_smo=None, spatial_smo=None, line = r'C$^{18}$O', key='', DC=False, FoF=False):
        Cube2map.__init__(self, cube, ext, getrms, rmssize, max_rms, snr, velocity_smo, spatial_smo, True, True)
        if velocity_smo is None:
            velocity_smo=3
        if snr is None:
            snr=3.0
        if spatial_smo is None:
            spatial_smo=3


        self._velocity_smo = float(velocity_smo)
        self._spatial_smo = float(spatial_smo)
        self._region = region
        self._line = line
        self._key = key
        #self._header = self.cube.header
        #self._data = self.cube.data
        _smdata = None
        #_px = None


    def save_smoothed_3dcube(self):
        if self._smdata is None:
            self._smdata = smooth3d(self.data, self._spatial_smo, self._velocity_smo)
        return self._smdata

    # def _px(self, d, r):
    #     return np.where(self.cube.comps[:, d, r] == 1)[0]


    def plot_dc_fit(self, fit, vr=None, yr=None, n_ch=None, nzp=None, vtype=None):
        """
        Interactive scan for fitting results.
        Plot fitting results over line profile on click position.
        Parameters:
            fit : astropy.table.Table
                Result table of decomposing.
                    Decompose.initial_fit_result
                    Decompose.second_fit_result
                    Decompose.final_fit_result
                    Decompose.decompose_result
            vr : list or tuple, optional
                Range of velocity axis for line plot.
            yr : list or tuple, optional
                Range of temperature axis for line plot.

            n_ch : number of pixels around the picked pixel including the center pixel.
               i.e., number of pixels in x axis and y axis.
               default: 3, available up to 9. (3, 5, 7, 9)

            nzp : number of pixels for the zoomed up img

        Returns:
            tuple (plt.figure, plt.subplot, plt.subplot)
        """

        plt.ion()
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        x = self.x
        if not isinstance(fit, Table):
            raise TypeError("'fit' is not astropy Table.".format(fit))
        if vr is None:
            vr = [x[self._rmssize], x[-self._rmssize]]
        elif isinstance(vr, list) or isinstance(vr, tuple):
            pass
        else:
            raise TypeError("'vr' is not list or tuple.".format(vr))
        if yr is None:
            yr = [-4*np.nanmedian(self.rms), np.nanmax(self.data)+np.nanmedian(self.rms)]
        elif isinstance(yr, list) or isinstance(yr, tuple):
            pass
        else:
            raise TypeError("'yr' is not list or tuple".format(yr))

        if nzp is None:
            nzp = 25
        elif isinstance(nzp, int):
            pass

        if vtype is None:
            vtype = 'v06'
        else:
            vtype = vtype



        # make map
        det = np.full(self.data.shape, 0.)
        det[np.nan_to_num(self.data) > 3*np.nanmedian(self.rms)] = 1.
        dsum = np.nansum(det, axis=(1, 2))
        imap = np.sum(self.data[dsum > 2*np.mean(dsum)], axis=0)
        # map drawing
        lsfig = plt.figure(figsize=(16, 9))
        lsmap = lsfig.add_axes([0.06, 0.4, 0.3, 0.58], projection=self.wcs2d)
        lsmap.set_xlabel('R.A. (J2000)')
        lsmap.set_ylabel('Dec. (J2000)')
        lsmap.imshow(imap, origin='lower', vmin=0, cmap='inferno')
        #self._pp = lsmap.scatter(-100, -100, s=60, color='lime')
        lsmap.set_xlim(-0.5, self.nr-0.5)
        lsmap.set_ylim(-0.5, self.nd-0.5)

        plt.gcf().text(0.01, 0.95, 'ssmo:'+str(int(self._spatial_smo))+' vsmo:'+str(int(self._velocity_smo)))
        plt.gcf().text(0.01, 0.92, 'snr:'+str(self._snr))
        plt.gcf().text(0.01, 0.89, r'line='+self._line)
        plt.gcf().text(0.01, 0.86, 'dc: '+self._key, color='green')
        plt.gcf().text(0.01, 0.83, 'srms*snr', color='red')
        plt.gcf().text(0.01, 0.80, 'rms*snr', color='gray')
        plt.gcf().text(0.01, 0.77, r'$\Sigma$(comps)', color='magenta')


        if n_ch is None:
            n_ch = 3
        if isinstance(n_ch, int):
            # line plot
            #lgs = lsfig.add_gridspec(nrows=n_ch, ncols=n_ch, left=0.37, right=0.98, wspace=0, hspace=0)
            lgs = lsfig.subplots(ncols=n_ch, nrows=n_ch, gridspec_kw={'left':0.37, 'right':0.98, 'bottom':0.07, 'top':0.98,'wspace':0, 'hspace':0})
            for a_ra in range(n_ch):
                for a_dec in range(n_ch):
                    lsfig.add_subplot(lgs[a_dec,a_ra])

            lin_index = int((n_ch/2)-0.5)

            def _onclick(event):
                if event.inaxes == lsmap.axes:
                    self._pr = int(round(event.xdata))
                    self._pd = int(round(event.ydata))
                elif event.inaxes:
                    for a_ra in range(n_ch):
                        for a_dec in range(n_ch):
                            add_ra = a_ra - lin_index
                            add_dec = a_dec - lin_index
                            if event.inaxes == lgs[a_dec, a_ra].axes:
                                self._pr += add_ra
                                self._pd -= add_dec
                else:
                    return
                r = self._pr
                d = self._pd

                # define the zoom up region
                zp = int(np.rint(nzp/2))
                dmin = d-zp
                dmax = d+zp
                rmin = r-zp
                rmax = r+zp

                num_pix_r = self.cube.header['NAXIS1']
                num_pix_d = self.cube.header['NAXIS2']
                if dmin < 0 : dmin = 0
                if dmax > num_pix_d-1 : dmax = num_pix_d-1
                if rmin < 0 : rmin = 0
                if rmax > num_pix_r-1 : rmax = num_pix_r-1
                # create a rectangle patch
                lsmap.patches=[]
                rect = Rectangle((rmin,dmin), nzp, nzp, linewidth=1,edgecolor='g',facecolor='none')
                lsmap.add_patch(rect)
                pcircle = Circle((r+0.5,d+0.5),radius=1, color='g', fill='True', alpha=0.7)
                lsmap.add_patch(pcircle)

                # local map zoomed up around the clicked pixel
                lsmap_zoom = lsfig.subplots(ncols=1, nrows=1, gridspec_kw={'left':0.06, 'right':0.36, 'bottom':0.07, 'top':0.34})

                lsmap_zoom.patches=[]
                imap_zoom = imap[dmin:dmax, rmin:rmax]
                imap_zoom_extent = [rmin, rmax, dmin, dmax]
                lsmap_zoom.set_xlabel('R.A.(pixel num.)')
                lsmap_zoom.set_ylabel('Dec.(pixel num.)')
                lsmap_zoom.imshow(imap_zoom, origin='lower', vmin=0, cmap='inferno', extent=imap_zoom_extent)
                zcircle = Circle((r+0.5,d+0.5),radius=0.5, color='g', fill='True', alpha=0.7)
                lsmap_zoom.add_patch(zcircle)

                for a_ra in range(n_ch):
                    for a_dec in range(n_ch):
                        ax = lgs[a_dec, a_ra]
                        ax.clear()
                        ax.set_xlim(*vr)
                        ax.set_ylim(*yr)
                        ax.tick_params(labelsize=7)
                        if a_ra != 0:
                            ax.set_yticklabels([])
                        if a_dec != n_ch-1:
                            ax.set_xticklabels([])
                        ri = a_ra - lin_index
                        di = a_dec - lin_index
                        rr = r+ri
                        dd = d-di
                        if 0 <= rr < self.nr and 0 <= dd < self.nd:
                            srms = self._snr*self.srms[dd, rr]
                            rms = self._snr*self.rms[dd, rr]
                            ax.plot([x[0], x[-1]], [0, 0], lw=1, alpha=0.5, color='black')
                            ax.plot([x[0], x[-1]], [srms, srms], ls='dotted', lw=1, alpha=0.3, color='red')
                            ax.plot([x[0], x[-1]], [rms, rms], ls='dotted', lw=1, alpha=0.5, color='gray')
                            #ax.step(x, self.mdata_for_finding[:, dd, rr], color='olive', lw=1, alpha=0.6)
                            ax.step(x, self.data[:, dd, rr], color='gray', lw=1)
                            ax.step(x, self.sdata[:, dd, rr], color='black', lw=1, alpha=0.7)
                            ct = fit[(fit['rp'] == rr) & (fit['dp'] == dd)]
                            if len(ct) == 0:
                                continue

                            cm = Gaussian1D(ct['tp'][0], ct['vp'][0], ct['sd'][0])

                            for c in range(1, len(ct)):
                                cm += Gaussian1D(ct['tp'][c], ct['vp'][c], ct['sd'][c])
                            if len(ct) == 1:
                                ax.plot(x, cm(x), color='green', lw=1, alpha=0.5)
                                ax.plot((ct['vp'][0], ct['vp'][0]), [-0.25,0], color='green', lw=1, alpha=0.5)
                            else:
                                for c in range(len(ct)):
                                    ax.plot(x, cm[c](x), color='green', lw=1, alpha=0.5)
                                    ax.plot((ct['vp'][c], ct['vp'][c]), [-0.25,0], color='green', lw=1, alpha=0.5)
                            ax.plot(x, cm(x), color='magenta', lw=2., alpha=0.5)


                lgs[n_ch-1, lin_index].set_xlabel(r'$V_\mathrm{LSR}$ [km/s]')
                lgs[lin_index,0].set_ylabel(r'$T_\mathrm{A}^* [K]$')
                c = self.wcs2d.pixel_to_world(r, d)
                c = c.to_string('hmsdms').split(' ')
                lgs[lin_index, lin_index].annotate('R.A.: {} ({})'.format(r, c[0]), (0.03, 0.92), xycoords='axes fraction', fontsize=6)
                lgs[lin_index, lin_index].annotate('Dec.: {} ({})'.format(d, c[1]), (0.03, 0.84), xycoords='axes fraction', fontsize=6)

                lgs[lin_index, lin_index].figure.canvas.draw()

                lsmap_zoom.remove()

            lsfig.canvas.mpl_connect('button_press_event', _onclick)
            return lsfig, lsmap, lgs[lin_index, lin_index]




    def plot_fof_profile(self, tdet=None, line='C18O', vr=None, yr=None, n_ch=None, nzp=None, vtype=None, nbeam=None, dckey=None, \
                          filsigmacut = None, dv0 = None, drd0 = None, fof_result = None, dcresult=None, dctype=' '):
        """
        Interactive scan for fitting results.
        Plot fitting results over line profile on click position.
        Parameters:
            line: the default line for decomposing and fof is C18O

            vr : list or tuple, optional
                Range of velocity axis for line plot.
            yr : list or tuple, optional
                Range of temperature axis for line plot.

            n_ch : number of pixels around the picked pixel including the center pixel.
               i.e., number of pixels in x axis and y axis.
               default: 3, available up to 9. (3, 5, 7, 9)

            nzp : number of pixels for the zoomed up img

        tdet : astropy.table.Table
               Result table of FoF
               tdet

        nbeam : number of beam

        dckey: the directory name which contains the table of the decomposed result

        filsigmacut: sigma level for FoF

        dv0: the minimal velocity difference between two filaments.
        if the velocity at the peak of gaussian component is less than dv0,
        it would be considerd as a same filmanet

        drd0: the maximum pixel distance between two filament components....

        fof_result: the output table of FoF result. It read the table and use it to draw the gaussian components identified as some filaments

        Returns:
            tuple (plt.figure, plt.subplot, plt.subplot
        """
        plt.ion()
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

    #--FF
        if dckey is None:
            dckey = 'sn3.0d3s1'
        else:
            pass
        if filsigmacut is None:
            filsigmacut = 4.0
        else:
            pass
        if dv0 is None:
            dv0 = 0.1
        else:
            pass
        if drd0 is None:
            drd0 = 2.0
        else:
            pass
        if vtype is None:
            vtype = 'v06'
        else:
            vtype = vtype

        cfilsigmacut = str(filsigmacut)
        cdv0 = str(dv0)
        cdrd0 = str(drd0)

        # read decomposed results
        dcresult_filename = dcresult
        fit = read(dcresult_filename)
        #FoF_save='/Users/khkim/Dropbox/Work5/FUNS-Serp/FOF/FoF_dc4/'
        if fof_result is not None:
            fof_result = fof_result
        else:
            raise TypeError('No FoF result table exist: Please run FoF first.')

        # read FoF results
        if tdet is not None:
            tdet = tdet
        elif os.path.exists(fof_result):
            tdet = read(fof_result)
        else:
            raise TypeError('No FoF result table exist: Please run FoF first.')


        # make map
        det = np.full(self.data.shape, 0.)
        det[np.nan_to_num(self.data) > filsigmacut*np.nanmedian(self.rms)] = 1.
        dsum = np.nansum(det, axis=(1, 2))
        imap = np.sum(self.data[dsum > 2*np.mean(dsum)], axis=0)
        # map drawing
        lsfig = plt.figure(figsize=(16, 9))
        lsmap = lsfig.add_axes([0.06, 0.4, 0.3, 0.58], projection=self.wcs2d)
        lsmap.set_xlabel('R.A. (J2000)')
        lsmap.set_ylabel('Dec. (J2000)')
        lsmap.imshow(imap, origin='lower', vmin=0, cmap='gist_yarg')
        lsmap.set_xlim(-0.5, self.nr-0.5)
        lsmap.set_ylim(-0.5, self.nd-0.5)


        plt.gcf().text(0.01, 0.95, 'ssmo:'+str(int(self._spatial_smo))+' vsmo:'+str(int(self._velocity_smo)))
        plt.gcf().text(0.01, 0.92, 'snr:'+str(self._snr))
        plt.gcf().text(0.01, 0.89, 'line='+self._line)
        plt.gcf().text(0.01, 0.86, 'dc: '+dctype, color='red')
        plt.gcf().text(0.01, 0.83, self._key, color='red')
        plt.gcf().text(0.01, 0.80, 'rms*snr', color='gray')
        plt.gcf().text(0.01, 0.77, 'srms*snr', color='darkblue')
        #plt.gcf().text(0.01, 0.77, r'$\Sigma$(comps)', color='magenta')

        plt.gcf().text(0.01, 0.47, 'FoF parameters:')
        plt.gcf().text(0.01, 0.44, 'sigcut:'+cfilsigmacut)
        plt.gcf().text(0.01, 0.41, 'dist_v0:'+cdv0)
        plt.gcf().text(0.01, 0.38, 'dist_rd0:'+cdrd0)



        prop_cycle = plt.rcParams['axes.prop_cycle']
        #colors = prop_cycle.by_key()['color']
        clr_list = ['dimgray', 'olive', 'peru', 'gold', 'coral', 'navy',\
                'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
                'darkkhaki', 'turquoise', 'pink', 'blue', 'silver', 'yellow', 'magenta', \
               'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
               'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
               'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
               'lightyellow', 'honeydew', 'lime', 'lightcyan', 'slateblue', \
               'violet', 'darkviolet', \
               'olive', 'peru', 'gold', 'coral', 'navy',\
               'blueviolet', 'crimson', 'darkorange', 'lawngreen', 'royalblue',\
               'darkkhaki', 'turquoise', 'pink', 'blue', 'silver', 'yellow', 'magenta', \
              'sienna', 'darkgoldenrod', 'lightseagreen', 'cyan', 'deeppink', \
              'salmon', 'khaki', 'darkseagreen', 'seagreen', 'teal', \
              'skyblue', 'dodgerblue', 'mediumpurple', 'brown', 'tomato', \
              'lightyellow', 'honeydew', 'lime', 'lightcyan', 'slateblue', \
              'violet', 'darkviolet']


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

        rp = tdet['rp']
        dp = tdet['dp']
        tp = tdet['tp']
        vp = tdet['vp']
        sd = tdet['sd']

        coloridx = np.zeros(len(tdet['fn']))
        cindex = 0
        for fil, cnt in f_count.items():
            if cnt > npix0:
                # # if 5 < cnt < npix0 :
                fof_result = tdet[(tdet['fn'] == fil)]
                # color = clr_list[np.mod(fof_result['fn'][0], len(clr_list))]
                # #color = clr_list[icolor]
                # print(fil, cnt, color)
                # cdic[fil] = color
                ctemp_index = np.where(tdet['fn']==fil)
                coloridx[ctemp_index] = cindex+1
                cindex += 1

                fil_array = np.zeros((self.nd, self.nr))
                for i in range(len(fof_result)):
                    fil_array[fof_result['dp'][i], fof_result['rp'][i]] = 1.
                kernel = Gaussian2DKernel(0.4)
                fil_array = convolve(fil_array, kernel)
                lsmap.scatter(fof_result['rp'], fof_result['dp'], marker='s', s=10, alpha=0.2, color = clr_list[cindex])
                lsmap.contour(fil_array, levels=[0.5], alpha=0.6, colors=clr_list[cindex])


        filcolor = Column(name='filcolor', data=coloridx, dtype='i4')
        tdet.add_column(filcolor)

    #--FF

        x = self.x
        if not isinstance(fit, Table):
            raise TypeError("'fit' is not astropy Table.".format(fit))
        if vr is None:
            vr = [x[self._rmssize], x[-self._rmssize]]
        elif isinstance(vr, list) or isinstance(vr, tuple):
            pass
        else:
            raise TypeError("'vr' is not list or tuple.".format(vr))
        if yr is None:
            yr = [-4*np.nanmedian(self.rms), np.nanmax(self.data)+np.nanmedian(self.rms)]
        elif isinstance(yr, list) or isinstance(yr, tuple):
            pass
        else:
            raise TypeError("'yr' is not list or tuple".format(yr))

        if nzp is None:
            nzp = 25
        elif isinstance(nzp, int):
            pass

        if n_ch is None:
            n_ch = 3
        if isinstance(n_ch, int):
            lgs = lsfig.subplots(ncols=n_ch, nrows=n_ch, gridspec_kw={'left':0.37, 'right':0.98, 'bottom':0.07, 'top':0.98,'wspace':0, 'hspace':0})
            for a_ra in range(n_ch):
                for a_dec in range(n_ch):
                    lsfig.add_subplot(lgs[a_dec,a_ra])

            lin_index = int((n_ch/2)-0.5)

            def _onclick(event):
                if event.inaxes == lsmap.axes:
                    self._pr = int(round(event.xdata))
                    self._pd = int(round(event.ydata))
                elif event.inaxes:
                    for a_ra in range(n_ch):
                        for a_dec in range(n_ch):
                            add_ra = a_ra - lin_index
                            add_dec = a_dec - lin_index
                            if event.inaxes == lgs[a_dec, a_ra].axes:
                                self._pr += add_ra
                                self._pd -= add_dec
                else:
                    return
                r = self._pr
                d = self._pd

                # define the zoom up region
                zp = int(np.rint(nzp/2))
                dmin = d-zp
                dmax = d+zp
                rmin = r-zp
                rmax = r+zp

                num_pix_r = self.cube.header['NAXIS1']
                num_pix_d = self.cube.header['NAXIS2']
                if dmin < 0 : dmin = 0
                if dmax > num_pix_d-1 : dmax = num_pix_d-1
                if rmin < 0 : rmin = 0
                if rmax > num_pix_r-1 : rmax = num_pix_r-1
                #
                # create a rectangle patch
                lsmap.patches=[]
                rect = Rectangle((rmin,dmin), nzp, nzp, linewidth=1,edgecolor='g',facecolor='none')
                lsmap.add_patch(rect)
                pcircle = Circle((r+0.5,d+0.5),radius=1, color='g', fill='True', alpha=0.7)
                lsmap.add_patch(pcircle)

                # local map zoomed up around the clicked pixel
                lsmap_zoom = lsfig.subplots(ncols=1, nrows=1, gridspec_kw={'left':0.06, 'right':0.36, 'bottom':0.07, 'top':0.34})

                lsmap_zoom.patches=[]
                imap_zoom = imap[dmin:dmax, rmin:rmax]
                imap_zoom_extent = [rmin, rmax, dmin, dmax]
                lsmap_zoom.set_xlabel('R.A.(pixel num.)')
                lsmap_zoom.set_ylabel('Dec.(pixel num.)')
                lsmap_zoom.imshow(imap_zoom, origin='lower', vmin=0, cmap='gist_yarg', extent=imap_zoom_extent)
                zcircle = Circle((r+0.5,d+0.5),radius=0.5, color='g', fill='True', alpha=0.7)
                lsmap_zoom.add_patch(zcircle)


                for a_ra in range(n_ch):
                    for a_dec in range(n_ch):
                        ax = lgs[a_dec, a_ra]
                        ax.clear()
                        ax.set_xlim(*vr)
                        ax.set_ylim(*yr)
                        ax.tick_params(labelsize=7)
                        if a_ra != 0:
                            ax.set_yticklabels([])
                        if a_dec != n_ch-1:
                            ax.set_xticklabels([])
                        ri = a_ra - lin_index
                        di = a_dec - lin_index
                        rr = r+ri
                        dd = d-di
                        if 0 <= rr < self.nr and 0 <= dd < self.nd:
                            srms = self._snr*self.srms[dd, rr]
                            rms = self._snr*self.rms[dd, rr]

                            ax.plot([x[0], x[-1]], [srms, srms], ls='dotted', lw=1, alpha=0.3, color='darkblue')
                            ax.plot([x[0], x[-1]], [rms, rms], ls='dotted', lw=1, alpha=0.5, color='gray')
                            ax.step(x, self.data[:, dd, rr], color='gray', lw=1)
                            ax.step(x, self.sdata[:, dd, rr], color='black', lw=1, alpha=0.7)
                            x_text_location = vr[0]+((vr[1]-vr[0])/12)
                            y_text_location = yr[1]-((yr[1]-yr[0])/10)
                            ax.text(x_text_location, y_text_location, "RA {}  Dec {}".format(str(rr), str(dd)), fontsize=6)
                            ct = tdet[(tdet['rp'] == rr) & (tdet['dp'] == dd)]

                            if len(ct) is not 0:
                                tct = ct
                                tcm = Gaussian1D(tct['tp'][0], tct['vp'][0], tct['sd'][0])
                            else:
                                tct = fit[(fit['rp'] == rr) & (fit['dp'] == dd)]

                            #tcm = Gaussian1D(tct['tp'][0], tct['vp'][0], tct['sd'][0])

                            if len(ct) == 0 and len(tct) != 0:
                                for c in range(len(tct)):
                                    rc = tct[c][0]
                                    dc = tct[c][1]
                                    tc = tct[c][4]
                                    vc = tct[c][5]
                                    sc = tct[c][6]
                                    #fc = tct[c][10]
                                    colc = 0

                                    cm = Gaussian1D(tc, vc, sc)

                                    ax.plot(x, cm(x), color=clr_list[colc], lw=2, alpha=0.5)
                                    filinfo2txt=str(int(colc))+' [ -- ]:'+str(round(vc, 3))
                                    ax.text(x_text_location, (y_text_location-0.5)-0.3*c, filinfo2txt, fontsize=6, color=clr_list[colc])


                            if len(ct) is not 0:
                                for c in range(1,len(ct)):
                                    tcm += Gaussian1D(tct['tp'][c], tct['vp'][c], tct['sd'][c])

                            if len(ct) == 1:
                                colc = ct[0][12]
                                fc = ct[0][10]
                                vc = ct[0][5]
                                ax.plot(x, tcm(x), color=clr_list[colc], lw=2, alpha=0.7)
                                filinfo2txt=str(int(colc))+' ['+str(fc)+']:'+str(round(vc, 3))
                                ax.text(x_text_location, (y_text_location-0.5), filinfo2txt, fontsize=6, color=clr_list[colc])


                            if len(ct) > 1:
                                ax.plot(x, tcm(x), color='red', lw=1, alpha=0.4)
                                for c in range(len(ct)):
                                    rc = ct[c][0]
                                    dc = ct[c][1]
                                    tc = ct[c][4]
                                    vc = ct[c][5]
                                    sc = ct[c][6]
                                    fc = ct[c][10]
                                    colc = ct[c][12]

                                    cm = Gaussian1D(tc, vc, sc)

                                    ax.plot(x, cm(x), color=clr_list[colc], lw=2, alpha=0.7)
                                    filinfo2txt=str(int(colc))+' ['+str(fc)+']:'+str(round(vc, 3))
                                    ax.text(x_text_location, (y_text_location-0.5)-0.3*c, filinfo2txt, fontsize=6, color=clr_list[colc])

                        else:
                            ax.plot([x[0],x[-1]], [0,0], color='k')

                lgs[n_ch-1, lin_index].set_xlabel(r'$V_\mathrm{LSR}$ [km/s]')
                lgs[lin_index,0].set_ylabel(r'$T_\mathrm{A}^* [K]$')
                c = self.wcs2d.pixel_to_world(r, d)
                c = c.to_string('hmsdms').split(' ')
                #lgs[lin_index, lin_index].annotate('R.A.: {} ({})'.format(r, c[0]), (0.03, 0.92), xycoords='axes fraction', fontsize=6)
                #lgs[lin_index, lin_index].annotate('Dec.: {} ({})'.format(d, c[1]), (0.03, 0.84), xycoords='axes fraction', fontsize=6)
                #lgs[lin_index, lin_index].annotate('R.A.: {}'.format(c[0]), (0.03, 0.92), xycoords='axes fraction', fontsize=6)
                #lgs[lin_index, lin_index].annotate('Dec.: {}'.format(c[1]), (0.03, 0.84), xycoords='axes fraction', fontsize=6)

                lgs[lin_index, lin_index].figure.canvas.draw()

                lsmap_zoom.remove()

        lsfig.canvas.mpl_connect('button_press_event', _onclick)
        return lsfig, lsmap, lgs[lin_index, lin_index]
