import funstools
from funstools import get_rms
from funstools import Cube2map
from funstools import save_fits
from matplotlib import pyplot as plt
from funstools import full_line_scan
from funstools import Decompose


#from .cube2map import Cube2map
import os
from glob import glob
from matplotlib import pyplot as plt
import numpy as np



#import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import csv
import re





###############################################################################
#data maps: rms map, m0 map, and m1 map
###############################################################################
matched_cubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'
org_cubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/org_cube/'

#output directory
default_outdir1='/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/outdir/default/'
v10_dir = 'v10/'
v06_dir = 'v06/'

vtype=''
side=''
snrchoice=''
vsmo=''
ssmo=''
ctype=''
csize=''
outdir = default_outdir1
v1=-22
v2=38
region = ''

#def get_maps4serp(line, vtype=vtype, side=side, snrchoice=snrchoice, csize=csize, vsmo=vsmo, ssmo=ssmo, outdir = outdir, v1=v1, v2=v2):
def get_maps4serp(line, vtype=vtype, side=side, snrchoice=snrchoice, csize=csize, vsmo=vsmo, ssmo=ssmo, outdir = outdir, v1=v1, v2=v2):
#     line = input("Enter line: '13CO', 'C18O', 'N2HP', 'HCOP', 'CS', 'SO', 'H13COP', 'NH2D', 'All' : ")
# #    ctype = input("Enter type of cell, matched or org : ")
#     vtype = input("Enter type of a file, v10 or v06 : ")
#     side = input("both, left, or right:")
#     snrchoice = input("snr=")
#     vsmo = input("number of channels for the velocity smoothing = ")
#     ssmo = input("number of channels for the spatical smoothing = ")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if line != 'All':
        if vtype == 'v10':
            csize = '200' # the number of channels
            outdir_v10 = outdir+'v10/'
            if not os.path.exists(outdir_v10):
                os.mkdir(outdir_v10)
            for i in range(2):
                if i == 0: # for the matched
                    cubename = matched_cubedir+'SerB_'+line+'_v10_match_cube.fits'
                    apdix = 'match_snr'+snrchoice+'_vsmo'+vsmo+'_ssmo'+ssmo

                elif i == 1: # for the org
                    cubename = org_cubedir+'SerB_'+line+'_v10_cube.fits'
                    apdix = 'snr'+snrchoice+'_vsmo'+vsmo+'_ssmo'+ssmo

                #print(i)
                #print(apdix)
                a = Cube2map(cubename, rmssize=float(csize), getrms=side, snr=float(snrchoice), velocity_smo=int(vsmo), spatial_smo=int(ssmo))
                save_fits(outdir_v10+'SerB_'+line+'_v10_'+apdix+'_srms.fits', a.srms, a.header2d, overwrite=True)
                #set the default velocity range of the interst from -2 to 18 km/s
                save_fits(outdir_v10+'SerB_'+line+'_v10_'+apdix+'_m0.fits', a.moment0(vr=[v1,v2]), a.header2d, overwrite=True)
                save_fits(outdir_v10+'SerB_'+line+'_v10_'+apdix+'_m1.fits', a.moment1(vr=[v1,v2]), a.header2d, overwrite=True)

        elif vtype == 'v06':
            csize = '300' # the number  channels
            outdir_v06 = outdir+'v06/'
            if not os.path.exists(outdir_v06):
                os.mkdir(outdir_v06)
            for i in range(2):
                if i == 0: # for the matched
                    cubename = matched_cubedir+'SerB_'+line+'_match_cube.fits'
                    apdix = 'match_snr'+snrchoice+'_vsmo'+vsmo+'_ssmo'+ssmo

                elif i == 1: # for the org
                    cubename = org_cubedir+'SerB_'+line+'_cube.fits'
                    apdix = 'snr'+snrchoice+'_vsmo'+vsmo+'_ssmo'+ssmo

                #print(i)
                #print(apdix)
                a = Cube2map(cubename, rmssize=float(csize), getrms=side, snr=float(snrchoice), velocity_smo=int(vsmo), spatial_smo=int(ssmo))
                save_fits(outdir_v06+'SerB_'+line+'_'+apdix+'_srms.fits', a.srms, a.header2d, overwrite=True)
                #set the default velocity range of the interst from -2 to 18 km/s
                save_fits(outdir_v06+'SerB_'+line+'_'+apdix+'_m0.fits', a.moment0(vr=[v1,v2]), a.header2d, overwrite=True)
                save_fits(outdir_v06+'SerB_'+line+'_'+apdix+'_m1.fits', a.moment1(vr=[v1,v2]), a.header2d, overwrite=True)


###############################################################################
# auto completion of the directory of a file name at which a use want to take a look
###############################################################################
def get_filename(line, vtype=vtype, ctype=ctype):
#     line = input("Enter line: '13CO', 'C18O', 'N2HP', 'HCOP', 'CS', 'SO', 'H13COP', 'NH2D', 'All' : ")
#     vtype = input("Enter the velocity resolution, 0.06 or 0.1 km/s: ")
#     ctype = input("Enter type of cell, matched or org : ")
    if vtype == 'v06':
        vtype = ''
    elif vtype == 'v10':
        vtype = '_v10'
    if ctype == 'matched':
        cubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'
        appendix = '_match_cube.fits'
    elif ctype == 'org':
        cubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/org_cube/'
        appendix = '_cube.fits'


    filename = cubedir+'SerB_'+line+vtype+appendix

    return filename


def get_crop_filename(line, vtype=vtype, region=region):
#     line = input("Enter line: '13CO', 'C18O', 'N2HP', 'HCOP', 'CS', 'SO', 'H13COP', 'NH2D', 'All' : ")
#     vtype = input("Enter the velocity resolution, 0.06 or 0.1 km/s: ")
#     ctype = input("Enter type of map, full or cropped : ")
#     region = input("Enter region: SerB_, SerB_A_, or SerB_B_")
    if vtype == 'v06':
        vtype = ''
    elif vtype == 'v10':
        vtype = '_v10'
    if region == 'SerB_':
        cubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'
    elif region == 'SerB_SF_':
        cubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropSerFil/FUNS/'
    else:
        cubedir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'

    filename = cubedir+region+line+vtype+'_match_cube.fits'

    return filename


###############################################################################
#take from https://thispointer.com/python-how-to-add-append-key-value-pairs-in-dictionary-using-dict-update/
###############################################################################

def append_value(dict_obj, key, value):
    # Check if key exist in dict or not
    if key in dict_obj:
        # Key exist in dict.
        # Check if type of value of key is list or not
        if not isinstance(dict_obj[key], list):
            # If type is not list then make it list
            dict_obj[key] = [dict_obj[key]]
        # Append the value in list
        dict_obj[key].append(value)
    else:
        # As key is not in dict,
        # so, add key-value pair
        dict_obj[key] = value


###############################################################################
# auto completion of the directory of a file name at which a use want to take a look
###############################################################################
def get_maxmin_maps(inputcubes):
    dic_maxs ={'rms': [ ], \
            'm0' : [ ], \
            'm1' : [ ], \
            'm2' : [ ],
            }

    dic_mins ={'rms': [ ], \
            'm0' : [ ], \
            'm1' : [ ], \
            'm2' : [ ],
            }

    for i in range(len(inputcubes)):
        cube = inputcubes[i]
        srms = cube.srms
        m0 = cube.m0
        m1 = cube.m1
        m2 = cube.m2
        dic_maxs['rms'].append(np.nanmax(srms))
        dic_mins['rms'].append(np.nanmin(srms))
        dic_maxs['m0'].append(np.nanmax(m0))
        dic_mins['m0'].append(np.nanmin(m0))
        dic_maxs['m1'].append(np.nanmax(m1))
        dic_mins['m1'].append(np.nanmin(m1))
        dic_maxs['m2'].append(np.nanmax(m2))
        dic_mins['m2'].append(np.nanmin(m2))
        # append_value(dic_maxs, 'rms', np.nanmax(srms))
        # append_value(dic_mins, 'rms', np.nanmin(srms))
        # append_value(dic_maxs, 'm0', np.nanmax(m0))
        # append_value(dic_mins, 'm0', np.nanmin(m0))
        # append_value(dic_maxs, 'm1', np.nanmax(m1))
        # append_value(dic_mins, 'm1', np.nanmin(m1))
        # append_value(dic_maxs, 'm2', np.nanmax(m2))
        # append_value(dic_mins, 'm2', np.nanmin(m2))

    return dic_maxs, dic_mins



###############################################################################
# Line Profile
# Modification of funstools/plot.py full_line_scan
#full_line_scan --> full_line_scan_vtype : to use the function full_line_scan easily for SerB data
###############################################################################


def full_line_scan_vtype(loc='', vr=None, yi=1.5, cut=None, ver='otfpro', vtype=''):
    """
    Interactive line scan for full data set of FUNS project.
    This recognizes only the 'fullcube' fits files in stored in the 'loc' directory.
    The C18O file is required and the remaining lines are optional.
    Parameters:
        loc : path
            Directory path stored fullcube fits files.
        vr : list or tuple, optional
            Velocity (x-axis) range of the plot for line profile.
        yi : float
            Temperature (y-axis) interval to plot multi-lines.
        cut : list [ra_start, dec_start, ra_size, dec_size], optional
            Position and size of map cut in axis fraction.
        vtyp : to designate the file name according to the velocity resolution, v10 or v06. default value is v10 (added by KHK)
    Returns:
        matplotlib.figure, matplotlib.axes, matplotlib.axes
            Tuple with 3 components.
            You can modify the detail elements to save and print the figure.
    """
    if ver == 'otfpro':
        release = False
    elif ver == 'release':
        release = True
    else:
        raise ValueError("ver = {}, full_line_scan only work for 'otfpro' or 'release' version.")
    if not len(loc) == 0:
        if not loc[-1] == '/':
            loc = loc+'/'

# Added by KHK
    if not len(vtype) == 0:
        vtype = vtype
    else:
        vtype = 'v10'
##

    if release:
        base = glob(loc+'*C18O*_match_cube.fits')
    else:
        base = glob(loc+'*C18O*_fullcube.fits')
    if len(base) == 0:
        raise IOError('{}: No such file.'.format(base))
    else:
        args = (base[0].split('/')[-1]).split('_')
    files = []
    # lines = ['13CO', 'C18O', 'CS', 'HCOP', 'N2HP', 'SO', 'NH2D', 'H13COP']
    dcw = [1, 0, 1, 0, 1, 1, 0, 0]
    acw = [0, 1, 0, 1, 0, 0, 1, 1]
    ors = [200, 330]
    lines = ['H13COP', 'NH2D', 'N2HP', 'SO', 'HCOP', 'CS', 'C18O', '13CO']
    tlin = [r'$\mathrm{H^{13}CO^+}$', r'$\mathrm{NH_2D}$', r'$\mathrm{N_2H^+}$', 'SO', r'$\mathrm{HCO^+}$', 'CS', r'$\mathrm{C^{18}O}$', r'$\mathrm{^{13}CO}$']
    lnam = []
    tnam = []
    rs = []
    for ln, lin in enumerate(lines):
        if release:
            #nam = [loc+args[0]+'_'+lin+'_v10_match_cube.fits', loc+args[0]+'_'+lin+'_match_cube.fits']
            nam = [loc+args[0]+'_'+lin+'_'+vtype+'_match_cube.fits', loc+args[0]+'_'+lin+'_match_cube.fits']
        else:
            nam = [loc+args[0]+'_'+lin+'_all_'+args[3]+'_v10_fullcube.fits', loc+args[0]+'_'+lin+'_all_'+args[3]+'_fullcube.fits']
        if os.path.exists(nam[dcw[ln]]):
            files.append(nam[dcw[ln]])
            rs.append(ors[dcw[ln]])
            lnam.append(lin)
            tnam.append(tlin[ln])
        elif os.path.exists(nam[acw[ln]]):
            files.append(nam[acw[ln]])
            rs.append(ors[acw[ln]])
            lnam.append(lin)
            tnam.append(tlin[ln])
        else:
            print('{}\n{}\n: No such files.'.format(nam[0], nam[1]))
    if 'C18O' in lnam:
        li = lnam.index('C18O')
        l0 = Cube2map(files[li], rmssize=rs[li], velocity_smo=3, spatial_smo=2)
    else:
        l0 = Cube2map(files[-1], rmssize=rs[-1], velocity_smo=3, spatial_smo=2)
    lset = []
    for ln, lin in enumerate(files):
        lset.append(Cube2map(lin, rmssize=rs[ln], smoothing=False, masking=False))

    plt.ion()
    x = l0.x
    if vr is None:
        vr = [x[rs[-1]], x[-rs[-1]]]
    elif isinstance(vr, list) or isinstance(vr, tuple):
        pass
    else:
        raise TypeError("'vr' is not list or tuple.".format(vr))
    if cut is None:
        cut = [0., 0., 1., 1.]
    elif isinstance(cut, list) and len(cut) == 4:
        pass
    else:
        raise TypeError("'cut' is not readable form of [xs, ys, xl, yl].")
    # cut data
    rcut = [round(l0.nr*cut[0]), round(l0.nr*(cut[0]+cut[2]))]
    dcut = [round(l0.nd*cut[1]), round(l0.nd*(cut[1]+cut[3]))]
    cdat = l0.data[:, dcut[0]:dcut[1], rcut[0]:rcut[1]]
    crms = l0.rms[dcut[0]:dcut[1], rcut[0]:rcut[1]]
    # make map
    cdet = np.full(cdat.shape, 0.)
    cdet[np.nan_to_num(cdat) > 3*np.nanmedian(crms)] = 1.
    csum = np.nansum(cdet, axis=(1, 2))
    imap = np.sum(cdat[csum > 2*np.mean(csum)], axis=0)
    # map drawing
    lsfig = plt.figure(figsize=(16, 9))
    # lsfig.subplots_adjust(0.1, 0.1, 0.95, 0.95, 0.2, 0.2)
    lsmap = lsfig.add_axes([0.06, 0.07, 0.55, 0.9], projection=l0.wcs2d)
    lsmap.set_xlabel('R.A.')
    lsmap.set_ylabel('Dec.')
    lsmap.imshow(imap, origin='lower', vmin=0, cmap='inferno')
    cp = lsmap.scatter(-100, -100, s=60, color='lime')
    lsmap.set_xlim(rcut[0]-0.5, rcut[1]-0.5)
    lsmap.set_ylim(dcut[0]-0.5, dcut[1]-0.5)
    # line plot
    # lslin = lsfig.add_subplot(122)
    lslin = lsfig.add_axes([0.67, 0.07, 0.28, 0.9])
    if lnam[-1] == '13CO':
        ymax = 0.5*np.nanmax(lset[-1].data)
    else:
        ymax = max([np.nanmax(i.data) for i in lset])
    yr = [-4*np.nanmedian(crms), yi*(len(lset)-1)+ymax+np.nanmedian(crms)]
    print(ymax, yr)
    wset = [i.wcs2d for i in lset]

    def _onclick(event):
        if event.inaxes != lsmap.axes: return
        ri = int(round(event.xdata))
        di = int(round(event.ydata))
        if not release:
            rd = l0.wcs2d.wcs_pix2world(ri, di, 0)
            px = [i.wcs_world2pix(rd[0], rd[1], 0) for i in wset]
        lslin.clear()
        lslin.set_xlabel(r'$V_\mathrm{LSR}$ [km/s]')
        lslin.set_ylabel(r'$T_\mathrm{A}^* [K]$')
        lslin.set_xlim(*vr)
        lslin.set_ylim(*yr)
        lslin.set_yticks(np.arange(0, yr[1], yi))
        for li, ld in enumerate(lset):
            if release:
                r, d = ri*1, di*1
            else:
                r = int(round(px[li][0].item()))
                d = int(round(px[li][1].item()))
            lslin.plot([x[0], x[-1]], [li*yi, li*yi], lw=1, alpha=0.5, color='black')
            # lslin.plot([x[0], x[-1]], [3*l0.srms[d, r], 3*l0.srms[d, r]], ls='dotted', lw=1, alpha=0.3, color='blue')
            # lslin.plot([x[0], x[-1]], [3*ld.rms[d, r]+li*yi, 3*ld.rms[d, r]+li*yi], ls='dotted', lw=1, alpha=0.5, color='blue')
            # lslin.step(x, l0.sdata[:, d, r], color='blue', lw=1, alpha=0.6)
            if (0 <= r < ld.nr) and (0 <= d < ld.nd):
                cp.set_offsets([r, d])
                if lnam[li] == 'N2HP':
                    lslin.step(ld.x+8., ld.data[:, d, r]+li*yi, color='red', lw=1)
                elif lnam[li] == '13CO':
                    lslin.step(ld.x, 0.5*ld.data[:, d, r]+li*yi, color='red', lw=1)
                    lslin.annotate(r'$(\times 0.5)$', (0.97, ((li+0.3)*yi-yr[0])/(yr[1]-yr[0])), xycoords='axes fraction', ha='right')
                else:
                    lslin.step(ld.x, ld.data[:, d, r]+li*yi, color='red', lw=1)
            else:
                pass
            lslin.annotate('[{}, {}]'.format(d, r), (0.03, ((li+0.3)*yi-yr[0])/(yr[1]-yr[0])), xycoords='axes fraction')
            lslin.annotate(tnam[li], (1.01, (li*yi-yr[0])/(yr[1]-yr[0])), xycoords='axes fraction')
        c = l0.wcs2d.pixel_to_world(ri, di)
        lslin.annotate('Coord.: {}'.format(c.to_string('hmsdms')), (0.03, 0.95), xycoords='axes fraction')
        lslin.figure.canvas.draw()

    lsfig.canvas.mpl_connect('button_press_event', _onclick)
    return lsfig, lsmap, lslin











##############################
### FoF-0.2 dist_v0 = 0.2 dist_rd0=2
### FoF (adopting Shinyoung's code but use snr istead of area) - modify to keep the previously assigned filament inforamtion
### modify to adjust dist_v0 for the components with large sd (> 0.5km/s) & snr >5 ?
##############################
def fof_seq(siglim=[2.5], maxvel=[0.06], maxdist=2):


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

    dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')
    dc_data_file = 'dc_results_combined_revise_4_w5-13.dat'
    dc_data = dc_dir+'/'+dc_data_file

    sigma_init=10
    sigma_cut_array = sigmalim
    dist_v0_array = maxvel
    dist_rd0 = maxdist
    for sigma_cut in sigma_cut_array:
        for dist_v0 in dist_v0_array:

            comps=read(dc_data)
            comps['fn'] = 0
            comps['ii'] = 0
            #comps['in'] = 0

            step_cut1=sigma_cut
            snr=sigma_init



            rms_area = (comps['tp']/comps['snr'])*np.sqrt(2.*np.pi*comps['sd']**2.)

            ifil = 1
            ii = 0
            snr_index = 0
            while snr >= step_cut1:
                cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
                #cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
                ok_index = np.where(cmask==1)
                if np.sum(cmask) == 0:
                    snr -= 0.2
                    print('snr:', snr)
                    snr_index += 1
                    continue
                ii += 1
                icomps = np.where(cmask)[0]
                for i in icomps:
                    di, ri = comps['dp'][i], comps['rp'][i]
                    vi = comps['vp'][i]
                    seps = np.sqrt((comps['dp']-di)**2+(comps['rp']-ri)**2)
                    vdif = np.abs(comps['vp']-vi)
                    fmask = (vdif < dist_v0) & (seps < dist_rd0) & (comps['fn'] > 0)
                    if np.sum(fmask) == 0:
                        comps['fn'][i] = ifil
                        #comps['in'][i] = k
                        comps['ii'][i] = ii
                        ifil += 1
                    else:
                        #k += 1
                        fils = np.unique(comps['fn'][fmask])
                        comps['fn'][i] = fils[0]
                        #comps['in'][i] = k
                        comps['ii'][i] = ii
                        if len(fils) > 1:
                            for fn in fils[1:]:
                                gmask = comps['fn'] == fn
                                comps['fn'][gmask] = fils[0]

                remains = np.sum(comps['fn'] == 0)
                print('snr: {},  ifil: {}, ii: {}, remains: {}'.format(snr, ifil, ii, remains))


            FoF_save = dc_dir+'/ff_sequence/ff_keep/'
        #    write(comps[comps['fn']>0], FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat', overwrite=True)
            write(comps, FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat', overwrite=True)





# ###############################################################################
# # Line profiles for n by n pixesl
# # Modification of funstools/Cube2map.line_scan()
# #line_scan --> line_scan_multi
# # compare 13CO, C18O, and N2HP
# ###############################################################################
# def line_scan_multi(self, vr=None, yr=None):
#     """
#     Interactive scan for line profiles.
#     Plot line profile on click position.
#     Parameters:
#         vr : list or tuple, optional
#             Range of velocity axis for line plot.
#         yr : list or tuple, optional
#             Range of temperature axis for line plot.
#         cut : list [ra_start, dec_start, ra_size, dec_size], optional
#             Position and size of map cut in axis fraction.
#     Returns:
#         tuple (plt.figure, plt.subplot, plt.subplot)
#     """
#     plt.ion()
#     plt.rcParams['xtick.direction'] = 'in'
#     plt.rcParams['ytick.direction'] = 'in'
#     x = self.x
#     if vr is None:
#         vr = [x[self._rmssize], x[-self._rmssize]]
#     elif isinstance(vr, list) or isinstance(vr, tuple):
#         pass
#     else:
#         raise TypeError("'vr' is not list or tuple.".format(vr))
#     if yr is None:
#         yr = [-4*np.nanmedian(self.rms), np.nanmax(self.data)+np.nanmedian(self.rms)]
#     elif isinstance(yr, list) or isinstance(yr, tuple):
#         pass
#     else:
#         raise TypeError("'yr' is not list or tuple".format(yr))
#     # make map
#     det = np.full(self.data.shape, 0.)
#     det[np.nan_to_num(self.data) > 3*np.nanmedian(self.rms)] = 1.
#     dsum = np.nansum(det, axis=(1, 2))
#     imap = np.sum(self.data[dsum > 2*np.mean(dsum)], axis=0)
#     # map drawing
#     lsfig = plt.figure(figsize=(16, 9))
#     lsmap = lsfig.add_axes([0.06, 0.07, 0.36, 0.9], projection=self.wcs2d)
#     lsmap.set_xlabel('R.A.')
#     lsmap.set_ylabel('Dec.')
#     lsmap.imshow(imap, origin='lower', vmin=0, cmap='inferno')
#     self._pp = lsmap.scatter(-100, -100, s=60, color='lime')
#     lsmap.set_xlim(-0.5, self.nr-0.5)
#     lsmap.set_ylim(-0.5, self.nd-0.5)
#     # line plot
#     lbl = lsfig.add_axes([0.47, 0.07, 0.17, 0.3])
#     lcl = lsfig.add_axes([0.47, 0.37, 0.17, 0.3])
#     ltl = lsfig.add_axes([0.47, 0.67, 0.17, 0.3])
#     lbc = lsfig.add_axes([0.64, 0.07, 0.17, 0.3])
#     lslin = lsfig.add_axes([0.64, 0.37, 0.17, 0.3])
#     ltc = lsfig.add_axes([0.64, 0.67, 0.17, 0.3])
#     lbr = lsfig.add_axes([0.81, 0.07, 0.17, 0.3])
#     lcr = lsfig.add_axes([0.81, 0.37, 0.17, 0.3])
#     ltr = lsfig.add_axes([0.81, 0.67, 0.17, 0.3])
#
#     def _onclick(event):
#         if event.inaxes == lsmap.axes:
#             self._pr = int(round(event.xdata))
#             self._pd = int(round(event.ydata))
#         elif event.inaxes == lbl.axes:
#             self._pr -= 1
#             self._pd -= 1
#         elif event.inaxes == lcl.axes:
#             self._pr -= 1
#         elif event.inaxes == ltl.axes:
#             self._pr -= 1
#             self._pd += 1
#         elif event.inaxes == lbc.axes:
#             self._pd -= 1
#         elif event.inaxes == ltc.axes:
#             self._pd += 1
#         elif event.inaxes == lbr.axes:
#             self._pr += 1
#             self._pd -= 1
#         elif event.inaxes == lcr.axes:
#             self._pr += 1
#         elif event.inaxes == ltr.axes:
#             self._pr += 1
#             self._pd += 1
#         else:
#             return
#         r = self._pr
#         d = self._pd
#         ri = [-1, -1, -1, 0, 0, 0, +1, +1, +1]
#         di = [-1, 0, +1, -1, 0, +1, -1, 0, +1]
#         for i, ax in enumerate([lbl, lcl, ltl, lbc, lslin, ltc, lbr, lcr, ltr]):
#             ax.clear()
#             ax.set_xlim(*vr)
#             ax.set_ylim(*yr)
#             if not i in [0, 1, 2]:
#                 ax.set_yticklabels([])
#             if not i in [0, 3, 6]:
#                 ax.set_xticklabels([])
#             rr = r+ri[i]
#             dd = d+di[i]
#             if 0 <= rr < self.nr and 0 <= dd < self.nd:
#                 self._pp.set_offsets([rr-1, dd-1])
#                 srms = self._snr*self.srms[dd, rr]
#                 rms = self._snr*self.rms[dd, rr]
#                 ax.plot([x[0], x[-1]], [0, 0], lw=1, alpha=0.5, color='black')
#                 ax.plot([x[0], x[-1]], [srms, srms], ls='dotted', lw=1, alpha=0.3, color='blue')
#                 ax.plot([x[0], x[-1]], [rms, rms], ls='dotted', lw=1, alpha=0.5, color='red')
#                 ax.step(x, self.sdata[:, dd, rr], color='blue', lw=1, alpha=0.6)
#                 ax.step(x, self.data[:, dd, rr], color='red', lw=1)
#
#         lbc.set_xlabel(r'$V_\mathrm{LSR}$ [km/s]')
#         lcl.set_ylabel(r'$T_\mathrm{A}^* [K]$')
#         c = self.wcs2d.pixel_to_world(r, d)
#         lslin.annotate('Coord.: {}'.format(c.to_string('hmsdms')), (0.03, 0.92), xycoords='axes fraction')
#         lslin.annotate('Pixel: [{}, {}]'.format(d, r), (0.03, 0.86), xycoords='axes fraction')
#         lslin.figure.canvas.draw()
#
#     lsfig.canvas.mpl_connect('button_press_event', _onclick)
#     return lsfig, lsmap, lslin
#

# ###############################################################################
# # Decomposing
# ###############################################################################
# snr=3
# rmssize=200
# ssm_find=5
# vsm_find=5
# ssm_fit=1
# vsm_fit=1
#
# def SerB_decompose(vtype='', snr=snr, rmssize=rmssize, vsm_find=vsm_find, ssm_find=ssm_find, vsm_fit=vsm_fit, ssm_fit=ssm_fit):
#     """
#     to use funstools.Decompose easily for decomposing multiple Gaussian components
#     from C18O line cube data of FUNS SerB data.
#     """
#
#     path = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'
#     wpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/'
#     if not len(vtype) == 0:
#         vtype = vtype
#     else:
#         vtype = 'v10'
#
#     if vtype == 'v06':
#         vtail = 'match_cube.fits'
#     else:
#         vtail = vtype+'_match_cube.fits'
#
#
#     filename = path+'SerB_C18O_'+vtail
#
#     wdirname = str(ssm_find)+str(vsm_find)+str(ssm_fit)+str(vsm_fit)+vtype+'snr'+str(snr)
#
#     if not os.path.exists(wpath+wdirname):
#         os.mkdir(wpath+wdirname)
#
#     cube = Decompose(filename, snr=snr, rmssize=rmssize, spasmo_find=ssm_find, velsmo_find=vsm_find, spasmo_fit = ssm_fit, velsmo_fit=vsm_fit)
#     cube.run_decompose(save=wpath+wdirname+'/SerB_C18O_decompose_result_'+wdirname+'.dat')
#



############################
# get_rms
############################
path='/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'
cube_C18O_v10 = path+'SerB_C18O_v10_match_cube.fits'
cube_C18O_v06 = path+'SerB_C18O_match_cube.fits'
cube_13CO_v10 = path+'SerB_13CO_v10_match_cube.fits'
cube_13CO_v06 = path+'SerB_13CO_match_cube.fits'
cube_N2HP_v10 = path+'SerB_N2HP_v10_match_cube.fits'
cube_N2HP_v06 = path+'SerB_N2HP_match_cube.fits'
cube_CS_v10 = path+'SerB_CS_v10_match_cube.fits'
cube_CS_v06 = path+'SerB_CS_match_cube.fits'
cube_HCOP_v10 = path+'SerB_HCOP_v10_match_cube.fits'
cube_HCOP_v06 = path+'SerB_HCOP_match_cube.fits'
cube_SO_v10 = path+'SerB_SO_v10_match_cube.fits'
cube_SO_v06 = path+'SerB_SO_match_cube.fits'
cube_H13COP_v10 = path+'SerB_H13COP_v10_match_cube.fits'
cube_H13COP_v06 = path+'SerB_H13COP_match_cube.fits'
cube_NH2D_v10 = path+'SerB_NH2D_v10_match_cube.fits'
cube_NH2D_v06 = path+'SerB_NH2D_match_cube.fits'

cube = cube_NH2D_v06
rms = get_rms(cube, where='both', size = 200)
np.nanmean(rms)
np.nanstd(rms)
