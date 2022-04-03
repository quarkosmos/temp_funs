from funstools import Cube2map
from funstools import full_line_scan
#from funstools import Decompose
from funstools import save_fits
from multiaxes import Multiaxes

from SerB_Cube import SerBCube2map
from SerB_PlotFitRes_FF import SerBDecompose
from run_funstools_Serp import get_crop_filename

from FUNS_FoF_tools import FUNS_FoF as FF
import os
from astropy.io.ascii import read, write
import pickle

m0file= '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/A/2521v06snr3mww0.060.034.0/SerB_A_C18O21v06snr3_m0_v0313.fits'
filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/SerB_A_C18O_match_cube.fits'
fitres3 = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/A/2521v06snr3mww0.060.034.0/SerB_A_C18O_fitres3_dc3.dat'

cube = SerBDecompose(filename)
dc = read(fitres3)
cube.plot_fit(dc, vr=(4, 12), n_ch=3, nzp=9)

cube = Cube2map(filename)
cube.line_scan(vr=(4,12))



#############
## gpy finalized data  to the funstools decomposed data table format

from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os


gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
dc_option = 'sn2.5d3s1'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option,'output', 'gpy_decomposed')
gpy_final = 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized.dat'.format(dc_option)
gpy2fitres = 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
#        print(j)
    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)

xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, amp, VLSR, vel_disp, snr, dv, int_tot], names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'))
fr.write(os.path.join(path_to_decomp_data, gpy2fitres), format='ascii', overwrite='True')




###############################
## 신영씨의 FoF
while snr > step_cut1:
    cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    if np.sum(cmask) == 0:
        snr -= 0.1
        print('snr:', snr)
        continue
    ii += 1
    icomps = np.where(cmask)[0]
    for i in icomps:
        di, ri = comps['dp'][i], comps['rp'][i]
        vi = comps['vp'][i]
        seps = np.sqrt((comps['dp']-di)**2+(comps['rp']-ri)**2)
        vdif = np.abs(comps['vp']-vi)
        fmask = (vdif < max_vel) & (seps < max_sep) & (comps['fn'] > 0)
        if np.sum(fmask) == 0:
            comps['fn'][i] = new_fil
            comps['ii'][i] = ii
            new_fil += 1
        else:
            fils = np.unique(comps['fn'][fmask])
            comps['fn'][i] = fils[0]
            comps['ii'][i] = ii
            if len(fils) > 1:
                for fn in fils[1:]:
                    gmask = comps['fn'] == fn
                    comps['fn'][gmask] = fils[0]
    remains = np.sum(comps['fn'] == 0)
    print('snr: {}, ii: {}, new_fil: {}, remains: {}'.format(snr, ii, new_fil, remains))





#######################################
## GPY FoF
from FUNS_FoF_tools import FUNS_FoF as FF
#from FUNS_FoF import run_FoF
#from FUNS_FoF import FoF_3D_plot
#from FUNS_FoF import FoF_line_profiles
#from FUNS_FoF import Fil_in_2D


dckey = 'sn2.5d2s1'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
snr=2.5
ssm_find=2
vsm_find=2
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.03
rd0=2.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_{}_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()

a.fof_3D_plot()


a.dcstat()



from FUNS_FoF_tools import FUNS_FoF as FF
dckey = 'sn2.5d2s2'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
snr=2.5
ssm_find=2
vsm_find=2
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.06
rd0=3.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()


from FUNS_FoF_tools import FUNS_FoF as FF
dckey = 'sn2.5d3s1'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
snr=2.5
ssm_find=2
vsm_find=2
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.06
rd0=3.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()


from FUNS_FoF_tools import FUNS_FoF as FF
dckey = 'sn2.5d3s2'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
snr=2.5
ssm_find=2
vsm_find=2
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.06
rd0=3.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()


from FUNS_FoF_tools import FUNS_FoF as FF
dckey = 'sn3.0d2s1'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
snr=3.0
ssm_find=2
vsm_find=2
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.06
rd0=3.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()


from FUNS_FoF_tools import FUNS_FoF as FF
dckey = 'sn3.0d2s2'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
snr=3.0
ssm_find=2
vsm_find=2
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.06
rd0=3.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()


from FUNS_FoF_tools import FUNS_FoF as FF
dckey = 'sn3.0d3s1'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
snr=3.0
ssm_find=2
vsm_find=2
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.06
rd0=3.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()



from FUNS_FoF_tools import FUNS_FoF as FF
dckey = 'sn3.0d3s2'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
snr=3.0
ssm_find=2
vsm_find=2
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.02
rd0=2.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()


###############reload
import importlib
import FUNS_FoF_tools
FUNS_FoF_tools = importlib.reload(FUNS_FoF_tools)
from FUNS_FoF_tools import FUNS_FoF as FF



############ Test Plot_FitRes code
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/SerB_A_C18O_match_cube.fits'
#wfilename = region+line+'_vsm'+str(self._velocity_smo)+'_ssm'+str(self._spatial_smo)+'_v06_cube.fits'
region = 'SerB_A_'
line = 'C18O'
dc_gpy_key = 'sn3.0d3s1'
a = PlotFitRes(input_cube, key=dc_gpy_key)

#smoothed file
b = a.save_smoothed_3dcube()
h=a.header
h['COMMENT'] = '= spatial smoothing by '+str(a._spatial_smo)
h['COMMENT'] = '= veolocity smoothing by '+str(a._velocity_smo)
wfilename = 'SerB_A_C18O_vsm'+str(a._velocity_smo)+'_ssm'+str(a._spatial_smo)+'_v06_cube.fits'
save_fits(wdir+wfilename, b, h, overwrite=True)

#gpy+_dc_result file
dc_gpy_key = 'sn3.0d3s1'
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/', dc_gpy_key, 'output', \
'gpy_decomposed')
dc_fin_filename = region+line+'_match_cube_'+dc_gpy_key+'_g+_fit_fin_sf-p2_finalized_fitres.dat'
dc_data_file = dc_dir+'/'+dc_fin_filename

fitres = read(dc_data_file)
a.plot_dc_fit(fitres, vr=(3,13), n_ch=5, nzp=15)

fof_file = dc_dir+'/fof/FoF_result_SerB_A_C18O_match_cube_'+dc_gpy_key+'_g+sig_3.0_dv0_0.06_rd0_1.5.dat'
dcresult = dc_data_file
a.plot_fof_profile(line=line, vr=(3,13), n_ch=5, nzp=15, vtype='v06', nbeam=5, dckey=dc_gpy_key, filsigmacut=3, dv0=0.06, drd0=1.5, fof_result=fof_file, dcresult=dcresult, dctype='gpy+')


import importlib
import Plot_FitRes_GPYFUNSTOOLs
Plot_FitRes_GPYFUNSTOOLs = importlib.reload(Plot_FitRes_GPYFUNSTOOLs)
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes

############






######## FUNSTOOLs  Decompose 8/07/21

#####################################################################
# Decomposing
import funstools
from funstools import get_rms
from funstools import Cube2map
from funstools import save_fits
from matplotlib import pyplot as plt
from funstools import full_line_scan
from funstools import Decompose
#from astroML.plotting import multiaxes
#from multiaxes import Multiaxes
from astropy.wcs import WCS
from run_funstools_Serp import get_filename
from run_funstools_Serp import get_maxmin_maps
from run_funstools_Serp import full_line_scan_vtype

#from SerB_Cube import SerBCube2map
#from SerB_Cube import Cube2map
from SerB_PlotFitRes_FF import SerBDecompose


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
from astropy.io.ascii import read

##### directory prep

funs_crop_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'
funs_full_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'

#dir = funs_crop_dir
dir = funs_full_dir

region_code_A = 'SerB_A_'
region_code_B = 'SerB_B_'
region_code_main = 'SerB_'

#region_code = region_code_A
#region_code = region_code_B
region_code = region_code_main

decom_wdir_A = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/A/'
decom_wdir_B = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/B/'
decom_wdir_main = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/Main/'

#decom_wdir = decom_wdir_A
#decom_wdir = decom_wdir_B
decom_wdir = decom_wdir_main

maps_out_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/figs/decompose_tests/'

funsmap_name = [r'$^{13}$CO',\
         r'C$^{18}$O',\
         r'N$_{2}$H$^{+}$', \
         r'HCO$^{+}$', \
         r'CS', \
         r'SO', \
         r'NH$_{2}$D', \
         r'H$^{13}$CO$^{+}$']

funslinename = ['13CO', 'C18O', 'N2HP', 'HCOP', 'CS', 'SO', 'NH2D', 'H13COP']

v10 = '_v10_match_cube.fits'
v06 = '_match_cube.fits'
#vtype = v10
vtype_file = v06
vtype = 'v06'

if vtype == 'v10':
    rs = 200
elif vtype == 'v06':
    rs = 300




#####################################################################
# Decomposing

decom_line = funslinename[1]

input_file = get_crop_filename(line=decom_line, vtype=vtype, region=region_code)

ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 2, 4, 2, 1, vtype, 3, 0.03, 0.03, 4.0 # 2421 0.03 0.03 4.0
# ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 2, 4, 2, 1, vtype, 3, 0.06, 0.06, 4.0 # 2421 0.06 0.06 4.0
# ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 1, 5, 1, 1, vtype, 3, 0.06, 0.03, 4.0 # 1511 0.06 0.03 4.0
# ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 2, 5, 2, 1, vtype, 3, 0.06, 0.03, 4.0 # 2421 0.06 0.03 4.0

decompose_key = str(ssm_find)+str(vsm_find)+str(ssm_fit)+str(vsm_fit)+vtype+'snr'+str(snr)+'mww'+str(mlim)+str(wmin)+str(wmax)
dc_filename_prefix = decom_wdir+decompose_key+'/'+region_code+decom_line

input_file = get_filename(line=decom_line, vtype='v06', ctype='matched')

# ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 2, 4, 2, 1, vtype, 3, 0.03, 0.03, 4.0 # 2421 0.03 0.03 4.0
# ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 2, 4, 2, 1, vtype, 3, 0.06, 0.06, 4.0 # 2421 0.06 0.06 4.0
# ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 1, 5, 1, 1, vtype, 3, 0.06, 0.03, 4.0 # 1511 0.06 0.03 4.0
# ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 2, 5, 2, 1, vtype, 3, 0.06, 0.03, 4.0 # 2421 0.06 0.03 4.0
ssm_find, vsm_find, ssm_fit, vsm_fit, vtype, snr, mlim, wmin, wmax = 4, 4, 1, 1, vtype, 5, 0.06, 0.06, 3.0  #4411 snr5 0.06 0.06 3.0

decompose_key = str(ssm_find)+str(vsm_find)+str(ssm_fit)+str(vsm_fit)+vtype+'snr'+str(snr)+'mww'+str(mlim)+str(wmin)+str(wmax)
dc_filename_prefix = decom_wdir+decompose_key+'/'+region_code+decom_line

cube = SerBDecompose(input_file, rmssize = rs, region=region_code, snr=snr, spasmo_find=ssm_find, velsmo_find=vsm_find, spasmo_fit=ssm_fit, velsmo_fit=vsm_fit, mlim=mlim, wmin=wmin, wmax=wmax)

fitres1 = read(dc_filename_prefix+'_fitres1_dc1.dat')
fitres2 = read(dc_filename_prefix+'_fitres2_dc2.dat')
fitres3 = read(dc_filename_prefix+'_fitres3_dc3.dat')
#cube.plot_fit(fitres3, vr=(6,11), n_ch=5, nzp=15)






############ FUNSTOOLs  FF #########################
#### run FoF with decomposed results by funstools

from FUNS_FoF_tools import FUNS_FoF as FF
dckey = '4411v06snr5mww0.060.063.0'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/Main/'+dckey+'/'
dcdata = dcpath+'SerB_C18O_fitres3_dc3.dat'
region='SerB_'
snr=3.0
ssm_find=2
vsm_find=4
ssm_fit=2
vsm_fit=1
sigma=3.0
v0=0.06
rd0=2.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_ft_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region=region, line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()


import importlib
import FUNS_FoF_tools
FUNS_FoF_tools = importlib.reload(FUNS_FoF_tools)
from FUNS_FoF_tools import FUNS_FoF as FF




############ Plot_FitRes : FUNSTOOLS decompose --> FoF
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/match_cube/SerB_C18O_match_cube.fits'
#wfilename = region+line+'_vsm'+str(self._velocity_smo)+'_ssm'+str(self._spatial_smo)+'_v06_cube.fits'
region = 'SerB_'
line = 'C18O'
sigma=3.0
v0=0.02
rd0=3.0
dc_ft_key = '4411v06snr5mww0.060.063.0'
a = PlotFitRes(input_cube, key=dc_ft_key)

# #smoothed file
# b = a.save_smoothed_3dcube()
# h=a.header
# h['COMMENT'] = '= spatial smoothing by '+str(a._spatial_smo)
# h['COMMENT'] = '= veolocity smoothing by '+str(a._velocity_smo)
# wfilename = 'SerB_A_C18O_vsm'+str(a._velocity_smo)+'_ssm'+str(a._spatial_smo)+'_v06_cube.fits'
# save_fits(wdir+wfilename, b, h, overwrite=True)

#gpy+_dc_result file
dc_ft_key = '4411v06snr5mww0.060.063.0'
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'decompose6_AB', 'Main', dc_ft_key)
dc_fin_filename = region+line+'_fitres3_dc3.dat'
dc_data_file = dc_dir+'/'+dc_fin_filename

fitres = read(dc_data_file)
a.plot_dc_fit(fitres, vr=(3,13), n_ch=5, nzp=15)

fof_file = dc_dir+'/fof/FoF_result_'+region+line+'_match_cube_'+dc_ft_key+'_ft_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'
dcresult = dc_data_file
a.plot_fof_profile(line=line, vr=(3,13), n_ch=5, nzp=15, vtype='v06', nbeam=5, dckey=dc_ft_key, filsigmacut=sigma, dv0=v0, drd0=rd0, fof_result=fof_file, dcresult=dcresult, dctype='funstools')


import importlib
import Plot_FitRes_GPYFUNSTOOLs
Plot_FitRes_GPYFUNSTOOLs = importlib.reload(Plot_FitRes_GPYFUNSTOOLs)
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes




#####################################################
############ Smoothing Data
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_match_cube.fits'
region = 'SerB_'
line = 'C18O'
snr=5.0
ssm_find=3
vsm_find=2
#sigma=3.0
#v0=0.03
#rd0=3.0
dc_gpy_key = 'sn5.0d2s1'
a = PlotFitRes(input_cube, key=dc_gpy_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)

#smoothed file
b = a.save_smoothed_3dcube()
h=a.header
h['COMMENT'] = '= spatial smoothing by '+str(a._spatial_smo)
h['COMMENT'] = '= veolocity smoothing by '+str(a._velocity_smo)
wfilename = 'SerB_C18O_vsm'+str(a._velocity_smo)+'_ssm'+str(a._spatial_smo)+'_v06_cube.fits'
save_fits(wdir+wfilename, b, h, overwrite=True)



############ Smoothing Data : Crop Main
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropMain/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropMain/SerB_C18O_match_cube_cut.fits'
region = 'SerB_'
line = 'C18O'
snr=5.0
ssm_find=1
vsm_find=2
#sigma=3.0
#v0=0.03
#rd0=3.0
dc_gpy_key = 'sn5.0d3s1'
a = PlotFitRes(input_cube, key=dc_gpy_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)

#smoothed file
b = a.save_smoothed_3dcube()
h=a.header
h['COMMENT'] = '= spatial smoothing by '+str(a._spatial_smo)
h['COMMENT'] = '= veolocity smoothing by '+str(a._velocity_smo)
wfilename = 'SerB_C18O_vsm'+str(a._velocity_smo)+'_ssm'+str(a._spatial_smo)+'_v06_cube_cut.fits'
save_fits(wdir+wfilename, b, h, overwrite=True)



############ Smoothing Data : rect_cube for GPY+
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/rect_cube/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/rect_cube/SerB_C18O_match_v06_cube_rect.fits'
region = 'SerB_'
line = 'C18O'
snr=5.0
ssm_find=1
vsm_find=1
#sigma=3.0
#v0=0.03
#rd0=3.0
dc_gpy_key = 'sn5.0d2s1'
a = PlotFitRes(input_cube, key=dc_gpy_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)

#smoothed file
b = a.save_smoothed_3dcube()
h=a.header
h['COMMENT'] = '= spatial smoothing by '+str(a._spatial_smo)
h['COMMENT'] = '= veolocity smoothing by '+str(a._velocity_smo)
wfilename = 'SerB_C18O_vsm'+str(a._velocity_smo)+'_ssm'+str(a._spatial_smo)+'_v06_cube_rect.fits'
save_fits(wdir+wfilename, b, h, overwrite=True)


#####################################################
####  PICKLE
import pickle
from gausspyplus.plotting import pickle_load_file
org = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/sn5.0d3s1_sv2521/output/gpy_prepared/SerB_C18O_match_cube.pickle'
smo1 = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/sn5.0d3s1_sv2521/output/gpy_prepared/SerB_C18O_vsm1.0_ssm2.0_v06_cube.pickle'
smo2 = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/sn5.0d3s1_sv2521/output/gpy_prepared/SerB_C18O_vsm5.0_ssm2.0_v06_cube.pickle'
step4 = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/sn5.0d3s1_sv2521/output/gpy_decomposed/SerB_C18O_vsm5.0_ssm2.0_v06_cube_sn5.0d3s1_g+_fit_fin.pickle'
step5 = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/sn5.0d3s1_sv2521/output/gpy_decomposed/SerB_C18O_vsm5.0_ssm2.0_v06_cube_sn5.0d3s1_g+_fit_fin_sf-p1.pickle'

data_org = pickle_load_file(org)
data_smo1 = pickle_load_file(smo1)
data_smo2 = pickle_load_file(smo2)
data_step4 = pickle_load_file(step4)
data_step5 = pickle_load_file(step5)

ymin, xmin = min(data['location'])



#################################################
#### beam convolution
import numpy as np
from funstools import Cube2map, save_fits, trao_beam, hgbs_beam
from astropy.io import fits
from astropy import units as u
from astropy.convolution.convolve import convolve
from astropy.convolution.kernels import Gaussian2DKernel
from FITS_tools.hcongrid import hcongrid

path = '/Volumes/DATA/shiny/Dropbox/Research/TRAO-FUNS/IC5146/'
path_c18o = path+'fits/IC5146_C18O_all_c20_v10_fullcube.fits'
path_hgbs = path+'HGBS_ic5146_'
maps = ['dust_temperature', 'column_density']
lines = ['C18O', 'N2HP']

cube_trao = Cube2map(path_c18o)
header_trao = cube_trao.header2d.copy()

for ln, line in enumerate(lines):
    for mn, map in enumerate(maps):
        hdu = fits.open(path_hgbs+map+'_map.fits')[0]
        pixel_size = hdu.header['CDELT2']*u.deg
        beam_hgbs = (hgbs_beam(500)/pixel_size).decompose()
        beam_trao = (trao_beam(line)/pixel_size).decompose()
        kernel_size = np.sqrt(beam_trao**2-beam_hgbs**2)/np.sqrt(8*np.log(2))
        kernel = Gaussian2DKernel(kernel_size)
        data = hdu.data.copy()
        data_convolved = convolve(data, kernel, 'extend', preserve_nan=True)
        header_old = hdu.header.copy()
        header_new = hdu.header.copy()
        for hnam in ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']:
            for hn in ['1', '2']:
                header_new[hnam+hn] = header_trao[hnam+hn]
        data_matched = hcongrid(data_convolved.copy(), header_old, header_new, True)
        save_fits(path_hgbs+map+'_matched_to_{}.fits'.format(line), data_matched, header_new, True)



######################################################
####### Sound speed, thermal velocity dispersion #####
"'To make the thermal velocity dispersion (sigma_th) map and isothermal sound speed map'"
import numpy as np
from funstools import Cube2map, save_fits, trao_beam, hgbs_beam
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
import math

line = 'C18O'
maps = ['Cs', 'Sig_th']
molw = [2.33, 30]
comm = ['isothermal_sound_speed (added by KHK)', 'thermal_velocity_dispersion (added by KHK)']

path = '/Users/khkim/Dropbox/Work5/FUNS-Serp/'
path_c18o = path+'datared_01/match_cube/SerB_C18O_match_cube.fits'
path_hgbs = path+'HGBS_cut/'
tdfitsfile = path_hgbs+'convolved/td_matched_to_{}.fits'.format(line)

for i, maps in enumerate(maps):
    hdu = fits.open(tdfitsfile)[0]
    data = hdu.data.copy()
    header_old = hdu.header.copy()
    header_new = hdu.header.copy()
    header_new['COMMENT'] = '{}'.format(comm[i])
    header_new['BUNIT'] = 'km/s'
    mol_weight = molw[i]
    cs2 = ((1*const.k_B * data*u.K)/(mol_weight *u.u)).to(u.km*u.km/u.s/u.s)
    cs1= (cs2.value)**(1/2)
    #cs1 = math.sqrt(cs2.value)
    cs = cs1
    save_fits(path_hgbs+'convolved/'+maps+'_matched_to_{}.fits'.format(line), cs, header_new, True)



#####################################################
##### draw maps of Cs and sigma_th

import funstools
from funstools import get_rms, get_mask, get_det, wcs2d, header2d, save_fits
from funstools import smooth3d, radsmo2d, boxsmo1d, make_velo
from funstools import Cube2map
import numpy as np
from astropy.table import Table, Column
from astropy.io.ascii import read, write
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle, Circle
import matplotlib.colors as mcolors
import matplotlib as mpl
import os
import csv
import re
from astropy.io.ascii import read

from astropy.wcs import WCS
from run_funstools_Serp import get_filename
from run_funstools_Serp import get_maxmin_maps

from astropy.io import fits



path = '/Users/khkim/Dropbox/Work5/FUNS-Serp/'
path_c18o = path+'datared_01/match_cube/SerB_C18O_match_cube.fits'
path_hgbs = path+'HGBS_cut/'

line = 'C18O'
maps = ['td', 'Cs', 'Sig_th']
titles = ['Herschel Dust Temperature', 'Isothermal Sound Speed', 'Thermal velocity dispersion']
ylabel = [r'T$_d$ (K)', r'C$_s$ (km/s)', r'$\sigma_{th}$ (km/s)']


for i, maps in enumerate(maps):
    fitsfile = path_hgbs+'convolved/{}_matched_to_{}.fits'.format(maps, line)
    hdu = fits.open(fitsfile)[0]
    data = hdu.data.copy()

    cnfig = plt.figure(figsize=(5, 8))
    ax = cnfig.add_axes([0.12, 0.02, 0.62, 0.5/(118/235)])
    #im = ax.imshow(m0, origin='lower', vmin=0, cmap='gist_yarg')
    #ax.xaxis.set_tick_params(labelsize=12)
    #ax.yaxis.set_tick_params(labelsize=12)
    ax.set_xlabel('RA (pixel num)', fontsize=10)
    ax.set_ylabel('Dec (pixel num)', fontsize=10)
    ax.set_title(titles[i])

    im_cn = ax.imshow(data, origin='lower', vmin= np.nanmin(data), vmax=np.nanmax(data), cmap=plt.cm.rainbow, alpha=0.4)
    cax = cnfig.add_axes([0.8, 0.15, 0.02, 0.72])
    cbar = cnfig.colorbar(im_cn, cax=cax)
    cax.set_ylabel(ylabel[i], fontsize=10, labelpad=10)


    #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.85, wspace=0.15, hspace=0.13)
    plt.savefig(path_hgbs+maps+'_convolved_map.pdf')


    plt.close()
    ### mo 그리고, 그 위에 color code로 number of components map overplot






#####################################################################
# Decomposing
import funstools
from funstools import get_rms
from funstools import Cube2map
from funstools import save_fits
from matplotlib import pyplot as plt
from funstools import full_line_scan
from funstools import Decompose
#from astroML.plotting import multiaxes
#from multiaxes import Multiaxes
from astropy.wcs import WCS
from run_funstools_Serp import get_filename
from run_funstools_Serp import get_maxmin_maps
from run_funstools_Serp import full_line_scan_vtype

#from SerB_Cube import SerBCube2map as Cube2map
#from SerB_Cube import Cube2map
from SerB_PlotFitRes_FF_faff import SerBDecompose
#from SerB_PlotFitRes import SerBDecompose


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
from astropy.io.ascii import read

##### directory prep

funs_crop_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'
funs_full_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'

#dir = funs_crop_dir
dir = funs_full_dir

region_code_A = 'SerB_A_'
region_code_B = 'SerB_B_'
region_code_main = 'SerB_'

#region_code = region_code_A
#region_code = region_code_B
region_code = region_code_main

decom_wdir_A = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/A/'
decom_wdir_B = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/B/'
decom_wdir_main = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/Main/'

#decom_wdir = decom_wdir_A
#decom_wdir = decom_wdir_B
decom_wdir = decom_wdir_main

maps_out_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/figs/decompose_tests/'

funsmap_name = [r'$^{13}$CO',\
         r'C$^{18}$O',\
         r'N$_{2}$H$^{+}$', \
         r'HCO$^{+}$', \
         r'CS', \
         r'SO', \
         r'NH$_{2}$D', \
         r'H$^{13}$CO$^{+}$']

funslinename = ['13CO', 'C18O', 'N2HP', 'HCOP', 'CS', 'SO', 'NH2D', 'H13COP']

v10 = '_v10_match_cube.fits'
v06 = '_match_cube.fits'
#vtype = v10
vtype_file = v06
vtype = 'v06'

if vtype == 'v10':
    rs = 200
elif vtype == 'v06':
    rs = 300




#####################################################################
##### for GPY+
"'modify a map (replace NaN pixel to a noisy pixel)'"

import funstools
from funstools import get_rms
from funstools import Cube2map
from funstools import save_fits
from matplotlib import pyplot as plt
from funstools import full_line_scan
from funstools import Decompose
#from astroML.plotting import multiaxes
#from multiaxes import Multiaxes
from astropy.wcs import WCS
from run_funstools_Serp import get_filename
from run_funstools_Serp import get_maxmin_maps
from run_funstools_Serp import full_line_scan_vtype

#from SerB_Cube import SerBCube2map as Cube2map
#from SerB_Cube import Cube2map
from SerB_PlotFitRes_FF_faff import SerBDecompose
#from SerB_PlotFitRes import SerBDecompose


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
from astropy.io.ascii import read
import math

##### directory prep

funs_crop_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'
funs_full_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'

#dir = funs_crop_dir
dir = funs_full_dir

region_code_A = 'SerB_A_'
region_code_B = 'SerB_B_'
region_code_main = 'SerB_'

#region_code = region_code_A
#region_code = region_code_B
region_code = region_code_main

wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/rect_cube/'

funsmap_name = [r'$^{13}$CO',\
         r'C$^{18}$O',\
         r'N$_{2}$H$^{+}$', \
         r'HCO$^{+}$', \
         r'CS', \
         r'SO', \
         r'NH$_{2}$D', \
         r'H$^{13}$CO$^{+}$']

funslinename = ['13CO', 'C18O', 'N2HP', 'HCOP', 'CS', 'SO', 'NH2D', 'H13COP']

v10 = '_v10_match_cube.fits'
v06 = '_match_cube.fits'
#vtype = v10
vtype_file = v06
vtype = 'v06'

if vtype == 'v10':
    rs = 200
elif vtype == 'v06':
    rs = 300


decom_line = funslinename[1]

input_file = get_filename(line=decom_line, vtype='v06', ctype='matched')

cube = Cube2map(input_file)
a = cube.data
h=cube.header
h['COMMENT'] = '= modified for GPY+ '

for i in range(cube.nd):
    for j in range(cube.nr):
        if math.isnan(a[0,i,j]):
            a[:,i,j] = a[:,230,3]


wfilename = wdir+'SerB_C18O_v06_rect_cube.fits'
save_fits(wfilename, a, h, overwrite=True)



#######################################
## to make default fitting parameter table

from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_option = 'sn5.0d2s1' # the decomposed paramter option title (i.e., directory name)
dc_appendix = '_sv2511'
vsm= '5.0'
ssm= '2.0'

#gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/rect'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option+dc_appendix,'output', 'gpy_decomposed')
gpy_final = 'SerB_C18O_vsm{}_ssm{}_v06_cube_rect_{}_g+_fit_fin_sf-p2_finalized.dat'.format(vsm, ssm, dc_option)
#gpy2fitres = 'SerB_C18O_vsm{}_ssm{}_v06_cube_rect_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(vsm, ssm, dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

refit_default = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/dc_refit_parameters.dat'

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

final = []
fitcode = []


a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
        final.append('n')
        fitcode.append()
#        print(j)

    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)



xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, amp, VLSR, vel_disp, final], names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'final'))
fr.write(refit_default, format='ascii', overwrite='True')



### To make a table of the dc method for each pixel
from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_option = 'sn5.0d2s1' # the decomposed paramter option title (i.e., directory name)
dc_appendix = '_sv2511'
vsm= '5.0'
ssm= '2.0'

#gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/rect'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option+dc_appendix,'output', 'gpy_decomposed')
gpy_final = 'SerB_C18O_vsm{}_ssm{}_v06_cube_rect_{}_g+_fit_fin_sf-p2_finalized.dat'.format(vsm, ssm, dc_option)
#gpy2fitres = 'SerB_C18O_vsm{}_ssm{}_v06_cube_rect_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(vsm, ssm, dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

refit_default = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/dc_fitting_table_info.dat'

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

final = []
fitcode = []


a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
        final.append('n')
        fitcode.append(1)
#        print(j)

    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)



xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, fitcode, final], names=('rp', 'dp', 'tn', 'cn', 'dc_code', 'final'))
fr.write(refit_default, format='ascii', overwrite='True')





#####
# 중간에 테이블 쓰기 하느라 만들어 본 걸까? _table_info.dat은 이우에 계속 업데이트되었음.
from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

data = read('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/dc_fitting_table_info_org.dat')
wdata = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/dc_refit/dc_fitting_table_info.dat'


rp = []
dp = []
dc_code = []
final = []

ct = data[data['cn']==0]
for i in range(len(ct)):
    rp.append(ct['rp'][i])
    dp.append(ct['dp'][i])
    dc_code.append(ct['dc_code'][i])
    final.append(ct['final'][i])

fr = Table([rp, dp, dc_code, final], names=('rp', 'dp', 'dc_code', 'final'))
fr.write(wdata, format='ascii', overwrite='True')



############
#find peaks besides 5-12km/s (11/22/21)
from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '2_111921')
dc_data_file = 'dc_results_combined_revise_2.dat'
dc_data = dc_dir+'/'+dc_data_file

dc_info = read(dc_data, header_start=0, data_start=1)#, delimiter=' ')
a = iter(range(len(dc_info)))

rev_rp = []
rev_dp = []
rev_tn = []
rev_cn = []
rev_tp = []
rev_vp = []
rev_sd = []
rev_snr = []
rev_dv = []
rev_area = []
rev_dc_code = []
rev_dc_final = []

for i in a:
    if ((dc_info['vp'][i] <= 5) or (dc_info['vp'][i] >=12)) and (dc_info['snr'][i] >= 1) :
        rev_rp.append(dc_info['rp'][i])
        rev_dp.append(dc_info['dp'][i])
        rev_tn.append(dc_info['tn'][i])
        rev_cn.append(dc_info['cn'][i])
        rev_tp.append(dc_info['tp'][i])
        rev_vp.append(dc_info['vp'][i])
        rev_sd.append(dc_info['sd'][i])
        rev_snr.append(dc_info['snr'][i])
        rev_dv.append(dc_info['dv'][i])
        rev_area.append(dc_info['area'][i])
        rev_dc_code.append(dc_info['dc_code'][i])
        rev_dc_final.append(dc_info['dc_final'][i])

dc_outer_bump = dc_dir+'/outer_bumps_check_snr1.dat'
fr = Table([rev_rp, rev_dp, rev_tn, rev_cn, rev_tp, rev_vp, rev_sd, rev_snr, rev_dv, rev_area, rev_dc_code, rev_dc_final], \
names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area', 'dc_code', 'dc_final'))
fr.write(dc_outer_bump, format='ascii', overwrite='True')


for i in a:
    if (dc_info['vp'][i] >= 5) and (dc_info['vp'][i] <= 13) :#and (dc_info['snr'][i] >= 0.6) :
        rev_rp.append(dc_info['rp'][i])
        rev_dp.append(dc_info['dp'][i])
        rev_tn.append(dc_info['tn'][i])
        rev_cn.append(dc_info['cn'][i])
        rev_tp.append(dc_info['tp'][i])
        rev_vp.append(dc_info['vp'][i])
        rev_sd.append(dc_info['sd'][i])
        rev_snr.append(dc_info['snr'][i])
        rev_dv.append(dc_info['dv'][i])
        rev_area.append(dc_info['area'][i])
        rev_dc_code.append(dc_info['dc_code'][i])
        rev_dc_final.append(dc_info['dc_final'][i])

dc_outer_bump = dc_dir+'/w5-13_bumps_check.dat'
fr = Table([rev_rp, rev_dp, rev_tn, rev_cn, rev_tp, rev_vp, rev_sd, rev_snr, rev_dv, rev_area, rev_dc_code, rev_dc_final], \
names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area', 'dc_code', 'dc_final'))
fr.write(dc_outer_bump, format='ascii', overwrite='True')




### check fof algorithm 1

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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 2

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0

step_cut1=sigma_cut
snr=sigma_init

max_vel=0.2
dp0 = 5
max_sep = np.sqrt(dp0**2 + dp0**2)

new_fil = 1
ii = 1

while snr > step_cut1:
    #cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
    if np.sum(cmask) == 0:
        snr -= 0.5
        print('snr:', snr)
        continue
    ii += 1
    icomps = np.where(cmask)[0] # 몇 번째 행에 있는지, 행 index
    for i in icomps:
        di, ri = comps['dp'][i], comps['rp'][i]
        vi = comps['vp'][i]
        seps = np.sqrt((comps['dp']-di)**2+(comps['rp']-ri)**2)
        vdif = np.abs(comps['vp']-vi)
        fmask = (vdif < max_vel) & (seps < max_sep) & (comps['fn'] > 0)
        if np.sum(fmask) == 0:
            comps['fn'][i] = new_fil
            comps['ii'][i] = ii
            new_fil += 1
        else:
            fils = np.unique(comps['fn'][fmask])
            comps['fn'][i] = fils[0]
            comps['ii'][i] = ii
            if len(fils) > 1:
                for fn in fils[1:]:
                    gmask = comps['fn'] == fn
                    comps['fn'][gmask] = fils[0]
    remains = np.sum(comps['fn'] == 0)
    print('snr: {}, ii: {}, new_fil: {}, remains: {}'.format(snr, ii, new_fil, remains))






### check fof algorithm 2

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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 7

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0
comps['in'] = 0

step_cut1=sigma_cut
snr=sigma_init

dist_v0 = 0.06
dist_rd0 = 3


ifil = 1
ii = 0
snr_index = 0
while snr > step_cut1:
    #cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
    ok_index = np.where(cmask==1)
    if np.sum(cmask) == 0:
        snr -= 0.5
        print('snr:', snr)
        snr_index += 1
        continue
    ii += 1

    if snr_index == 0:
        remain_max = np.max(comps['tp'][ok_index])
        kseed = np.argwhere((comps['fn'] == 0) & (comps['tp'] == remain_max))[0]
        ir, id, iv = comps['rp'][kseed], comps['dp'][kseed], comps['vp'][kseed]
        k = 1
        comps['fn'][kseed] = ifil
        comps['in'][kseed] = k
        nfc = 1
        while nfc > 0:
            nfc = 0
            k += 1
            if k == 2 :
                fcond = ((comps['rp'][ok_index] - ir)**2 + (comps['dp'][ok_index] - id)**2 <= dist_rd0**2) \
                        & ~((comps['rp'][ok_index] == ir) & (comps['dp'][ok_index] == id)) & (np.abs(comps['vp'][ok_index]- iv) <= dist_v0)
                nfc += np.sum(fcond)
                comps['fn'][ok_index][fcond] = ifil
                comps['in'][ok_index][fcond] = k
            else :
                seeds = comps[(comps['fn'] == ifil) & (comps['in'] == k-1)]
                for j in range(len(seeds)):
                    ir, id, iv = seeds['rp'][j], seeds['dp'][j], seeds['vp'][j]
                    fcond = (comps['fn'][ok_index] == 0) & ((comps['rp'][ok_index] - ir) ** 2 + (comps['dp'][ok_index] - id) ** 2 <= dist_rd0 ** 2) \
                            & ~((comps['rp'][ok_index] == ir) & (comps['dp'][ok_index] == id)) & (np.abs(comps['vp'][ok_index] - iv) <= dist_v0)
                    nfc += np.sum(fcond)
                    comps['fn'][fcond] = ifil
                    comps['in'][fcond] = k
        ifil += 1

        remains = np.sum(comps['fn'] == 0)
        print('snr: {},  ifil: {}, ii: {}, new_fil: {}, remains: {}'.format(snr, ifil, ii, nfc, remains))

    # else:
    #     nfc = 1
    #     for i in ok_index:
    #         di, ri = comps['dp'][i], comps['rp'][i]
    #         vi = comps['vp'][i]
    #         seps = np.sqrt((comps['dp'][ok_index]-di)**2+(comps['rp'][ok_index]-ri)**2)
    #         vdif = np.abs(comps['vp'][ok_index]-vi)
    #         fmask = (vdif < dist_v0) & (seps < dist_rd0) & (comps['fn'][ok_index] > 0)
    #         if np.sum(fmask) == 0:
    #             comps['fn'][i] = ifil
    #             comps['in'][i] = k
    #             comps['ii'][i] = ii
    #             ifil += 1
    #         else:
    #             k += 1
    #             fils = np.unique(comps['fn'][fmask])
    #             comps['fn'][i] = fils[0]
    #             comps['in'][i] = k
    #             comps['ii'][i] = ii
    #             if len(fils) > 1:
    #                 for fn in fils[1:]:
    #                     gmask = comps['fn'] == fn
    #                     comps['fn'][gmask] = fils[0]
    #
    #     remains = np.sum(comps['fn'] == 0)
    #     print('snr: {},  ifil: {}, ii: {}, new_fil: {}, remains: {}'.format(snr, ifil, ii, nfc, remains))
    #

FoF_save = dc_dir+'/fof/'
write(comps, FoF_save+'FoF_result_'+str(float(sigma_cut))+'sig_'+str(dist_v0)+'_'+str(dist_rd0)+'_test4.dat')#, overwrite=True)







### check fof algorithm 3

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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 7

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0
#comps['in'] = 0

step_cut1=sigma_cut
snr=sigma_init

dist_v0 = 0.06
dist_rd0 = 3


ifil = 1
ii = 0
snr_index = 0
while snr > step_cut1:
    #cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
    ok_index = np.where(cmask==1)
    if np.sum(cmask) == 0:
        snr -= 0.5
        print('snr:', snr)
        snr_index += 1
        continue
    ii += 1
    for i in ok_index:
        di, ri = comps['dp'][i], comps['rp'][i]
        vi = comps['vp'][i]
        seps = np.sqrt((comps['dp'][ok_index]-di)**2+(comps['rp'][ok_index]-ri)**2)
        vdif = np.abs(comps['vp'][ok_index]-vi)
        fmask = (vdif < dist_v0) & (seps < dist_rd0) & (comps['fn'][ok_index] > 0)
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


FoF_save = dc_dir+'/fof/'
write(comps[comps['fn']>0], FoF_save+'FoF_result_'+str(float(sigma_cut))+'sig_'+str(dist_v0)+'_'+str(dist_rd0)+'_test3.dat', overwrite=True)


##############################
### check fof algorithm 4
### FoF (adopting Shinyoung's code)
##############################


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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 2

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0
#comps['in'] = 0

step_cut1=sigma_cut
snr=sigma_init

dist_v0 = 0.03
dist_rd0 = 2


ifil = 1
ii = 0
snr_index = 0
while snr >= step_cut1:
    #cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
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


FoF_save = dc_dir+'/fof3/'
write(comps[comps['fn']>0], FoF_save+'FoF_result_C18O_dc_comb_rev3_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat', overwrite=True)




##############################
### check fof algorithm 5
### FoF (adopting Shinyoung's code) - modify to keep the previously assigned filament inforamtion
##############################


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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 2

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0
#comps['in'] = 0

step_cut1=sigma_cut
snr=sigma_init

dist_v0 = 0.03
dist_rd0 = 3

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

#FoF_save = dc_dir+'/fof1/'
FoF_save = dc_dir+'/fof4/'
write(comps[comps['fn']>0], FoF_save+'FoF_result_C18O_dc_comb_rev3_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat', overwrite=True)


##############################
### check fof algorithm 6
### FoF (adopting Shinyoung's code but use snr istead of area) - modify to keep the previously assigned filament inforamtion
### modify to adjust dist_v0 for the components with large sd (> )
##############################


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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 1.5

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0
#comps['in'] = 0

step_cut1=sigma_cut
snr=sigma_init

sdv = 0.5
dist_v0r = [0.03, 0.1]
dist_rd0 = 2

rms_area = (comps['tp']/comps['snr'])*np.sqrt(2.*np.pi*comps['sd']**2.)

ifil = 1
ii = 0
snr_index = 0
while snr >= step_cut1:
    #cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
    ok_index = np.where(cmask==1)
    if np.sum(cmask) == 0:
        snr -= 0.2
        print('snr:', snr)
        snr_index += 1
        continue
    ii += 1
    icomps = np.where(cmask)[0]
    for i in icomps:
        di, ri, sd = comps['dp'][i], comps['rp'][i], comps['sd'][i]
        vi = comps['vp'][i]
        seps = np.sqrt((comps['dp']-di)**2+(comps['rp']-ri)**2)
        vdif = np.abs(comps['vp']-vi)
        if sd <= sdv:
            dist_v0 = dist_v0r[0]
        else:
            dist_v0 = dist_v0r[1]
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
                    comps['fn'][gmask] = fils[0]  #이부분은 고민의 여지가..

    remains = np.sum(comps['fn'] == 0)
    print('snr: {},  ifil: {}, ii: {}, remains: {}'.format(snr, ifil, ii, remains))


FoF_save = dc_dir+'/fof5/'
write(comps[comps['fn']>0], FoF_save+'FoF_result_C18O_dc_comb_rev3_'+str(float(sigma_cut))+'sig_v'+'0.03-0.10'+'_d'+str(dist_rd0)+'.dat', overwrite=True)




##############################
### check fof algorithm 7
### FoF (adopting Shinyoung's code but use snr istead of area) - modify to keep the previously assigned filament inforamtion
### modify to adjust dist_v0 for the components with large sd (> 0.5km/s) & snr >5 ?
##############################


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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 2

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0
#comps['in'] = 0

step_cut1=sigma_cut
snr=sigma_init

sdv = 0.5
dist_v0r = [0.03, 0.1]
dist_rd0 = 2

rms_area = (comps['tp']/comps['snr'])*np.sqrt(2.*np.pi*comps['sd']**2.)

ifil = 1
ii = 0
snr_index = 0
while snr >= step_cut1:
    #cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
    ok_index = np.where(cmask==1)
    if np.sum(cmask) == 0:
        snr -= 0.2
        print('snr:', snr)
        snr_index += 1
        continue
    ii += 1
    icomps = np.where(cmask)[0]
    for i in icomps:
        di, ri, sd = comps['dp'][i], comps['rp'][i], comps['sd'][i]
        vi = comps['vp'][i]
        seps = np.sqrt((comps['dp']-di)**2+(comps['rp']-ri)**2)
        vdif = np.abs(comps['vp']-vi)
        if sd <= sdv:
            dist_v0 = dist_v0r[0]
        elif (sd > sdv) & (comps['snr'][i] >= 5):
            dist_v0 = dist_v0r[1]
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
                    comps['fn'][gmask] = fils[0]  #이부분은 고민의 여지가..

    remains = np.sum(comps['fn'] == 0)
    print('snr: {},  ifil: {}, ii: {}, remains: {}'.format(snr, ifil, ii, remains))


FoF_save = dc_dir+'/fof8/'
write(comps[comps['fn']>0], FoF_save+'FoF_result_C18O_dc_comb_rev3_'+str(float(sigma_cut))+'sig_v'+'0.03-0.10s5'+'_d'+str(dist_rd0)+'.dat', overwrite=True)



##############################
### check fof algorithm 8
### FoF (adopting Shinyoung's code but use snr istead of area) - modify to keep the previously assigned filament inforamtion
### modify to adjust dist_v0 for the components with large sd (> 0.5km/s) & snr >5 , dist_rd0=3 if snr < 3?
##############################


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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 2

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0
#comps['in'] = 0

step_cut1=sigma_cut
snr=sigma_init

sdv = 0.5
dist_v0r = [0.03, 0.1]
dist_rd0r = [2, 3]

rms_area = (comps['tp']/comps['snr'])*np.sqrt(2.*np.pi*comps['sd']**2.)

ifil = 1
ii = 0
snr_index = 0
while snr >= step_cut1:
    #cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
    ok_index = np.where(cmask==1)
    if np.sum(cmask) == 0:
        snr -= 0.2
        print('snr:', snr)
        snr_index += 1
        continue
    ii += 1
    icomps = np.where(cmask)[0]
    for i in icomps:
        di, ri, sd = comps['dp'][i], comps['rp'][i], comps['sd'][i]
        vi = comps['vp'][i]
        seps = np.sqrt((comps['dp']-di)**2+(comps['rp']-ri)**2)
        vdif = np.abs(comps['vp']-vi)
        if (sd <= sdv):
            dist_v0 = dist_v0r[0]
            dist_rd0 = dist_rd0r[0]

        elif (sd > sdv) & (comps['snr'][i] >= 5):
            dist_v0 = dist_v0r[1]
            dist_rd0 = dist_rd0r[0]

        elif (comps['snr'][i] < 3):
            dist_v0 = dist_v0r[0]
            dist_rd0 = dist_rd0r[1]

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
                    comps['fn'][gmask] = fils[0]  #이부분은 고민의 여지가..

    remains = np.sum(comps['fn'] == 0)
    print('snr: {},  ifil: {}, ii: {}, remains: {}'.format(snr, ifil, ii, remains))


FoF_save = dc_dir+'/fof9/'
write(comps[comps['fn']>0], FoF_save+'FoF_result_C18O_dc_comb_rev3_'+str(float(sigma_cut))+'sig_v'+'0.03-0.10s5'+'_d2-3s3.dat', overwrite=True)





##############################
### check fof algorithm 9
### FoF (adopting Shinyoung's code but use snr istead of area) - modify to keep the previously assigned filament inforamtion
### modify to adjust dist_v0 for the components with large sd (> )
### grouping 에 대한 처리...
##############################


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


dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '3_112221')
dc_data_file = 'dc_results_combined_revise_3.dat'
dc_data = dc_dir+'/'+dc_data_file

sigma_init=10
sigma_cut = 1.5

comps=read(dc_data)
comps['fn'] = 0
comps['ii'] = 0
#comps['in'] = 0
comps['gmask'] = []

step_cut1=sigma_cut
snr=sigma_init

sdv = 0.5
dist_v0r = [0.03, 0.1]
dist_rd0 = 2

rms_area = (comps['tp']/comps['snr'])*np.sqrt(2.*np.pi*comps['sd']**2.)

ifil = 1
ii = 0
snr_index = 0
while snr >= step_cut1:
    #cmask = (comps['area'] > snr*rms_area) & (comps['fn'] == 0)
    cmask = (comps['snr'] > snr) & (comps['fn'] == 0)
    ok_index = np.where(cmask==1)
    if np.sum(cmask) == 0:
        snr -= 0.2
        print('snr:', snr)
        snr_index += 1
        continue
    ii += 1
    icomps = np.where(cmask)[0]
    for i in icomps:
        di, ri, sd = comps['dp'][i], comps['rp'][i], comps['sd'][i]
        vi = comps['vp'][i]
        seps = np.sqrt((comps['dp']-di)**2+(comps['rp']-ri)**2)
        vdif = np.abs(comps['vp']-vi)
        if sd <= sdv:
            dist_v0 = dist_v0r[0]
        else:
            dist_v0 = dist_v0r[1]
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
                    comps['fn'][gmask] = fils[0]  #이부분은 고민의 여지가..
                    comps['gmask'][gmaks] = fils

    remains = np.sum(comps['fn'] == 0)
    print('snr: {},  ifil: {}, ii: {}, remains: {}'.format(snr, ifil, ii, remains))


FoF_save = dc_dir+'/fof10/'
write(comps[comps['fn']>0], FoF_save+'FoF_result_C18O_dc_comb_rev3_'+str(float(sigma_cut))+'sig_v'+'0.03-0.10'+'_d'+str(dist_rd0)+'.dat', overwrite=True)


###
## 12/5

from Plot_FitRes_dc_combined_FF import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
#input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_match_cube.fits'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_v10_match_cube.fits'
#wfilename = region+line+'_vsm'+str(self._velocity_smo)+'_ssm'+str(self._spatial_smo)+'_v06_cube.fits'
region = 'SerB_'
line = 'C18O'
snr=3.0
ssm_find=2
vsm_find=5
ssm_fit=1
vsm_fit=1

#dc_ft_key = 'combined_2'
dc_ft_key = 'w5-13 (rev4)'
aa = PlotFitRes(input_cube, key=dc_ft_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')
#dc_data_file = 'dc_results_combined_revise_2.dat'
dc_data_file = 'dc_results_combined_revise_4_w5-13.dat'
dc_data = dc_dir+'/'+dc_data_file

#plot decomposed line profiles
fitres = read(dc_data)
aa.plot_dc_fit(fitres, vr=(4,12), yr=[-0.5,3.5], n_ch=5, nzp=15)



#####################################
#### Save a cube to a pickle
#### 1/07/22
#####################################

import os
from Plot_FitRes_dc_combined_FF import PlotFitRes

vtype = 'v06'

if vtype == 'v10':
    rs = 200
    #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
elif vtype == 'v06':
    #rs = 300
    rs = 300

dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')

from Plot_FitRes_dc_combined_FF import PlotFitRes
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_match_cube.fits'
#input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_v10_match_cube.fits'
#wfilename = region+line+'_vsm'+str(self._velocity_smo)+'_ssm'+str(self._spatial_smo)+'_v06_cube.fits'
region = 'SerB_'
line = 'C18O'
snr=2.5
ssm_find=2
vsm_find=5
ssm_fit=1
vsm_fit=1
dc_ft_key = 'w5-13 (rev5)'
a = PlotFitRes(input_cube, key=dc_ft_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)
awfile = dc_dir+'/cube_pickle/'+region+line+'_{}_snr{}'.format(vtype,snr)

import pickle
with open(file=awfile+'.pickle', mode='wb') as fa:
    pickle.dump(a, fa)

pickled_file = awfile+'.pickle'
if os.path.exists(pickled_file):
    with open(file=pickled_file, mode='rb') as fr:
        acube = pickle.load(fr)




dc_data_file = 'dc_results_combined_revise_4_w5-13.dat'
dc_data = dc_dir+'/'+dc_data_file

### To check components in each filament and velocity distributions of adjacent filaments
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.5
dist_v0=0.2
dist_rd0=2
nbeam = 2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
# working_dir = dc_dir+'/SeqFF_rev1/'
# fof_wfile = working_dir+'SeqFF_C18O_rev4_010722.dat'
working_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/Filfinder/01252022/SeqFf_rev/'
fof_wfile = working_dir+'SeqFF_C18O_rev5.dat'
acube.plot_fof_profile(fntype='finfn', note='SeqFF stack finfn',line=line, vr=(5,11), yr=[-0.5,3.5], n_ch=9, nzp=15, vtype='v06', nbeam=nbeam, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')


from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.5
dist_v0=0.2
dist_rd0=2
nbeam = 10
snr=2.5

# working_dir = dc_dir+'/SeqFF_rev1/'
# fof_wfile = working_dir+'SeqFF_C18O_rev4_010722.dat'
working_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/Filfinder/01252022/SeqFf_rev/'
fof_wfile = working_dir+'SeqFF_C18O_rev5.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, rmssize=300, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)

b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=nbeam, note='SeqFF stack finfn', fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=nbeam, note='SeqFF stack finfn', fnum_note=False)
b.SeqFF_2D_plot_cores_fnum(overplot_m0=True, nbeam=nbeam, note='SeqFF stack finfn', fnum_note=False)

b.SeqFF_2D_plot_fnum_selfil(overplot_m0=True, nbeam=nbeam, note='SeqFF stack finfn selected', fnum_note=True, selfil=[48, 21, 64, 55, 28, 48, 82 ])

b.SeqFF_3D_plot_selfil(five_fn='finfn', nbeam=nbeam, selfil=[64, 55, 28, 43, 82], fnum_note=True)
b.SeqFF_3D_plot_selfil(five_fn='finfn', nbeam=nbeam, selfil=[64, 55, 28, 43, 82], fnum_note=False)


b.SeqFF_3D_plot(five_fn='finfn', nbeam=nbeam)
#######

import importlib
import Plot_FitRes_dc_combined_FF
Plot_FitRes_dc_combined_FF = importlib.reload(Plot_FitRes_dc_combined_FF)
from Plot_FitRes_dc_combined_FF import PlotFitRes


import importlib
import FUNS_FoF_tools_PltFFRes
FUNS_FoF_tools_PltFFRes = importlib.reload(FUNS_FoF_tools_PltFFRes)
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF


#####################################
#### Combine SeqFF selected table and remained table
#### 1/07/22
#####################################

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io.ascii import read, write
from astropy.table import vstack, Table, sort
from astropy import units as u
import os

dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')
working_dir = dc_dir+'/SeqFF_rev1/'

selected_table = read(working_dir+'SeqFF_selected_C18O_dc_comb_rev4.dat')
left_table = read(working_dir+'SeqFF_remain_C18O_dc_comb_rev4_5.5sig_v0.2_d2_step21.dat')
left_table['step'] = 0
left_table['finfn'] = 0
selected_table['dc_code'] = selected_table['dc_code'].astype(str)
st = vstack([selected_table, left_table])
st.sort(['dp', 'rp', 'tn', 'cn'])
stack_table_name = working_dir+'SeqFF_C18O_rev4_010722.dat'
write(st, stack_table_name, overwrite=True)




#####################################
#### Save a cube (vr=[3,13]) to a pickle
#### 1/25/22
#####################################

import os
from Plot_FitRes_dc_combined_FF import PlotFitRes

vtype = 'v06'

if vtype == 'v10':
    rs = 200
    #filename = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'+self._region+self._line+'_'+self._vtype+'_match_cube.fits'
elif vtype == 'v06':
    #rs = 300
    rs = 50

dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/', 'Filfinder')

from Plot_FitRes_dc_combined_FF import PlotFitRes
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/vrcut_cube/SerB_C18O_0313_match_cube.fits'
#input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_v10_match_cube.fits'
#wfilename = region+line+'_vsm'+str(self._velocity_smo)+'_ssm'+str(self._spatial_smo)+'_v06_cube.fits'
region = 'SerB_'
line = 'C18O'
snr=3.0
ssm_find=2
vsm_find=5
ssm_fit=1
vsm_fit=1
dc_ft_key = 'v0313'
a = PlotFitRes(input_cube, key=dc_ft_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)
awfile = dc_dir+'/cube_pickle/'+region+line+'_{}_{}_snr{}'.format(vtype,dc_ft_key,snr)

import pickle
with open(file=awfile+'.pickle', mode='wb') as fa:
    pickle.dump(a, fa)


pickled_file = awfile+'.pickle'
if os.path.exists(pickled_file):
    with open(file=pickled_file, mode='rb') as fr:
        acube = pickle.load(fr)



#######################
## save dc or ff file to pickle
##
from astropy.io.ascii import read, write
from astropy.table import Table
import pickle

a = read('/Users/khkim/Dropbox/Work5/FUNS-Serp/Filfinder/01252022/SeqFF_rev/SeqFF_C18O_rev5.dat')
with open(file='/Users/khkim/Dropbox/Work5/FUNS-Serp/Filfinder/01252022/SeqFF_rev/SeqFF_C18O_rev5.pickle', mode='wb') as fa:
    pickle.dump(a, fa)



###########################
### test Save_Filaments ###
from Save_filaments import SaveFil

sf = SaveFil()
for i in range(82):
    sf.write2fits_cube(fnum=i+1)
    sf.write2fits_m0(fnum=i+1)
    #sf.write2pickle(fnum=i+1)


#############################
## save N2H+ m0 maps: SerB_A & SerB_B
## snr >=1; snr >= 1.5; snr >= 2
## ssmo = 1, vsmo =1
## ssmo = 2, vsmo = 5

import funstools
from funstools import get_rms
from funstools import Cube2map
from funstools import save_fits
from matplotlib import pyplot as plt
# from funstools import full_line_scan
# from funstools import Decompose
# from astropy.wcs import WCS
# from run_funstools_Serp import get_filename
# from run_funstools_Serp import get_maxmin_maps
# from run_funstools_Serp import full_line_scan_vtype
# from SerB_PlotFitRes_FF import SerBDecompose

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
from astropy.io.ascii import read
#
# from Plot_FitRes_dc_combined_FF import PlotFitRes
# from funstools import save_fits
# import os
# from astropy.io.ascii import read, write

wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/N2HP_m0s/'
region = 'SerB_A_'
line = 'N2HP'
snr=2.0
ssm_find=2
vsm_find=5
ssm_fit=1
vsm_fit=1
rs = 300

cubeA = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/SerB_A_N2HP_match_cube.fits'
cubeB = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/SerB_B_N2HP_match_cube.fits'
m0A_wname = wdir+'SerB_A_N2HP_v06_m0_snr{}_{}{}{}{}.fits'.format(snr, ssm_find, vsm_find, ssm_fit, vsm_fit)
m0B_wname = wdir+'SerB_B_N2HP_v06_m0_snr{}_{}{}{}{}.fits'.format(snr, ssm_find, vsm_find, ssm_fit, vsm_fit)

a = Cube2map(cubeA, rmssize=rs, velocity_smo=vsm_find, spatial_smo=ssm_find, snr=snr)
b = Cube2map(cubeB, rmssize=rs, velocity_smo=vsm_find, spatial_smo=ssm_find, snr=snr)
a.moment0(verbose=True)
b.moment0(verbose=True)
save_fits(m0A_wname, a.m0, a.header2d, overwrite=True)
save_fits(m0B_wname, b.m0, b.header2d, overwrite=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection=a.wcs2d)
ax.imshow(a.m0, origin='lower', vmin=0, vmax=np.nanpercentile(a.m0, 99.9))
fig.savefig(m0A_wname+'.pdf')
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection=b.wcs2d)
ax.imshow(b.m0, origin='lower', vmin=0, vmax=np.nanpercentile(b.m0, 99.9))
fig.savefig(m0B_wname+'.pdf')
plt.close()


#################
## nh2 mask pixel > 2e+22 to reamin faint streamers only
###
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
from astropy.io.ascii import read, write
from astropy.table import Table
from astropy.io import fits
import funstools
from funstools import save_fits

nh2_A = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/HGBS/hr_nH_cutA.fits'
nh2_A_mask='/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/HGBS/hr_nH_cutA_mask.fits'
hdu = fits.open(nh2_A)
data = hdu[0].data
header = hdu[0].header
data2 = np.where(data > 1e+22, np.nan, data)
save_fits(nh2_A_mask, data2, header, overwrite=True)



###########################
#plot fellwalker n2h+ cores on n2h+ m0 map
###########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
import os
from astropy.io.ascii import read
from astropy.wcs import WCS
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
from matplotlib.patches import Rectangle, Circle, Ellipse
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
from astropy.visualization.stretch import SinhStretch, LinearStretch, SqrtStretch, HistEqStretch
from astropy.visualization import ImageNormalize, MinMaxInterval
from astropy.convolution import convolve, Gaussian2DKernel

from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles


clr_list = ['orange', 'crimson', 'green', 'magenta', 'navy']
plt.rcParams["axes.grid"] = False

n2hp_A_m0_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/N2HP_m0s/SerB_A_N2HP_v06_m0_snr1p0_1111.fits'
n2hp_B_m0_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/N2HP_m0s/SerB_B_N2HP_v06_m0_snr1p0_1111.fits'
FW_A_log = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/A_FW19.dat'
FW_B_log = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/B_FW21.dat'
FW_A = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/SerB_A_N2HP_v06_m0_snr1p0_1111_outmap19.fits'
FW_B = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/SerB_B_N2HP_v06_m0_snr1p0_1111_outmap21.fits'
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/'

files = [n2hp_A_m0_file, n2hp_B_m0_file]
logdata = [FW_A_log, FW_B_log]
fwresults = [FW_A, FW_B]
title = ['A', 'B']

for i in range(np.size(files)):
    hdul = fits.open(files[i])
    hdr = hdul[0].header
    hdata = hdul[0].data
    hwcs = WCS(hdr)
    ny = np.shape(hdata)[0]
    nx = np.shape(hdata)[1]

    fw = read(logdata[i])
    nclump = np.size(fw)
    pra = fw['Peak1']*u.deg
    pdec = fw['Peak2']*u.deg
    psky = SkyCoord(pra, pdec, frame='icrs')
    pwidth = fw['Size1']/19.99912 # in pixel unit
    pheight = fw['Size2']/19.99912 # in pixel unit

    fwdul = fits.open(fwresults[i])
    fw_hdr = fwdul[0].header
    fw_data = fwdul[0].data
    fw_wcs = WCS(fw_hdr)

    fig = plt.figure(figsize=(0.1*nx, 0.1*ny))
    ax = fig.add_axes([0.2, 0.07, 0.62,0.9], projection=hwcs)
    im = ax.imshow(hdata, origin='lower', vmin=0, vmax=np.nanmax(hdata)*0.9, cmap='binary')
    #ax.xaxis.set_tick_params(labelsize=12)
    #ax.yaxis.set_tick_params(labelsize=12)
    ax.set_xlabel('RA (J2000)', fontsize=12)
    ax.set_ylabel('Dec (J2000)', fontsize=12)
    ax.set_title(title[i])

    for j in range(nclump):
        fw_array = np.zeros((ny, nx))
        fwarr_y = np.where(fw_data == j+1)[0]
        fwarr_x = np.where(fw_data == j+1)[1]
        fw_array[fwarr_y, fwarr_x] = 1.

        kernel = Gaussian2DKernel(0.4)
        fw_array = convolve(fw_array, kernel)
        #ax.scatter(fwarr_x, fwarr_y, marker='s', s=10, alpha=0.2, color = clr_list[j])
        ax.contour(fw_array, levels=[0.5], alpha=0.6, colors=clr_list[j])

        peak_x, peak_y = fw_wcs.world_to_pixel(psky[j])
        ellipse = Ellipse((peak_x, peak_y), pwidth[j], pheight[j], ec='cyan', fill=False, linewidth=1)
        ax.add_artist(ellipse)


    plt.savefig(wdir+'figs/SerB_{}_N2H+_core_fellwalker.pdf'.format(title[i]))

    # rp = np.arange(nx)
    # dp = np.arange(ny)
    # zz = np.transpose(fw_data)
    ax.contour(fw_data, level=[0.5, 1.5, 2.5, 3.5], colors='r')




###########################
#plot fellwalker n2h+ cores on n2h+ m0 map by reading .FIT output data
###########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
import os
from astropy.io.ascii import read
from astropy.wcs import WCS
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
from matplotlib.patches import Rectangle, Circle, Ellipse
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
from astropy.visualization.stretch import SinhStretch, LinearStretch, SqrtStretch, HistEqStretch
from astropy.visualization import ImageNormalize, MinMaxInterval
from astropy.convolution import convolve, Gaussian2DKernel

from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import matplotlib.patches as patches

clr_list = ['orange', 'crimson', 'green', 'magenta', 'navy']
plt.rcParams["axes.grid"] = False

n2hp_A_m0_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/N2HP_m0s/SerB_A_N2HP_v06_m0_snr1p0_1111.fits'
n2hp_B_m0_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/N2HP_m0s/SerB_B_N2HP_v06_m0_snr1p0_1111.fits'
FW_A_cat = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/SerB_A_N2HP_v06_m0_snr1p0_1111_outcat19.FIT'
FW_B_cat = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/SerB_B_N2HP_v06_m0_snr1p0_1111_outcat21.FIT'
FW_A = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/SerB_A_N2HP_v06_m0_snr1p0_1111_outmap19.fits'
FW_B = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/SerB_B_N2HP_v06_m0_snr1p0_1111_outmap21.fits'
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/findclump/FWres_selected/'

files = [n2hp_A_m0_file, n2hp_B_m0_file]
catdata = [FW_A_cat, FW_B_cat]
fwresults = [FW_A, FW_B]
title = ['A', 'B']

for i in range(np.size(files)):
    hdul = fits.open(files[i])
    hdr = hdul[0].header
    hdata = hdul[0].data
    hwcs = WCS(hdr)
    ny = np.shape(hdata)[0]
    nx = np.shape(hdata)[1]

    # fw = read(logdata[i])
    # nclump = np.size(fw)
    # pra = fw['Peak1']*u.deg
    # pdec = fw['Peak2']*u.deg
    # psky = SkyCoord(pra, pdec, frame='icrs')
    # pwidth = fw['Size1']/19.99912 # in pixel unit
    # pheight = fw['Size2']/19.99912 # in pixel unit

    fit = fits.open(catdata[i])
    nclump = np.shape(fit[1].data['shape'])[0]




    fwdul = fits.open(fwresults[i])
    fw_hdr = fwdul[0].header
    fw_data = fwdul[0].data
    fw_wcs = WCS(fw_hdr)

    fig = plt.figure(figsize=(0.1*nx, 0.1*ny))
    ax = fig.add_axes([0.2, 0.07, 0.62,0.9], projection=hwcs)
    im = ax.imshow(hdata, origin='lower', vmin=0, vmax=np.nanmax(hdata)*0.9, cmap='binary')
    #ax.xaxis.set_tick_params(labelsize=12)
    #ax.yaxis.set_tick_params(labelsize=12)
    ax.set_xlabel('RA (J2000)', fontsize=12)
    ax.set_ylabel('Dec (J2000)', fontsize=12)
    ax.set_title(title[i])

    for j in range(nclump):
        #fw_array = np.zeros((ny, nx))
        #fwarr_y = np.where(fw_data == j+1)[0]
        #fwarr_x = np.where(fw_data == j+1)[1]
        #fw_array[fwarr_y, fwarr_x] = 1.

        #kernel = Gaussian2DKernel(0.4)
        #fw_array = convolve(fw_array, kernel)
        #ax.scatter(fwarr_x, fwarr_y, marker='s', s=10, alpha=0.2, color = clr_list[j])
        #ax.contour(fw_array, levels=[0.5], alpha=0.6, colors=clr_list[j])

        poly0 =  np.array((fit[1].data['shape'][j]).split(' '))[3:].reshape(-1, 2)
        upoly = np.array(poly0, dtype=np.float64)*u.deg
        polysky =SkyCoord(upoly, frame='icrs')
        peak_x, peak_y = fw_wcs.world_to_pixel(polysky)
        ppoint = []
        for k in range(np.size(peak_x)):
            ppoint.append([peak_x[k], peak_y[k]])

        polyshape = patches.Polygon(ppoint, ec='magenta', fill=None, lw=1)

        #Ellipse((peak_x, peak_y), pwidth[j], pheight[j], ec='cyan', fill=False, linewidth=1)
        ax.add_artist(polyshape)


    plt.savefig(wdir+'figs/SerB_{}_N2H+_core_fellwalker_test.pdf'.format(title[i]))
