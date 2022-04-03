
######## use FUNSTOOLs to draw decomposed line profile and FIVE 12/22/21
######## dc data ver4 (4_121921)

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

# ##### directory prep



import importlib
import Plot_FitRes_dc_combined_FF
Plot_FitRes_dc_combined_FF = importlib.reload(Plot_FitRes_dc_combined_FF)
from Plot_FitRes_dc_combined_FF import PlotFitRes


import importlib
import FUNS_FoF_tools_PltFFRes
FUNS_FoF_tools_PltFFRes = importlib.reload(FUNS_FoF_tools_PltFFRes)
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF



############################## 12/22/21 (fof with adjustable dist_v0 for the bright compoents)
############ Plot_FitRes : dc_combined_revise --> FoF
from Plot_FitRes_dc_combined_FF import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_match_cube.fits'
#input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_v10_match_cube.fits'
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
a = PlotFitRes(input_cube, key=dc_ft_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')
#dc_data_file = 'dc_results_combined_revise_2.dat'
dc_data_file = 'dc_results_combined_revise_4_w5-13.dat'
dc_data = dc_dir+'/'+dc_data_file

#plot decomposed line profiles
#fitres = read(dc_data)
#a.plot_dc_fit(fitres, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15)

from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.5
dist_v0=0.06
dist_rd0=2
FoF_save = dc_dir+'/five_test/'
fof_wfile = FoF_save+'FIVE_result_C18O_dc_comb_rev4_step1.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)

b.five_step1_seed_plot(overplot_m0=True, nbeam=2, note='vel_grad <= 1.36km/s/pc, max_dist <= 2 pixels')
b.five_step1_fils_plot(overplot_m0=True, nbeam=2, note='vel_grad <= 1.36km/s/pc, max_dist <= 2 pixels')

fof_wfile2 = FoF_save+'FIVE_result_C18O_dc_comb_rev4_step2.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile2)
c.five_step2_fils_plot(overplot_m0=True, nbeam=5, note='vel_grad <= 1.36km/s/pc, max_dist <= 2 pixels')

#FoF_save = dc_dir+'/five_test_sig6-4_vg1.5_1227code1/'
FoF_save = dc_dir+'/five_test/'
fof_wfile3 = FoF_save+'FIVE_result_C18O_dc_comb_rev4_step3.dat'
d = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile3)
d.five_step3_fils_plot(overplot_m0=True, nbeam=5, note='vel_grad <= 1.36km/s/pc, max_dist <= 2 pixels')

d.five_3D_plot(five_fn='fn3')

a.plot_fof_profile(note='vel_grad:1.5',line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=0, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile2, dcresult=dc_data, dctype='rev4_w5-13')

b.save_fil_by_fil()
b.plot_all_fils_on_chmap()
b.plot_vc_fils_on_chmap()
b.fof_3D_plot()
b.dc_3D_plot()


FoF_save = dc_dir+'/five_test_sig6-4_vg1.5_1227code1/'
#FoF_save = dc_dir+'/five_test/'
fof_wfile3 = FoF_save+'FIVE_result_C18O_dc_comb_rev4_step3.dat'
dd = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile3)
dd.five_step2_fils_plot(overplot_m0=True, nbeam=5, note='vel_grad <= 1.5km/s/pc, max_dist <= 2 pixels')
dd.five_step3_fils_plot(overplot_m0=True, nbeam=0, note='vel_grad <= 1.5km/s/pc, max_dist <= 2 pixels')
dd.five_3D_plot(five_fn='fn3')
dd.five_3D_plot_selfil(five_fn='fn3', selfil=[6, 9])
a.plot_fof_profile(note='sig6-4_vg1.5_1227code1', fntype='fn2', line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=0, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile3, dcresult=dc_data, dctype='rev4_w5-13')



FoF_save = dc_dir+'/five_test_sig6-4_vg1.5_ss0.5_1228_code2/'
#FoF_save = dc_dir+'/five_test/'
fof_wfile3 = FoF_save+'FIVE_result_C18O_dc_comb_rev4_step3_a.dat'
d = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile3)

a.plot_fof_profile(note='sig6-4_vg1.5_ss0.5_1228_code2', fntype='fn3', line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=0, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile3, dcresult=dc_data, dctype='rev4_w5-13')




############## v10
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

FoF_save = dc_dir+'/five_test_sig6-4_vg1.5_ss0.5_1228_code2/'
#FoF_save = dc_dir+'/five_test/'
fof_wfile3 = FoF_save+'FIVE_result_C18O_dc_comb_rev4_step3_a.dat'

aa.plot_fof_profile(note='sig6-4_vg1.5_ss0.5_1228_code2', fntype='fn2', line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v10', nbeam=0, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile3, dcresult=dc_data, dctype='rev4_w5-13')
