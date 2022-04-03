
############ GaussPyPlus  #########################
#### run FoF with decomposed results by GPY+

from FUNS_FoF_tools import FUNS_FoF as FF
dckey = 'sn3.0d2s1'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'+dckey+'/output/gpy_decomposed/'
dcdata = dcpath+'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(dckey)
region='SerB_A_'
# snr=3.0
# ssm_find=2
# vsm_find=2
# ssm_fit=1
# vsm_fit=1
# sigma=3.0
# v0=0.02
# rd0=2.0
snr=3.0
ssm_find=2
vsm_find=4
ssm_fit=1
vsm_fit=1
sigma=3.0
v0=0.06
rd0=2
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=v0, rd0=rd0, fof_wfile=fof_wfile)
a.run_fof()
a.fof_2D_plot()


import importlib
import FUNS_FoF_tools
FUNS_FoF_tools = importlib.reload(FUNS_FoF_tools)
from FUNS_FoF_tools import FUNS_FoF as FF




############ Plot_FitRes : GPY+ --> FoF
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/SerB_A_C18O_match_cube.fits'
#wfilename = region+line+'_vsm'+str(self._velocity_smo)+'_ssm'+str(self._spatial_smo)+'_v06_cube.fits'
region = 'SerB_A_'
line = 'C18O'
snr=3.0
ssm_find=3
vsm_find=2
sigma=3.0
v0=0.03
rd0=3.0
dc_gpy_key = 'sn3.0d2s1'
a = PlotFitRes(input_cube, key=dc_gpy_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)

# #smoothed file
# b = a.save_smoothed_3dcube()
# h=a.header
# h['COMMENT'] = '= spatial smoothing by '+str(a._spatial_smo)
# h['COMMENT'] = '= veolocity smoothing by '+str(a._velocity_smo)
# wfilename = 'SerB_A_C18O_vsm'+str(a._velocity_smo)+'_ssm'+str(a._spatial_smo)+'_v06_cube.fits'
# save_fits(wdir+wfilename, b, h, overwrite=True)

#gpy+_dc_result file
dc_gpy_key = dc_gpy_key
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/', dc_gpy_key, 'output', \
'gpy_decomposed')
dc_fin_filename = region+line+'_match_cube_'+dc_gpy_key+'_g+_fit_fin_sf-p2_finalized_fitres.dat'
dc_data_file = dc_dir+'/'+dc_fin_filename

fitres = read(dc_data_file)
a.plot_dc_fit(fitres, vr=(3,13), yr=[-0.5,3.5], n_ch=5, nzp=15)


fof_file = dc_dir+'/fof/FoF_result_'+region+line+'_match_cube_'+dc_gpy_key+'_g+sig_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'
dcresult = dc_data_file
a.plot_fof_profile(line=line, vr=(3,13), yr=[-0.5,3.5], n_ch=5, nzp=15, vtype='v06', nbeam=5, dckey=dc_gpy_key, filsigmacut=sigma, dv0=v0, drd0=rd0, fof_result=fof_file, dcresult=dcresult, dctype='gpy+')


import importlib
import Plot_FitRes_GPYFUNSTOOLs
Plot_FitRes_GPYFUNSTOOLs = importlib.reload(Plot_FitRes_GPYFUNSTOOLs)
from Plot_FitRes_GPYFUNSTOOLs import PlotFitRes





############ FUNSTOOLs  #########################
#### run FoF with decomposed results by funstools

from FUNS_FoF_tools import FUNS_FoF as FF
dckey = '2421v06snr3mww0.060.064.0'
dcpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose/decompose6_AB/A/'+dckey+'/'
dcdata = dcpath+'SerB_A_C18O_fitres3_dc3.dat'
region='SerB_A_'
snr=3.0
ssm_find=2
vsm_find=4
ssm_fit=2
vsm_fit=1
sigma=3.0
v0=0.06
rd0=2.0
fof_wfile = dcpath+'fof/'+'FoF_result_'+region+'C18O_match_cube_'+dckey+'_ft_'+str(sigma)+'_dv0_'+str(v0)+'_rd0_'+str(rd0)+'.dat'

a = FF(decompose_key=dckey, decompose_path=dcpath, decompose_result=dcdata, region='SerB_A_', line='C18O', vtype='v06', \
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
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/SerB_A_C18O_match_cube.fits'
#wfilename = region+line+'_vsm'+str(self._velocity_smo)+'_ssm'+str(self._spatial_smo)+'_v06_cube.fits'
region = 'SerB_'
line = 'C18O'
sigma=3.0
v0=0.02
rd0=3.0
dc_ft_key = '2421v06snr3mww0.060.064.0'
a = PlotFitRes(input_cube, key=dc_ft_key)

# #smoothed file
# b = a.save_smoothed_3dcube()
# h=a.header
# h['COMMENT'] = '= spatial smoothing by '+str(a._spatial_smo)
# h['COMMENT'] = '= veolocity smoothing by '+str(a._velocity_smo)
# wfilename = 'SerB_A_C18O_vsm'+str(a._velocity_smo)+'_ssm'+str(a._spatial_smo)+'_v06_cube.fits'
# save_fits(wdir+wfilename, b, h, overwrite=True)

#gpy+_dc_result file
dc_ft_key = '2421v06snr3mww0.060.064.0'
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'decompose6_AB', 'A', dc_ft_key)
dc_fin_filename = region+line+'_fitres3.dat'
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
