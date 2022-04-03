

######## use FUNSTOOLs to draw decomposed line profile and FoF  11/18/21
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



############################## 12/19/21 (fof with adjustable dist_v0 for the bright compoents)
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
fitres = read(dc_data)
a.plot_dc_fit(fitres, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15)

from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.03
dist_rd0=3
FoF_save = dc_dir+'/fof1/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
#a.run_fof()
b.fof_2D_plot(overplot_m0=True, nbeam=2)
b.save_fil_by_fil()
b.plot_all_fils_on_chmap()
b.plot_vc_fils_on_chmap()
b.fof_3D_plot()
b.dc_3D_plot()




from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=6.0
dist_v0=0.03
dist_rd0=2
FoF_save = dc_dir+'/fof3/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
#a.run_fof()
c.fof_2D_plot(overplot_m0=True, nbeam=1.5)
c.save_fil_by_fil()
c.plot_all_fils_on_chmap()
c.plot_vc_fils_on_chmap()




from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.10
dist_rd0=2
FoF_save = dc_dir+'/fof5/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+'0.03-0.10'+'_d'+str(dist_rd0)+'.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
#a.run_fof()
b.fof_2D_plot(overplot_m0=True, nbeam=5)

b.save_fil_by_fil()
b.plot_all_fils_on_chmap()
b.plot_vc_fils_on_chmap()
b.fof_3D_plot()





from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.0
dist_v0=0.1
dist_rd0=2
FoF_save = dc_dir+'/fof8/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev3_'+str(float(sigma))+'sig_v'+'0.03-0.10s5'+'_d'+str(dist_rd0)+'.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
#a.run_fof()
c.fof_2D_plot(overplot_m0=True)
c.save_fil_by_fil()
c.plot_all_fils_on_chmap()
c.plot_vc_fils_on_chmap()


from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.0
dist_v0=0.1
dist_rd0=2
FoF_save = dc_dir+'/fof9/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev3_'+str(float(sigma))+'sig_v'+'0.03-0.10s5'+'_d2-3s3.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
#a.run_fof()
c.plot_chmap()
c.fof_2D_plot(overplot_m0=True)
c.save_fil_by_fil()
c.plot_all_fils_on_chmap()
c.plot_vc_fils_on_chmap()




from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.1
dist_rd0=2
FoF_save = dc_dir+'/fof10/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+'0.06-0.10'+'_d2.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot(overplot_m0=True, nbeam=0, note='fof10: 0.06-1.0')
c.fof_2D_plot(overplot_m0=True, nbeam=2, note='fof10: 0.06-1.0')
c.fof_2D_plot(overplot_m0=True, nbeam=5, note='fof10: 0.06-1.0')
#a.run_fof()
c.plot_chmap()

c.save_fil_by_fil()
c.plot_all_fils_on_chmap()
c.plot_vc_fils_on_chmap()




from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.06
dist_rd0=2
FoF_save = dc_dir+'/fof11/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot(overplot_m0=True, nbeam=0, note='fof11: 0.06')
c.fof_2D_plot(overplot_m0=True, nbeam=2, note='fof11: 0.06')
c.fof_2D_plot(overplot_m0=True, nbeam=5, note='fof11: 0.06')
#a.run_fof()
c.plot_chmap()

c.save_fil_by_fil()
c.plot_all_fils_on_chmap()
c.plot_vc_fils_on_chmap()



from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.0
dist_v0=0.1
dist_rd0=2
FoF_save = dc_dir+'/fof8/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+'0.03-0.10s5'+'_d'+str(dist_rd0)+'.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot(overplot_m0=True, nbeam=0, note='fof8: 0.03-0.1 sd0.5 snr5')
c.fof_2D_plot(overplot_m0=True, nbeam=2, note='fof8: 0.03-0.1 sd0.5 snr5')
c.fof_2D_plot(overplot_m0=True, nbeam=5, note='fof8: 0.03-0.1 sd0.5 snr5')
#a.run_fof()
c.plot_chmap()

c.save_fil_by_fil()
c.plot_all_fils_on_chmap()
c.plot_vc_fils_on_chmap()



from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma_array=[2.5, 3.5]
dist_v0_array=[0.1, 0.12, 0.16, 0.2]
dist_rd0=2
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        FoF_save = dc_dir+'/fof-delv/fof-'+str(dist_v0)+'/'
        fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
        c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        c.fof_2D_plot(overplot_m0=True, nbeam=0, note='fof-dv: {}'.format(dist_v0))
        c.fof_2D_plot(overplot_m0=True, nbeam=2, note='fof-dv: {}'.format(dist_v0))
        c.fof_2D_plot(overplot_m0=True, nbeam=5, note='fof-dv: {}'.format(dist_v0))
        plt.close('all')

#a.run_fof()

c.plot_chmap()

c.save_fil_by_fil()
c.plot_all_fils_on_chmap()
c.plot_vc_fils_on_chmap()


from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma_array=[2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]
dist_v0_array=[0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12, 0.16, 0.2]
dist_rd0=2
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        FoF_save = dc_dir+'/fof-delv/fof-'+str(dist_v0)+'/'
        fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
        c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        c.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='fof-dv: {}; d <= 2'.format(dist_v0))
        c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='fof-dv: {}; d <= 2'.format(dist_v0))
        c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='fof-dv: {}; d <= 2'.format(dist_v0))
        plt.close('all')


from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma_array=[2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]
dist_v0_array=[0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12, 0.16, 0.2]
dist_rd0=2
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        FoF_save = dc_dir+'/fof-delv/fof-'+str(dist_v0)+'/'
        fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
        c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        c.fof_2D_plot(overplot_m0=True, nbeam=1, note='fof-dv: {}; d <= 2'.format(dist_v0))
        c.fof_2D_plot(overplot_m0=True, nbeam=2, note='fof-dv: {}; d <= 2'.format(dist_v0))
        c.fof_2D_plot(overplot_m0=True, nbeam=5, note='fof-dv: {}; d <= 2'.format(dist_v0))
        plt.close('all')



from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma_array=[2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0]
dist_v0_array=[0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12, 0.16, 0.2]
dist_rd0=2
dist_v0=dist_v0_array[6]
sigma = sigma_array[0]
nb=2
FoF_save = dc_dir+'/fof-delv/fof-'+str(dist_v0)+'/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'

a.plot_fof_profile(note='fof-dv:{}, {}sig'.format(dist_v0, sigma),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=nb, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='rev4_w5-13')



#####################################################################
# FOF 12/27/21
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


############################## 12/27/21 (fof with velocity gradient 1.5)
############ Plot_FitRes : dc_combined_revise --> FoF
from Plot_FitRes_dc_combined_FF import PlotFitRes
from funstools import save_fits
import os
from astropy.io.ascii import read, write
wdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/smoothed_cube/'
input_cube = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/SerB_C18O_match_cube.fits'
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
a = PlotFitRes(input_cube, key=dc_ft_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')
#dc_data_file = 'dc_results_combined_revise_2.dat'
dc_data_file = 'dc_results_combined_revise_4_w5-13.dat'
dc_data = dc_dir+'/'+dc_data_file

#plot decomposed line profiles
fitres = read(dc_data)
a.plot_dc_fit(fitres, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15)


from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=7.0
dist_v0=0.03
dist_rd0=2
vg = 1.5
FoF_save = dc_dir+'/fof_vg1.5_d2/'
fof_wfile = FoF_save+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_vg'+str(1.5)+'_d'+str(dist_rd0)+'.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
#a.run_fof()
b.fof_2D_plot(overplot_m0=True, nbeam=2, note='vel_grad <= 1.5km/s/pc, max_dist <= 2 pixels')

a.plot_fof_profile(note='vel_grad:1.5',line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=0, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='funstools')

b.save_fil_by_fil()
b.plot_all_fils_on_chmap()
b.plot_vc_fils_on_chmap()
b.fof_3D_plot()
b.dc_3D_plot()




#####################################################################
# FF Sequence 12/31/21
import funstools
from funstools import get_rms
from funstools import Cube2map
from funstools import save_fits
from matplotlib import pyplot as plt
from funstools import full_line_scan
from funstools import Decompose
from astropy.wcs import WCS
from run_funstools_Serp import get_filename
from run_funstools_Serp import get_maxmin_maps
from run_funstools_Serp import full_line_scan_vtype
from SerB_PlotFitRes_FF import SerBDecompose

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
from astropy.io.ascii import read

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

dc_ft_key = 'w5-13 (rev4)'
a = PlotFitRes(input_cube, key=dc_ft_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')
dc_data_file = 'dc_results_combined_revise_4_w5-13.dat'
dc_data = dc_dir+'/'+dc_data_file


#step 1 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.12
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step1.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#step1 selected
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.12
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))


#after step 1, process checking for step2
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.1
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step1_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st2: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))


a.plot_fof_profile(note='FF check for st2: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')



#step 2 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.1
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step2.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=0, note='ff seq st2 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 2
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.1
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step2.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=0, note='ff seq st2 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 2
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.1
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=0, note='ff seq (selec.) up to step 2'.format(sigma, dist_v0))



#after step 2, process checking for step3
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.05
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step2_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st3: deld < 2, snr>={}, delV < {}'.format(sigma, dist_v0))



#step 3 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.08
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step3.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='ff seq st3 (remain): deld < 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 3
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step3.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='ff seq st3 (selec.): deld < 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 2
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=1, note='ff seq (selec.) up to step 3'.format(sigma, dist_v0))


#after step 3, process checking for step4
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.05
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step3_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st4: deld < 2, snr>={}, delV < {}'.format(sigma, dist_v0))

from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.06
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step3_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st4: deld < 2, snr>={}, delV < {}'.format(sigma, dist_v0))

from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.05
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step3_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st4: deld < 2, snr>={}, delV < {}'.format(sigma, dist_v0))

from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.06
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step3_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st4: deld < 2, snr>={}, delV < {}'.format(sigma, dist_v0))

from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.08
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step3_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st4: deld < 2, snr>={}, delV < {}'.format(sigma, dist_v0))



##########################################


#####################################################################
# FF Sequence 1/3/22 (with FoF maxd <= 2)
import funstools
from funstools import get_rms
from funstools import Cube2map
from funstools import save_fits
from matplotlib import pyplot as plt
from funstools import full_line_scan
from funstools import Decompose
from astropy.wcs import WCS
from run_funstools_Serp import get_filename
from run_funstools_Serp import get_maxmin_maps
from run_funstools_Serp import full_line_scan_vtype
from SerB_PlotFitRes_FF import SerBDecompose

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
from astropy.io.ascii import read

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

dc_ft_key = 'w5-13 (rev4)'
a = PlotFitRes(input_cube, key=dc_ft_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)
dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')
dc_data_file = 'dc_results_combined_revise_4_w5-13.dat'
dc_data = dc_dir+'/'+dc_data_file


#step 1 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.2
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step1.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#step1 selected
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.2
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))


#after step 1, process checking for step2
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
sigma_array=[2.5]
dist_v0_array=[0.16, 0.12, 0.1]
dist_rd0=2
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step1_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st2: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st2: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st2: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))
plt.close('all')

a.plot_fof_profile(note='FF check for st2: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')



#step 2 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.16
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step2.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st2 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 2
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.16
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step2.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st2 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 2
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.16
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step 2'.format(sigma, dist_v0))



#after step 2, process checking for step3
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
sigma_array=[2.5, 3.0]
dist_v0_array=[0.12, 0.1, 0.08, 0.06]
dist_rd0=2
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step2_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st3: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st3: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st3: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))
plt.close('all')


sigma=2.5
dist_v0=0.1
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step2_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st3: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')


#step 3 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.1
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step3.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st3 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 3
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.1
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step3.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st3 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 2
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.1
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step 3'.format(sigma, dist_v0))


#after step 3, process checking for step4
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=3
step2=4
sigma_array=[2.5, 3.0]
dist_v0_array=[0.1, 0.08, 0.06]
dist_rd0=2
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
plt.close('all')

sigma=2.5
dist_v0=0.08
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step3_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st4: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=3, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')



#step 4 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.08
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step4.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st4 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 4
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st4 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 4
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step 4'.format(sigma, dist_v0))


#after step 4, process checking for step5
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=4
step2=5
sigma_array = [2.5, 3.0]
dist_v0_array = [0.08, 0.07, 0.06]
dist_rd0=2
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
plt.close('all')



#step 5 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.07
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step5.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st5 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 5
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.07
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step5.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st5 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 5
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=2.5
dist_v0=0.07
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step 5'.format(sigma, dist_v0))


#after step 5, process checking for step6
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=5
step2=6
sigma_array = [3.0, 3.5]
dist_v0_array = [0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
plt.close('all')


sigma=3.0
dist_v0=0.07
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step5_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st6: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')




#step 6 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.08
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step6.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st6 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 6
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step6.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st6 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 6
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step 6'.format(sigma, dist_v0))



#after step 6, process checking for step7
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=6
step2=7
sigma_array = [3.0, 3.5, 4.0]
dist_v0_array = [0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
plt.close('all')


sigma=3.0
dist_v0=0.06
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step6_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st7: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')




#step 7 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.06
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step7.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st7 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 7
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.06
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step7.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st7 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 7
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.06
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step 7'.format(sigma, dist_v0))


#after step 7, process checking for step8
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=7
step2=8
sigma_array = [3.0, 3.5, 4.0]
dist_v0_array = [0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
plt.close('all')


sigma=3.0
dist_v0=0.08
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step7_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st8: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')




#step 8 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.08
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step8.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st8 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 8
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step8.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st8 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 8
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step 8'.format(sigma, dist_v0))

#after step 8, process checking for step9
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=8
step2=9
sigma_array = [3.0, 3.5, 4.0]
dist_v0_array = [0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
plt.close('all')


sigma=3.0
dist_v0=0.05
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step8_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st9: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')




#step 9 remained
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.05
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step9.dat'
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st9 (remain): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected at Step 9
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.05
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step9.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st9 (selec.): deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0))

#selected up to Step 9
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.0
dist_v0=0.05
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step 9'.format(sigma, dist_v0))


#after step 9, process checking for step 10
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=9
step2=10
sigma_array = [3.5, 4.0]
dist_v0_array = [0.12, 0.1, 0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
plt.close('all')



#step 10 remained
step=10
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.12
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 10
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.12
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 10
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.12
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step))


#after step 10, process checking for step 11
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=10
step2=11
sigma_array = [3.5, 4.0]
dist_v0_array = [0.12, 0.1, 0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
plt.close('all')




#step 11 remained
step=11
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.08
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 10
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 10
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step))


#after step 11, process checking for step 12
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=11
step2=12
sigma_array = [3.5, 4.0, 4.5]
dist_v0_array = [0.12, 0.1, 0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')


sigma=3.5
dist_v0=0.06
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step11_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st11: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')





#step 12 remained
step=12
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.06
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 10
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.06
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 10
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=3.5
dist_v0=0.06
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step))


#after step 11, process checking for step 12
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=12
step2=13
sigma_array = [3.5, 4.0, 4.5]
dist_v0_array = [0.2, 0.16]#0.12, 0.1, 0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')



sigma=4.0
dist_v0=0.12
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step12_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st11: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')





#step 13 remained
step=13
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.12
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 13
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.12
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 13
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.12
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step))


#after step 13, process checking for step 14
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=13
step2=14
sigma_array = [4.0, 4.5, 5.0]
dist_v0_array = [0.2, 0.16, 0.12, 0.1, 0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')


sigma=4.0
dist_v0=0.1
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step13_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st14: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')




#step 14 remained
step=14
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.1
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 14
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.1
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 14
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.1
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step))


#after step 14, process checking for step 15
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=14
step2=15
sigma_array = [4.0]
dist_v0_array = [0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')




#step 15 remained
step=15
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.07
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 15
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.07
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 15
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.0
dist_v0=0.07
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)


#after step 15, process checking for step 16
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=15
step2=16
sigma_array = [4.5, 5.0, 5.5]
dist_v0_array = [0.2, 0.16, 0.12, 0.1, 0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')




#step 16 remained
step=16
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.08
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 16
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 16
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.08
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)

plt.close('all')

#after step 16, process checking for step 17
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=16
step2=17
sigma_array = [3.0, 3.5, 4.5, 5.0]
dist_v0_array = [0.2, 0.12, 0.1]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')




#step 17 remained
step=17
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.07
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 17
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.07
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 17
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.07
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)

plt.close('all')

#after step 17, process checking for step 18
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=17
step2=18
sigma_array = [4.5, 5.0, 5.5]
dist_v0_array = [0.2, 0.12, 0.1, 0.08, 0.07, 0.06, 0.05]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')




#step 18 remained
step=18
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.06
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 18
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.06
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 18
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.06
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)

plt.close('all')

#after step 18, process checking for step 19
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=18
step2=19
sigma_array = [2.5, 3.0, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5]
dist_v0_array = [0.2, 0.12, 0.1, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')



sigma=5.0
dist_v0=0.06
dist_rd0=2
fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step19_remain_{}sig_v{}_d{}.dat'.format(sigma, dist_v0, dist_rd0)
a.plot_fof_profile(note='FF check for st18: deld <= 2, snr>={}, delV < {}'.format(sigma, dist_v0),line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=2, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')



#step 19 remained
step=19
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.05
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 19
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.05
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 19
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=4.5
dist_v0=0.05
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)

plt.close('all')

#after step 18, process checking for step 19
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=19
step2=20
sigma_array = [2.5, 3.0, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5]
dist_v0_array = [0.2, 0.16, 0.12, 0.1, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')


#step 20 remained
step=20
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.0
dist_v0=0.06
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 20
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.0
dist_v0=0.06
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 20
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.0
dist_v0=0.06
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=5, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)

plt.close('all')

#after step 18, process checking for step 19
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=20
step2=21
sigma_array = [2.5, 3.0, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5]
dist_v0_array = [0.2, 0.16, 0.12, 0.1, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')




#step 21 remained
step=21
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.5
dist_v0=0.2
dist_rd0=2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step)
c = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
c.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (remain): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected at Step 21
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.5
dist_v0=0.2
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step)
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq st{} (selec.): deld <= 2, snr>={}, delV < {}'.format(step, sigma, dist_v0))

#selected up to Step 21
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.5
dist_v0=0.2
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=2, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)

plt.close('all')

#after step 21, process checking for step 21
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
from matplotlib import pyplot as plt
step1=21
step2='left over'
sigma_array = [2.5, 3.0, 3.5, 4.5, 5.0]
dist_v0_array = [0.2, 0.16, 0.12, 0.1, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03]
for dist_v0 in dist_v0_array:
    for sigma in sigma_array:
        fof_wfile = dc_dir+'/ff_sequence/'+'process/'+'FoF_result_C18O_SeqFF_step{}_remain_{}sig_v{}_d{}.dat'.format(step1, sigma, dist_v0, dist_rd0)
        b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=1, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=2, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        b.fof_2D_plot_fnum(overplot_m0=True, nbeam=5, note='FF check for st{}: deld <= 2, snr>={}, delV < {}'.format(step2, sigma, dist_v0))
        plt.close('all')



#after Step 21: accumulated

#selected at Step 21
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.5
dist_v0=0.2
dist_rd0=2
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=7, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=7, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=8, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=8, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=9, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=9, note='ff seq (selec.) up to step {}'.format(step), fnum_note=False)

plt.close('all')





#####################################################################
# FF Sequence: accumulated filaments 1/06/21
import funstools
from funstools import get_rms
from funstools import Cube2map
from funstools import save_fits
from matplotlib import pyplot as plt
from funstools import full_line_scan
from funstools import Decompose
from astropy.wcs import WCS
from run_funstools_Serp import get_filename
from run_funstools_Serp import get_maxmin_maps
from run_funstools_Serp import full_line_scan_vtype
from SerB_PlotFitRes_FF import SerBDecompose

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
from astropy.io.ascii import read

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
dc_ft_key = 'w5-13 (rev4)'
a = PlotFitRes(input_cube, key=dc_ft_key, snr=snr, velocity_smo=vsm_find, spatial_smo=ssm_find)

dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')
dc_data_file = 'dc_results_combined_revise_4_w5-13.dat'
dc_data = dc_dir+'/'+dc_data_file

from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
sigma=5.5
dist_v0=0.2
dist_rd0=2
nbeam = 2
FoF_save = dc_dir+'/ff_sequence/'
fof_wfile = FoF_save+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'

working_dir = dc_dir+'/SeqFF_rev1/'
fof_wfile = working_dir+'SeqFF_C18O_rev4_010722.dat'
b = FF(decompose_key=dc_ft_key, decompose_path=dc_dir, decompose_result=dc_data, region=region, line='C18O', vtype='v06', snr=snr, ssm_find=ssm_find, vsm_find=vsm_find, ssm_fit=ssm_fit, vsm_fit=vsm_fit, sigma=sigma, v0=dist_v0, rd0=dist_rd0, fof_wfile=fof_wfile)
b.SeqFF_2D_plot_fnum(overplot_m0=True, nbeam=nbeam, note='ff seq (selec.) up to step {}'.format(step), fnum_note=True)

a.plot_fof_profile(fntype='finfn', note='SeqFF finfn (after step 21): deld <= 2',line=line, vr=(2,14), yr=[-0.5,3.5], n_ch=7, nzp=15, vtype='v06', nbeam=nbeam, dckey=dc_ft_key, filsigmacut=sigma, dv0=dist_v0, drd0=dist_rd0, fof_result=fof_wfile, dcresult=dc_data, dctype='comb_rev4')


b.SeqFF_3D_plot(five_fn='finfn', nbeam=nbeam)
b.SeqFF_3D_plot_selfil(five_fn='finfn', nbeam=nbeam, selfil=[1,2,3], fnum_note='True')

nbeam=2
b.SeqFF_save_fil_by_fil(nbeam=nbeam)
b.SeqFF_plot_vc_fils_on_chmap(nbeam=nbeam)


import importlib
import Plot_FitRes_dc_combined_FF
Plot_FitRes_dc_combined_FF = importlib.reload(Plot_FitRes_dc_combined_FF)
from Plot_FitRes_dc_combined_FF import PlotFitRes


import importlib
import FUNS_FoF_tools_PltFFRes
FUNS_FoF_tools_PltFFRes = importlib.reload(FUNS_FoF_tools_PltFFRes)
from FUNS_FoF_tools_PltFFRes import FUNS_FoF as FF
